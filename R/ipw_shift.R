#' IPW Estimates of the Causal Effects of Stochatic Shift Interventions
#'
#' @param W A \code{matrix}, \code{data.frame}, or similar containing a set of
#'  baseline covariates.
#' @param A A \code{numeric} vector corresponding to a exposure variable. The
#'  parameter of interest is defined as a location shift of this quantity.
#' @param Y A \code{numeric} vector of the observed outcomes.
#' @param delta A \code{numeric} value indicating the shift in the exposure to
#'  be used in defining the target parameter. This is defined with respect to
#'  the scale of the exposure (A).
#' @param n_bins A \code{numeric}, scalar or vector, indicating the number of
#'  bins into which the support of A is to be partitioned for constructing
#'  conditional density estimates.
#' @param cv_folds A \code{numeric} giving the number of folds to be used for
#'  cross-validation. Note that this form of sample splitting is used for the
#'  selection of tuning parameters by empirical risk minimization, not for the
#'  estimation of nuisance parameters (i.e., to relax regularity conditions).
#' @param lambda_seq A \code{numeric} sequence of the regularization parameter
#'  (L1 norm of HAL coefficients) to be used in fitting HAL models.
#' @param ... Additional arguments for model fitting to be passed directly to
#'  \code{\link[haldensify]{haldensify}}.
#' @param bin_type A \code{character} indicating the strategy to be used in
#'  creating bins along the observed support of \code{A}. For bins of equal
#'  range, use \code{"equal_range"}; to ensure each bin has the same number of
#'  observations, use instead \code{"equal_mass"}. For more information, see
#'  documentation of \code{grid_type} in \code{\link[haldensify]{haldensify}}.
#' @param trim_density A \code{logical} indicating whether estimates of the
#'  conditional density should be trimmed. Refer to the \code{trim} argument of
#'  the \code{predict} method of \code{haldensify} for details. The default is
#'  \code{FALSE} since propensity score truncation can lead to estimation bias.
#' @param undersmooth_type A \code{character} indicating the selection strategy
#'  to be used in identifying an efficent IPW estimator. The choices include
#'  \code{"gcv"} for global cross-validation, \code{"dcar"} for solving the
#'  IPW representation of the EIF through, and \code{"plateau"} for an approach
#'  that balances changes in the parameter estimate and its mean squared error,
#'  based on Lepski's method. The option \code{"all"} produces results based on
#'  all three selection strategies, sharing redundant computation between each.
#' @param bootstrap A \code{logical} indicating whether the estimator variance
#'  should be approximated using the nonparametric bootstrap. The default is
#'  \code{FALSE}, in which case the empirical variances of the IPW estimating
#'  function and the EIF are used for for estimator selection and for variance
#'  estimation, respectively. When set to \code{TRUE}, the bootstrap variance
#'  is used for both of these purposes instead. Note that the bootstrap is very
#'  computationally intensive and scales relatively poorly.
#' @param n_boot A \code{numeric} giving the number of bootstrap re-samples to
#'  be used in computing the mean squared error as part of the plateau selector
#'  criterion. Ignored when \code{undersmooth_type} is not \code{"plateau"}.
#'
#' @importFrom hal9001 fit_hal
#' @importFrom stats predict
#'
#' @export
#'
#' @examples
#' # simulate data
#' n_obs <- 50
#' W1 <- rbinom(n_obs, 1, 0.6)
#' W2 <- rbinom(n_obs, 1, 0.2)
#' A <- rnorm(n_obs, (2 * W1 - W2 - W1 * W2), 2)
#' Y <- rbinom(n_obs, 1, plogis(3 * A + W1 + W2 - W1 * W2))
#'
#' # fit the IPW estimator
#' est_ipw_shift <- ipw_shift(
#'   W = cbind(W1, W2), A = A, Y = Y,
#'   delta = 0.5, n_bins = 3L, cv_folds = 2L,
#'   lambda_seq = exp(seq(-1, -10, length = 100L)),
#'   # arguments passed to hal9001::fit_hal()
#'   max_degree = 1,
#'   # ...continue arguments for IPW
#'   undersmooth_type = "gcv"
#' )
ipw_shift <- function(W, A, Y,
                      delta,
                      n_bins = make_bins(A, "hist"),
                      cv_folds = 10L,
                      lambda_seq,
                      ...,
                      bin_type = c("equal_range", "equal_mass"),
                      trim_density = FALSE,
                      undersmooth_type = c("dcar", "plateau", "gcv", "all"),
                      bootstrap = FALSE,
                      n_boot = 1000L) {
  # catch dot args for convenience
  dot_args <- list(...)

  # outcome family for selectors that require outcome regression Qn
  outcome_levels <- length(unique(Y))
  outcome_family <- ifelse(outcome_levels > 2L, "gaussian", "binomial")
  n_obs <- length(Y)

  # fit haldensify for generalized propensity score
  gn_fit_haldensify <- haldensify(
    A = A, W = W,
    n_bins = n_bins,
    cv_folds = cv_folds,
    lambda_seq = lambda_seq,
    grid_type = bin_type,
    ## arguments passed to hal9001::fit_hal()
    ...
  )
  cv_lambda_idx <- gn_fit_haldensify$cv_tuning_results$lambda_loss_min_idx
  cv_n_bins <- gn_fit_haldensify$n_bins_cvselect

  # generalized propensity score predictions for...
  ## 1) "natural" A
  gn_pred_natural <- stats::predict(
    gn_fit_haldensify,
    new_A = A,
    new_W = W,
    trim = trim_density,
    lambda_select = "undersmooth"
  )

  ## 2) counterfactual (down)shifted A
  gn_pred_shifted <- stats::predict(
    gn_fit_haldensify,
    new_A = (A - delta),
    new_W = W,
    trim = trim_density,
    lambda_select = "undersmooth"
  )

  # fit outcome mechanism Qn via CV-HAL
  Qn_fit <- hal9001::fit_hal(
    X = cbind(A, W),
    Y = Y,
    family = outcome_family,
    fit_control = list(
      cv_select = TRUE,
      n_folds = cv_folds
    )
  )

  # outcome predictions for "natural" A and counterfactual (up)shifted A
  Qn_pred_natural <- stats::predict(Qn_fit, new_data = cbind(A, W))
  Qn_pred_shifted <- stats::predict(Qn_fit, new_data = cbind(A + delta, W))

  if (undersmooth_type == "dcar") {
    # select an IPW estimate based on ~minimization of D_CAR projection
    ipw_est <- dcar_selector(
      W = W, A = A, Y = Y, delta = delta,
      gn_pred_natural = gn_pred_natural,
      gn_pred_shifted = gn_pred_shifted,
      Qn_pred_natural = Qn_pred_natural,
      Qn_pred_shifted = Qn_pred_shifted
    )
  } else if (undersmooth_type == "plateau") {
    # select an IPW estimate based on changes in psi and MSE
    ipw_est <- plateau_selector(W, A, Y, delta,
      gn_pred_natural = gn_pred_natural,
      gn_pred_shifted = gn_pred_shifted,
      gn_fit_haldensify = gn_fit_haldensify,
      Qn_pred_natural = Qn_pred_natural,
      Qn_pred_shifted = Qn_pred_shifted,
      cv_folds = cv_folds,
      bootstrap = bootstrap,
      n_boot = n_boot,
      ## arguments passed to hal9001::fit_hal()
      ...
    )
  } else if (undersmooth_type == "gcv") {
    # select an IPW estimate based on global cross-validation
    ipw_est <- gcv_selector(
      W = W, A = A, Y = Y, delta = delta,
      gn_pred_natural = gn_pred_natural,
      gn_pred_shifted = gn_pred_shifted,
      Qn_pred_natural = Qn_pred_natural,
      Qn_pred_shifted = Qn_pred_shifted
    )
  } else if (undersmooth_type == "all") {
    # select an IPW estimate based on global cross-validation
    ipw_est_gcv <- gcv_selector(
      W = W, A = A, Y = Y, delta = delta,
      gn_pred_natural = gn_pred_natural,
      gn_pred_shifted = gn_pred_shifted,
      Qn_pred_natural = Qn_pred_natural,
      Qn_pred_shifted = Qn_pred_shifted
    )

    # select an IPW estimate based on ~minimization of D_CAR projection
    ipw_est_dcar <- dcar_selector(
      W = W, A = A, Y = Y, delta = delta,
      gn_pred_natural = gn_pred_natural,
      gn_pred_shifted = gn_pred_shifted,
      Qn_pred_natural = Qn_pred_natural,
      Qn_pred_shifted = Qn_pred_shifted
    )

    # select an IPW estimate based on changes in psi and MSE
    ipw_est_plateau <- plateau_selector(W, A, Y, delta,
      gn_pred_natural = gn_pred_natural,
      gn_pred_shifted = gn_pred_shifted,
      gn_fit_haldensify = gn_fit_haldensify,
      Qn_pred_natural = Qn_pred_natural,
      Qn_pred_shifted = Qn_pred_shifted,
      cv_folds = cv_folds,
      bootstrap = bootstrap,
      n_boot = n_boot,
      ## arguments passed to hal9001::fit_hal()
      ...
    )

    # organize output across different selectors
    est_mat <- rbind(ipw_est_gcv$est, ipw_est_dcar$est, ipw_est_plateau$est)
    eif_mat <- cbind(
      gcv = ipw_est_gcv$eif,
      ipw_est_dcar$eif,
      ipw_est_plateau$eif
    )
    ipw_est <- list(est = est_mat, eif = eif_mat)
  }

  # set reported indices of selected lambdas to match input sequence index
  ipw_est$est$lambda_idx <- ipw_est$est$lambda_idx + cv_lambda_idx - 1
  ipw_est$est$gn_nbins <- cv_n_bins

  # add hidden slots with useful diagnostics
  ipw_est$.delta <- delta
  ipw_est$.outcome_levels <- outcome_levels
  class(ipw_est) <- "ipw_haldensify"
  return(ipw_est)
}
