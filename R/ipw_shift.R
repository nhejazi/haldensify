utils::globalVariables(c("lambda_idx", "se_est", "l1_norm", "type"))

#' IPW Estimator of the Causal Effects of Additive Modified Treatment Policies
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
#'  \code{\link{haldensify}}.
#' @param bin_type A \code{character} indicating the strategy to be used in
#'  creating bins along the observed support of \code{A}. For bins of equal
#'  range, use \code{"equal_range"}; to ensure each bin has the same number of
#'  observations, use instead \code{"equal_mass"}. For more information, see
#'  documentation of \code{grid_type} in \code{\link{haldensify}}.
#' @param selector_type A \code{character} indicating the selection strategy
#'  for identifying an efficent IPW estimator. The choices include \code{"gcv"}
#'  for global cross-validation, \code{"dcar"} for solving the EIF equation,
#'  and \code{"plateau"} for agnostic approaches (1) balancing changes in the
#'  IPW estimate and its standard error (adapting Lepski's method) and (2) a
#'  plateau detector for inflection points in the IPW estimator's trajectory.
#'  The option \code{"all"} runs all three selection strategies while sharing
#'  redundant computation between each.
#'
#' @importFrom hal9001 fit_hal
#' @importFrom matrixStats colSums2 colMeans2
#' @importFrom stats predict weighted.mean
#' @importFrom data.table as.data.table data.table setnames
#' @importFrom tibble add_column as_tibble tibble
#' @importFrom dplyr mutate relocate "%>%"
#' @importFrom stringr str_remove
#'
#' @export
#'
#' @examples
#' # simulate data
#' set.seed(11249)
#' n_obs <- 50
#' W1 <- rbinom(n_obs, 1, 0.6)
#' W2 <- rbinom(n_obs, 1, 0.2)
#' W3 <- rpois(n_obs, 3)
#' A <- rpois(n_obs, 3 * W1 - W2 + 2 * W1 * W2 + 4)
#' Y <- rbinom(n_obs, 1, plogis(A + W1 + W2 - W3 - W1 * W3))
#'
#' # fit the IPW estimator
#' est_ipw <- ipw_shift(
#'   W = cbind(W1, W2, W3), A = A, Y = Y,
#'   delta = 0.5, cv_folds = 3L,
#'   n_bins = 4L, bin_type = "equal_range",
#'   lambda_seq = exp(seq(-1, -10, length = 100L)),
#'   # arguments passed to hal9001::fit_hal()
#'   max_degree = 1L,
#'   smoothness_orders = 0,
#'   reduce_basis = 1 / sqrt(n_obs)
#' )
ipw_shift <- function(W, A, Y,
                      delta = 0,
                      n_bins = make_bins(A, "hist"),
                      cv_folds = 10L,
                      lambda_seq,
                      ...,
                      bin_type = c("equal_range", "equal_mass"),
                      selector_type = c("dcar", "plateau", "gcv", "all")) {
  # default: return IPW estimators from all three selectors
  selector_type <- match.arg(selector_type)

  # outcome family for selectors that require outcome regression Qn
  outcome_family <- ifelse(length(unique(Y)) > 2L, "gaussian", "binomial")
  outcome_levels <- length(unique(Y))
  n_obs <- length(Y)

  # set truncation level for density predictions
  gn_trunc <- 1 / sqrt(n_obs) / log(n_obs)

  # fit haldensify for generalized propensity score
  gn_fit_haldensify <- haldensify(
    A = A, W = W,
    n_bins = n_bins,
    cv_folds = cv_folds,
    lambda_seq = lambda_seq,
    grid_type = bin_type,
    # the following arguments are passed to hal9001::fit_hal()
    ...
  )

  # get L1 norms and lasso lambdas for all HAL fits + CV-selected bin number
  l1_norm_grid <- matrixStats::colSums2(abs(gn_fit_haldensify$hal_fit$coefs))
  lambda_grid <- gn_fit_haldensify$cv_tuning_results$lambda_seq
  cv_lambda_idx <- gn_fit_haldensify$cv_tuning_results$lambda_loss_min_idx
  cv_nbins <- gn_fit_haldensify$n_bins_cvselect

  # generalized propensity score predictions for "natural" A
  # NOTE: this is the denominator in the IP weights and can be truncated
  gn_pred_natural <- stats::predict(
    gn_fit_haldensify,
    new_A = A,
    new_W = W,
    lambda_select = "undersmooth"
  )
  gn_pred_natural <- apply(gn_pred_natural, 2, pmax, gn_trunc)

  # generalized propensity score predictions for counterfactual (down)shifted A
  # NOTE: this is the numerator in the IP weights and should _NOT_ be truncated
  gn_pred_downshift <- stats::predict(
    gn_fit_haldensify,
    new_A = (A - delta),
    new_W = W,
    lambda_select = "undersmooth"
  )

  # fit outcome mechanism Qn via CV-HAL
  Qn_fit <- hal9001::fit_hal(
    X = cbind(A, W), Y = Y,
    max_degree = 3L,
    smoothness_orders = 0L,
    reduce_basis = 1 / sqrt(n_obs),
    family = outcome_family,
    fit_control = list(
      cv_select = TRUE,
      nfolds = cv_folds
    ),
    yolo = FALSE
  )

  # outcome predictions for "natural" and counterfactual (shifted) A
  Qn_pred_natural <- stats::predict(
    Qn_fit,
    new_data = cbind(A, W)
  )
  Qn_pred_shifted <- stats::predict(
    Qn_fit,
    new_data = cbind(A + delta, W)
  )

  if (selector_type == "dcar") {
    # select IPW estimate based on ~minimization of D_CAR projection
    est_out <- dcar_selector(
      W = W, A = A, Y = Y, delta = delta,
      gn_pred_natural = gn_pred_natural,
      gn_pred_shifted = gn_pred_downshift,
      Qn_pred_natural = Qn_pred_natural,
      Qn_pred_shifted = Qn_pred_shifted
    )
  } else if (selector_type == "plateau") {
    # select IPW estimate based on changes in psi and MSE
    est_out <- plateau_selector(W, A, Y, delta,
      gn_pred_natural = gn_pred_natural,
      gn_pred_shifted = gn_pred_downshift,
      gn_fit_haldensify = gn_fit_haldensify,
      Qn_pred_natural = Qn_pred_natural,
      Qn_pred_shifted = Qn_pred_shifted,
      cv_folds = cv_folds
    )
  } else if (selector_type == "gcv") {
    # select IPW estimate based on global cross-validation
    est_out <- gcv_selector(
      W = W, A = A, Y = Y, delta = delta,
      gn_pred_natural = gn_pred_natural,
      gn_pred_shifted = gn_pred_downshift,
      Qn_pred_natural = Qn_pred_natural,
      Qn_pred_shifted = Qn_pred_shifted
    )
  } else if (selector_type == "all") {
    # select IPW estimate based on global cross-validation
    ipw_est_gcv <- gcv_selector(
      W = W, A = A, Y = Y, delta = delta,
      gn_pred_natural = gn_pred_natural,
      gn_pred_shifted = gn_pred_downshift,
      Qn_pred_natural = Qn_pred_natural,
      Qn_pred_shifted = Qn_pred_shifted
    )

    # select IPW estimate based on ~minimization of D_CAR projection
    ipw_est_dcar <- dcar_selector(
      W = W, A = A, Y = Y, delta = delta,
      gn_pred_natural = gn_pred_natural,
      gn_pred_shifted = gn_pred_downshift,
      Qn_pred_natural = Qn_pred_natural,
      Qn_pred_shifted = Qn_pred_shifted
    )

    # select IPW estimate based on changes in psi and MSE
    ipw_est_plateau <- plateau_selector(W, A, Y, delta,
      gn_pred_natural = gn_pred_natural,
      gn_pred_shifted = gn_pred_downshift,
      gn_fit_haldensify = gn_fit_haldensify,
      Qn_pred_natural = Qn_pred_natural,
      Qn_pred_shifted = Qn_pred_shifted,
      cv_folds = cv_folds
    )

    # organize output across different selectors
    est_mat <- rbind(
      ipw_est_gcv$est, ipw_est_dcar$est, ipw_est_plateau$est
    )
    eif_mat <- cbind(
      gcv = ipw_est_gcv$eif,
      ipw_est_dcar$eif,
      ipw_est_plateau$eif
    )
    est_out <- list(est = est_mat, eif = eif_mat)
  }

  # set reported indices of selected lambdas to match input sequence index
  # and add other useful metrics (e.g., L1 norms, histogram binning)
  est_out$est <- est_out$est %>%
    dplyr::mutate(
      lambda_idx = lambda_idx + cv_lambda_idx - 1,
      l1_norm = l1_norm_grid[lambda_idx],
      gn_nbins = cv_nbins
    ) %>%
    dplyr::relocate(type, .after = se_est) %>%
    dplyr::relocate(l1_norm, .before = lambda_idx)

  # add hidden slots with useful diagnostics
  est_out$.delta <- delta
  est_out$.outcome_levels <- outcome_levels
  class(est_out) <- "ipw_haldensify"
  return(est_out)
}
