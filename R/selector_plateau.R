#' Agnostic IPW Estimator Selector via Lepski's and Variance-Blind Methods
#'
#' @param W A \code{matrix}, \code{data.frame}, or similar containing a set of
#'  baseline covariates.
#' @param A A \code{numeric} vector corresponding to a exposure variable. The
#'  parameter of interest is defined as a location shift of this quantity.
#' @param Y A \code{numeric} vector of the observed outcomes.
#' @param delta A \code{numeric} value indicating the shift in the exposure to
#'  be used in defining the target parameter. This is defined with respect to
#'  the scale of the exposure (A).
#' @param gn_pred_natural A \code{matrix} of conditional density estimates of
#'  the exposure mechanism g(A|W) along a grid of the regularization parameter,
#'  at the natural (observed, actual) values of the exposure.
#' @param gn_pred_shifted A \code{matrix} of conditional density estimates of
#'  the exposure mechanism g(A+delta|W) along a grid of the regularization
#'  parameter, at the shifted (counterfactual) values of the exposure.
#' @param gn_fit_haldensify An object of class \code{haldensify} of the fitted
#'  conditional density model for the natural exposure mechanism. This should
#'  be the fit object returned by \code{\link{haldensify}[haldensify]} as part
#'  of a call to \code{\link{ipw_shift}}.
#' @param Qn_pred_natural A \code{numeric} of the outcome mechanism estimate at
#'  the natural (i.e., observed) values of the exposure. HAL regression is used
#'  for the estimate, with the regularization term chosen by cross-validation.
#' @param Qn_pred_shifted A \code{numeric} of the outcome mechanism estimate at
#'  the shifted (i.e., counterfactual) values of the exposure. HAL regression
#'  is used for the estimate, with the regularization term chosen by
#'  cross-validation.
#' @param cv_folds A \code{numeric} giving the number of folds to be used for
#'  cross-validation. Note that this form of sample splitting is used for the
#'  selection of tuning parameters by empirical risk minimization, not for the
#'  estimation of nuisance parameters (i.e., to relax regularity conditions).
#' @param ci_level A \code{numeric} indicating the confidence level to be used
#'  in determining the cutoff used by the Lepski-type selector. This is only
#'  exposed for the sake of accommodating experimentation.
#' @param l1norm_mult A \code{numeric} indicating the multipler to be used by
#'  the plateau-based selector in reducing the candidate set of L1 norms
#'  relative to the choice made by the cross-validation selector.
#'
#' @importFrom matrixStats colVars colMeans2 colSums2 diff2
#' @importFrom stats predict qnorm weighted.mean loess
#' @importFrom tibble tibble as_tibble
#' @importFrom dplyr between
#'
#' @keywords internal
plateau_selector <- function(W, A, Y,
                             delta = 0,
                             gn_pred_natural,
                             gn_pred_shifted,
                             gn_fit_haldensify,
                             Qn_pred_natural,
                             Qn_pred_shifted,
                             cv_folds = 10L,
                             ci_level = 0.95,
                             l1norm_mult = 10L) {
  #   # constants and IP weights
  n_obs <- length(Y)
  ci_mult <- abs(stats::qnorm(p = (1 - ci_level) / 2))

  # shorten sequence of lambdas to just the undersmoothed values
  cv_lambda <- gn_fit_haldensify$cv_tuning_results$lambda_loss_min
  lambda_seq <- gn_fit_haldensify$cv_tuning_results$lambda_seq
  lambda_usm_seq <- lambda_seq[lambda_seq <= cv_lambda]
  nlambda_usm <- length(lambda_usm_seq)
  l1norm_grid <- matrixStats::colSums2(abs(gn_fit_haldensify$hal_fit$coefs))
  l1norm_usm_grid <- l1norm_grid[lambda_seq <= cv_lambda]

  # compute IPW estimators for IP weights along the regularization sequence
  ip_wts_mat <- gn_pred_shifted / gn_pred_natural
  ipw_est <- apply(ip_wts_mat, 2, function(ipw) {
    psi_ipw <- stats::weighted.mean(Y, ipw)
    dipw_score <- ipw * (Y - psi_ipw)
    return(list(psi = psi_ipw, dipw = dipw_score))
  })
  psi_ipw_lambda <- do.call(c, lapply(ipw_est, `[[`, "psi"))
  dipw_est_mat <- do.call(cbind, lapply(ipw_est, `[[`, "dipw"))

  # compute the estimated DCAR and EIF along the regularization sequence
  dcar_est_mat <- est_dcar(psi_ipw_lambda, gn_pred_natural, gn_pred_shifted,
                           Qn_pred_natural, Qn_pred_shifted)
  eif_est_mat <- dipw_est_mat - dcar_est_mat

  # compute empirical mean of DCAR for the hybrid plateau selector
  dcar_min_idx <- which.min(abs(matrixStats::colMeans2(dcar_est_mat)))

  # empirical variances of the estimating equation and EIF
  var_ipw_dipw <- matrixStats::colVars(dipw_est_mat) / n_obs
  var_ipw_eif <- matrixStats::colVars(eif_est_mat) / n_obs

  # NOTE: use the conservative estimating equation variance (instead of the
  #       EIF variance) since we only need *changes in SE* for selectors
  psi_ipw_selector <- unname(psi_ipw_lambda[seq_along(lambda_usm_seq)])
  var_ipw_selector <- var_ipw_dipw[seq_along(lambda_usm_seq)]
  se_ipw <- sqrt(var_ipw_selector)

  # NOTE: selector #1 -- Lepski-type plateau of changes in psi vs. SE
  # trading off changes in psi vs. Z_{1-alpha/2}/log(n) * change in SE, where
  # the inclusion of sqrt(n) allows optimal selection wrt MSE asymptotically
  psi_ipw_diff <- matrixStats::diff2(psi_ipw_selector)
  se_ipw_diff <- matrixStats::diff2(se_ipw)
  lepski_crit <- abs(psi_ipw_diff) <=
    (ci_mult / log10(n_obs)) * abs(se_ipw_diff)
  lepski_idx <- which.max((lepski_crit))

  # Lepski-select IPW estimator and get SE based on EIF
  psi_lepski <- psi_ipw_selector[lepski_idx]
  se_lepski <- sqrt(var_ipw_eif)[lepski_idx]
  lambda_lepski <- lambda_usm_seq[lepski_idx]
  ip_wts_lepski <- ip_wts_mat[, lepski_idx]
  eif_lepski <- eif_est_mat[, lepski_idx]

  # NOTE: selector #2 -- plateau in the point estimate, _loosely_ inspired by
  # ideas for sieve variance estimation (cf. Molly Davies's dissertation ch 3)
  # 1) narrow down to only those regions with L1 norm near CV-selected L1 norm
  #    and within a standard error or so of the CV-selected IPW point estimate
  cv_psi_bounds <- psi_ipw_selector[1] + c(-1, 1) * ci_mult * se_ipw[1]
  cv_l1norm_check <- l1norm_usm_grid <= (l1norm_usm_grid[1] * l1norm_mult)
  psi_bounds_check <- dplyr::between(
    psi_ipw_selector, min(cv_psi_bounds), max(cv_psi_bounds)
  )
  psi_plateau_region <- as.logical(cv_l1norm_check * psi_bounds_check)

  # add region for hybrid selector, using DCAR minimizer as stopping point
  hybrid_plateau_region <- psi_plateau_region
  if (dcar_min_idx > max(which(psi_plateau_region))) {
    hybrid_plateau_region[seq_len(dcar_min_idx)] <- TRUE
  }

  # compute IPW selection for both psi-based and D_CAR-based plateau grids
  plateau_region_grid <- cbind(psi_plateau_region, hybrid_plateau_region)
  plateau_results <- apply(plateau_region_grid, 2, function(plateau_region) {
    # 2) smoothen undersmoothing trajectory within selected window
    loess_plateau_region <- stats::loess(
      psi_ipw_selector[plateau_region] ~ l1norm_usm_grid[plateau_region],
      span = 0.4, degree = 2L
    )
    psi_ipw_smooth <- stats::predict(loess_plateau_region)

    # 3) select first inflection point of the smoothed regional trajectory
    ddx_psi_ipw <- matrixStats::diff2(psi_ipw_smooth)
    ddx2_psi_ipw <- matrixStats::diff2(psi_ipw_smooth, differences = 2L)
    inflects <- intersect(
      # NOTE: 1st derivative should change sign
      which(matrixStats::diff2(sign(ddx_psi_ipw)) != 0),
      # NOTE: 2nd derivative should be non-zero
      which(abs(ddx2_psi_ipw) > 0)
    )
    # NOTE: if no inflection point in window, default to CV-selected choice
    if (length(inflects) > 0) {
      # set chosen L1-norm value to inflection point
      plateau_idx <- min(inflects)

      # plateau-select IPW estimator and get SE based on EIF
      psi_plateau <- psi_ipw_selector[plateau_idx]
      se_plateau <- sqrt(var_ipw_eif)[plateau_idx]
    } else {
      # set CV-selected value of L1-norm
      plateau_idx <- 1L

      # CV-selected point and SE estimates
      psi_plateau <- psi_ipw_selector[plateau_idx]
      se_plateau <- sqrt(var_ipw_eif)[plateau_idx]
      #psi_plateau <- psi_plateau + ci_mult * se_plateau
    }

    # set regularization value, weights, and EIF based on plateau selection
    lambda_plateau <- lambda_usm_seq[plateau_idx]
    ip_wts_plateau <- ip_wts_mat[, plateau_idx]
    eif_plateau <- eif_est_mat[, plateau_idx]

    # organize and return output
    out <- list(lambda = lambda_plateau, idx = plateau_idx, eif = eif_plateau,
                ip_wts = ip_wts_plateau, psi = psi_plateau, se = se_plateau)
    return(out)
  })
  psi_plateau_results <- plateau_results[["psi_plateau_region"]]
  hybrid_plateau_results <- plateau_results[["hybrid_plateau_region"]]

  # bundle each selector's choice together for return object
  est_mat <- tibble::as_tibble(list(
    psi = c(psi_lepski, psi_plateau_results$psi, hybrid_plateau_results$psi),
    se_est = c(se_lepski, psi_plateau_results$se, hybrid_plateau_results$se),
    lambda_idx = c(lepski_idx, psi_plateau_results$idx,
                   hybrid_plateau_results$idx),
    type = c("lepski_plateau", "smooth_plateau", "hybrid_plateau")
  ))
  eif_mat <- cbind(eif_lepski, psi_plateau_results$eif,
                   hybrid_plateau_results$eif)
  colnames(eif_mat) <- c("lepski_plateau", "smooth_plateau", "hybrid_plateau")
  eif_mat <- tibble::as_tibble(eif_mat)

  # output
  out <- list(est = est_mat, eif = eif_mat)
  return(out)
}
