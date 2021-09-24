#' IPW Estimator Selector via Global Cross-Validation
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
#' @param Qn_pred_natural A \code{numeric} of the outcome mechanism estimate at
#'  the natural (i.e., observed) values of the exposure. HAL regression is used
#'  for the estimate, with the regularization term chosen by cross-validation.
#' @param Qn_pred_shifted A \code{numeric} of the outcome mechanism estimate at
#'  the shifted (i.e., counterfactual) values of the exposure. HAL regression
#'  is used for the estimate, with the regularization term chosen by
#'  cross-validation.
#'
#' @importFrom stats predict var weighted.mean
#' @importFrom tibble as_tibble
#' @importFrom dplyr "%>%"
#'
#' @keywords internal
gcv_selector <- function(W, A, Y,
                         delta = 0,
                         gn_pred_natural,
                         gn_pred_shifted,
                         Qn_pred_natural,
                         Qn_pred_shifted) {
  # useful constants and IP weights
  n_obs <- length(Y)
  ip_wts <- gn_pred_shifted / gn_pred_natural

  # compute IPW estimator from CV-selected HAL fits
  psi_ipw_cv <- stats::weighted.mean(Y, ip_wts[, 1])
  dipw_cv <- ip_wts[, 1] * (Y - psi_ipw_cv)
  # NOTE: this is technically the variance of the IPW estimator based on its
  #       estimating function, but we ignore and use instead the variance
  #       estimate based on the EIF, which ought to be _conservative_ since the
  #       undersmoothing procedure should improve efficiency.
  # var_ipw_cv <- var(dipw_cv) / n_obs

  # compute D_CAR projection and find minimizing lambda
  # D_CAR = Q(a+delta,w) - [g*/g]Q(a,w) - psi [(g-g*)/g]
  gn_resid <- psi_ipw_cv * ((gn_pred_natural[, 1] - gn_pred_shifted[, 1]) /
    gn_pred_natural[, 1])
  dcar_cv <- Qn_pred_shifted - (ip_wts[, 1] * Qn_pred_natural) - gn_resid
  eif_cv <- dipw_cv - dcar_cv
  # NOTE: the EIF variance should be _conservative_ for the CV-IPW
  var_ipw_cv <- stats::var(eif_cv) / n_obs

  # organize output into tibble
  est <- list(
    psi = psi_ipw_cv,
    se_est = sqrt(var_ipw_cv),
    lambda_idx = 1,
    type = "gcv"
  ) %>% tibble::as_tibble()

  # output simple list
  out <- list(est = est, eif = eif_cv)
  return(out)
}
