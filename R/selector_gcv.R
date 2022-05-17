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
#' @importFrom stats var weighted.mean
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

  # CV-selected HAL fit is the first by definition
  cv_idx <- 1L

  # compute IPW estimator from CV-selected HAL fits
  psi_ipw <- apply(ip_wts, 2, function(x) stats::weighted.mean(Y, x))
  dipw_cv <- ip_wts[, cv_idx] * (Y - psi_ipw[cv_idx])

  # NOTE: this is technically the variance of the IPW estimator based on its
  #       estimating function, but we ignore it and use instead the variance
  #       estimate based on the EIF, which ought to be _conservative_ since
  #       undersmoothing debiases to improve efficiency
  # var_ipw_cv <- stats::var(dipw_cv) / n_obs

  # compute the D_CAR projection, the EIF, and variance from EIF
  # NOTE: D_CAR = Q(a + delta, w) - [g* / g] Q(a, w) - psi [(g - g*) / g]
  dcar_est <- est_dcar(
    psi_ipw, gn_pred_natural, gn_pred_shifted,
    Qn_pred_natural, Qn_pred_shifted
  )
  dcar_cv <- dcar_est[, cv_idx]
  eif_cv <- dipw_cv - dcar_cv
  var_ipw_cv <- stats::var(eif_cv) / n_obs

  # organize output into tibble
  est <- list(
    psi = psi_ipw[cv_idx],
    se_est = sqrt(var_ipw_cv),
    lambda_idx = cv_idx,
    type = "gcv"
  ) %>% tibble::as_tibble()

  # output simple list
  out <- list(est = est, eif = eif_cv)
  return(out)
}
