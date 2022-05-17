#' Targeted IPW Estimator Selector via Solving the Efficient Influence Function
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
#'  at the natural (i.e., observed) values of the exposure.
#' @param gn_pred_shifted A \code{matrix} of conditional density estimates of
#'  the exposure mechanism g(A+delta|W) along a grid of the regularization
#'  parameter, at the shifted (i.e., counterfactual) values of the exposure.
#' @param Qn_pred_natural A \code{numeric} of the outcome mechanism estimate at
#'  the natural (i.e., observed) values of the exposure. HAL regression is used
#'  for the estimate, with the regularization term chosen by cross-validation.
#' @param Qn_pred_shifted A \code{numeric} of the outcome mechanism estimate at
#'  the shifted (i.e., counterfactual) values of the exposure. HAL regression
#'  is used for the estimate, with the regularization term chosen by
#'  cross-validation.
#'
#' @importFrom matrixStats colMeans2
#' @importFrom stats predict var weighted.mean
#' @importFrom tibble as_tibble
#' @importFrom dplyr "%>%" bind_rows
#' @importFrom rlang set_names
#'
#' @keywords internal
dcar_selector <- function(W, A, Y,
                          delta = 0,
                          gn_pred_natural,
                          gn_pred_shifted,
                          Qn_pred_natural,
                          Qn_pred_shifted) {
  # useful constants and IP weights
  n_obs <- length(Y)
  ip_wts <- gn_pred_shifted / gn_pred_natural

  # compute IPW estimators for IP weights along the regularization sequence
  ipw_est <- apply(ip_wts, 2, function(ipw) {
    psi_ipw <- stats::weighted.mean(Y, ipw)
    dipw_score <- ipw * (Y - psi_ipw)
    return(list(psi = psi_ipw, dipw = dipw_score))
  })
  psi_ipw <- do.call(c, lapply(ipw_est, `[[`, "psi"))
  dipw_est_eqn <- do.call(cbind, lapply(ipw_est, `[[`, "dipw"))

  # extract IPW estimator for CV-selected choice of regularization parameter
  psi_ipw_cv <- unname(psi_ipw[1])
  dipw_cv <- dipw_est_eqn[, 1]

  # compute D_CAR projection and use this to compute the EIF
  # NOTE: D_CAR = Q(a+delta,w) - [g*/g]Q(a,w) - psi [(g-g*)/g]
  dcar_est <- est_dcar(psi_ipw, gn_pred_natural, gn_pred_shifted,
                       Qn_pred_natural, Qn_pred_shifted)
  eif_est <- dipw_est_eqn - dcar_est
  var_ipw_cv <- stats::var(eif_est[, 1]) / n_obs

  # NOTE: using the D_CAR minimizer is too stringent so we instead use a
  #       criterion like |Pn D_CAR| <= sqrt(var(EIF_cv) / n) / log(n)
  tol_crit <- sqrt(var_ipw_cv) / log(n_obs)
  dcar_crit_idx <- which.max(abs(matrixStats::colMeans2(dcar_est)) <= tol_crit)
  dcar_min_idx <- which.min(abs(matrixStats::colMeans2(dcar_est)))
  dcar_selector_idx <- c(dcar_crit_idx, dcar_min_idx)

  # compute undersmoothed IPW estimator
  dcar_est <- lapply(dcar_selector_idx, function(us_idx) {
    psi_ipw_us <- stats::weighted.mean(Y, ip_wts[, us_idx])
    dipw_est_us <- ip_wts[, us_idx] * (Y - psi_ipw_us)
    eif_est_us <- dipw_est_us - dcar_est[, us_idx]
    var_ipw_us <- stats::var(eif_est_us) / n_obs

    # organize output
    est <- list(
      psi = psi_ipw_us,
      se_est = sqrt(var_ipw_us),
      lambda_idx = us_idx,
      type = "dcar"
    ) %>% tibble::as_tibble()
    return(list(est = est, eif = eif_est_us))
  })

  # output simple list
  est <- lapply(dcar_est, `[[`, "est") %>%
    rlang::set_names(c("dcar_tol", "dcar_min")) %>%
    dplyr::bind_rows(.id = "type")
  eif_us <- do.call(cbind, lapply(dcar_est, `[[`, "eif"))
  colnames(eif_us) <- c("dcar_tol", "dcar_min")
  out <- list(est = est, eif = eif_us)
  return(out)
}

#' DCAR Projection of the Efficient Influence Function
#'
#' @param psi_ipw A \code{numeric} vector of IPW point estimates corresponding
#'  to the conditional density estimates of the exposure mechanism provided in
#'  the arguments \code{gn_pred_natural} and \code{gn_pred_shifted}.
#' @param gn_pred_natural A \code{matrix} of conditional density estimates of
#'  the exposure mechanism g(A|W) along a grid of the regularization parameter,
#'  at the natural (i.e., observed) values of the exposure.
#' @param gn_pred_shifted A \code{matrix} of conditional density estimates of
#'  the exposure mechanism g(A+delta|W) along a grid of the regularization
#'  parameter, at the shifted (i.e., counterfactual) values of the exposure.
#' @param Qn_pred_natural A \code{numeric} of the outcome mechanism estimate at
#'  the natural (i.e., observed) values of the exposure. HAL regression is used
#'  for the estimate, with the regularization term chosen by cross-validation.
#' @param Qn_pred_shifted A \code{numeric} of the outcome mechanism estimate at
#'  the shifted (i.e., counterfactual) values of the exposure. HAL regression
#'  is used for the estimate, with the regularization term chosen by
#'  cross-validation.
#'
#' @keywords internal
est_dcar <- function(psi_ipw,
                     gn_pred_natural, gn_pred_shifted,
                     Qn_pred_natural, Qn_pred_shifted) {
  # compute the IP weights
  ip_wts <- gn_pred_shifted / gn_pred_natural

  # compute D_CAR projection for use in finding the minimizing lambda
  # NOTE: D_CAR = Q(d(a,w),w) - [g*/g]Q(a,w) - psi [(g-g*)/g]
  gn_resid <- (gn_pred_natural - gn_pred_shifted) / gn_pred_natural
  gn_resid_rescaled <- lapply(seq_along(psi_ipw), function(lambda_idx) {
    psi_ipw[lambda_idx] * gn_resid[, lambda_idx]
  })
  gn_resid_rescaled <- do.call(cbind, gn_resid_rescaled)
  dcar_est <- Qn_pred_shifted - (ip_wts * Qn_pred_natural) - gn_resid_rescaled
  return(dcar_est)
}
