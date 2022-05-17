#' Confidence Intervals for IPW Estimates of the Causal Effects of Stochatic
#' Shift Interventions
#'
#' @details Compute confidence intervals for estimates produced by
#'  \code{\link{ipw_shift}}.
#'
#' @param object An object of class \code{ipw_haldensify}, produced by invoking
#'  the function \code{\link{ipw_shift}}, for which a confidence interval is to
#'  be computed.
#' @param parm A \code{numeric} vector indicating indices of \code{object$est}
#'  for which to return confidence intervals.
#' @param level A \code{numeric} indicating the nominal level of the confidence
#'  interval to be computed.
#' @param ... Other arguments. Not currently used.
#'
#' @method confint ipw_haldensify
#'
#' @importFrom stats qnorm plogis qlogis
#' @importFrom tibble as_tibble add_column
#'
#' @return A named \code{numeric} vector containing the parameter estimate from
#'  a \code{ipw_haldensify} object, alongside lower/upper Wald-style confidence
#'  intervals at a specified coverage level.
#'
#' @export
#'
#' @examples
#' # simulate data
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
#'   delta = 0.5, cv_folds = 2L,
#'   n_bins = 5L, bin_type = "equal_range",
#'   lambda_seq = exp(seq(-1, -10, length = 100L)),
#'   # arguments passed to hal9001::fit_hal()
#'   max_degree = 3,
#'   smoothness_orders = 0,
#'   num_knots = NULL,
#'   reduce_basis = 1 / sqrt(n_obs)
#' )
#' confint(est_ipw)
confint.ipw_haldensify <- function(object,
                                   parm = seq_len(object$psi),
                                   level = 0.95,
                                   ...) {
  # first, let's get Z_(1 - alpha/2)
  ci_mult <- abs(stats::qnorm(p = (1 - level) / 2))

  if (object$.outcome_levels > 2) { # assume continuous outcome
    # compute the interval around the point estimate
    ci_upr_psi <- object$est$psi + object$est$se_est * ci_mult
    ci_lwr_psi <- object$est$psi - object$est$se_est * ci_mult
  } else if (object$.outcome_levels == 2) { # binary outcome
    # for binary outcome case, compute on the logit scale
    psi_ratio <- stats::qlogis(object$est$psi)
    grad_ratio_delta <- (1 / object$est$psi) + (1 / (1 - object$est$psi))
    se_eif_logit <- sqrt(grad_ratio_delta^2 * object$est$se_est^2)

    # compute the lower/upper confidence limits, then back-transform
    ci_upr_psi <- stats::plogis(psi_ratio + se_eif_logit * ci_mult)
    ci_lwr_psi <- stats::plogis(psi_ratio - se_eif_logit * ci_mult)
  } else {
    stop("The outcome has fewer than 2 levels: this case is not supported.")
  }

  # append computed CIs to custom class object
  est_with_cis <- object$est %>%
    tibble::add_column(
      lwr_ci = ci_lwr_psi,
      .before = "psi"
    ) %>%
    tibble::add_column(
      upr_ci = ci_upr_psi,
      .after = "psi"
    )
  return(est_with_cis)
}
