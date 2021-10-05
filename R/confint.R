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
#' @importFrom tibble as_tibble
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
#' confint(est_ipw_shift)
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

  # set up output CI object
  ci_out <- cbind(ci_lwr_psi, object$est$psi, ci_upr_psi)
  colnames(ci_out) <- c("lwr_ci", "est", "upr_ci")
  ci_out <- tibble::as_tibble(ci_out)
  return(ci_out)
}
