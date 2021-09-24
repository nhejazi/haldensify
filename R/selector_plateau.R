#' IPW Estimator Selector Using Lepski's Plateau Method for the MSE
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
#' @param gcv_mult TODO
#' @param bootstrap A \code{logical} indicating whether the estimator variance
#'  should be approximated using the nonparametric bootstrap. The default is
#'  \code{FALSE}, in which case the empirical variances of the IPW estimating
#'  function and the EIF are used for for estimator selection and for variance
#'  estimation, respectively. When set to \code{TRUE}, the bootstrap variance
#'  is used for both of these purposes instead. Note that the bootstrap is very
#'  computationally intensive and scales relatively poorly.
#' @param n_boot A \code{numeric} giving the number of bootstrap re-samples to
#'  be used in computing the plateau estimator selection criterion. The default
#'  uses 1000 bootstrap samples, though it may be appropriate to use fewer such
#'  samples for experimentation purposes. This is ignored when \code{bootstrap}
#'  is set to \code{FALSE} (its default).
#' @param ... Additional arguments for model fitting to be passed directly to
#'  \code{\link[haldensify]{haldensify}}.
#'
#' @importFrom future.apply future_lapply
#' @importFrom matrixStats colMeans2 colVars diff2
#' @importFrom stats predict qnorm weighted.mean
#' @importFrom dplyr "%>%" contains select
#' @importFrom tibble tibble as_tibble
#' @importFrom rsample bootstraps
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
                             gcv_mult = 50L,
                             bootstrap = FALSE,
                             n_boot = 1000L,
                             ...) {
  # useful constants and IP weights
  n_obs <- length(Y)
  ip_wts_mat <- gn_pred_shifted / gn_pred_natural
  ci_level <- 0.95 # NOTE: *hard-coded* :(
  ci_mult <- abs(stats::qnorm(p = (1 - ci_level) / 2))
  plateau_cutoff <- c(0.20, 0.15, 0.10, 0.05, 0.01) # NOTE: *hard-coded* :(

  # shorten sequence of lambdas to just the undersmoothed values
  cv_lambda_idx <- gn_fit_haldensify$cv_tuning_results$lambda_loss_min_idx
  cv_lambda <- gn_fit_haldensify$cv_tuning_results$lambda_loss_min
  lambda_seq <- gn_fit_haldensify$cv_tuning_results$lambda_seq
  lambda_usm_seq <- lambda_seq[cv_lambda_idx:length(lambda_seq)]
  lambda_usm_seq <- lambda_usm_seq[lambda_usm_seq > (cv_lambda / gcv_mult)]

  # compute the estimated EIF across the regularization sequence
  score_eif_est <- apply(ip_wts_mat, 2, function(ipw_est) {
    dcar_ipw_est <- (ipw_est * Qn_pred_natural) - Qn_pred_shifted
    psi_ipw_est <- stats::weighted.mean(Y, ipw_est)
    dipw_score_est <- ipw_est * (Y - psi_ipw_est)
    eif_est <- dipw_score_est - dcar_ipw_est
    return(list(psi = psi_ipw_est, dipw = dipw_score_est, eif = eif_est))
  })
  psi_ipw_lambda <- do.call(cbind, lapply(score_eif_est, `[[`, "psi"))
  dipw_est_mat <- do.call(cbind, lapply(score_eif_est, `[[`, "dipw"))
  eif_est_mat <- do.call(cbind, lapply(score_eif_est, `[[`, "eif"))

  if (bootstrap) {
    # construct bootstrap samples
    data_obs <- tibble::tibble(W, A, Y)
    boot_samples <- rsample::bootstraps(data = data_obs, times = n_boot)

    # fit haldensify on each of the bootstrap re-samples
    # NOTE: pass in CV-selected tuning parameters and basis list to HAL
    haldensify_pred_boot <-
      future.apply::future_lapply(seq_along(boot_samples$splits),
        function(boot_idx) {
          # get bootstrap sample
          boot_samp <- boot_samples$splits[[boot_idx]]
          data_boot <- tibble::as_tibble(boot_samp)

          # fit HAL model on bootstrap sample
          haldensify_boot <- haldensify(
            A = data_boot$A,
            W = data_boot %>% dplyr::select(dplyr::contains("W")),
            cv_folds = cv_folds,
            n_bins = gn_fit_haldensify$n_bins_cvselect,
            grid_type = gn_fit_haldensify$grid_type_cvselect,
            lambda_seq = gn_fit_haldensify$hal_fit$lambda_star,
            hal_basis_list = gn_fit_haldensify$hal_fit$basis_list,
            ## arguments passed to hal9001::fit_hal()
            ...
          )

          # extract conditional density estimates on bootstrap sample
          gn_pred_natural_boot <- stats::predict(
            haldensify_boot,
            new_A = data_boot$A,
            new_W = data_boot %>% dplyr::select(dplyr::contains("W")),
            lambda_select = "all"
          )
          gn_pred_shifted_boot <- stats::predict(
            haldensify_boot,
            new_A = (data_boot$A - delta),
            new_W = data_boot %>% dplyr::select(dplyr::contains("W")),
            lambda_select = "all"
          )

          # return predicted CDE from HAL model on given bootstrap resample
          gn_preds_out <- list(
            gn_natural = gn_pred_natural_boot,
            gn_shifted = gn_pred_shifted_boot
          )
          return(gn_preds_out)
        },
        future.seed = TRUE
      )

    # reduce predicted CDE from HAL fits on bootstrap samples to only those
    # lambdas selected as part of the undersmoothed sequence
    trim_boot_preds <- function(hal_boot_preds, hal_orig_preds,
                                type = c("natural", "shifted")) {
      # extract predicted CDE
      hal_boot_preds_typed <- hal_boot_preds[[paste("gn", type, sep = "_")]]

      # get subset of lambdas by matching column names
      lambda_col_idx <- which(colnames(hal_boot_preds_typed) %in%
        colnames(hal_orig_preds))
      return(hal_boot_preds_typed[, lambda_col_idx])
    }
    gn_pred_natural_boot <- lapply(
      haldensify_pred_boot, trim_boot_preds, gn_pred_natural, "natural"
    )
    gn_pred_shifted_boot <- lapply(
      haldensify_pred_boot, trim_boot_preds, gn_pred_shifted, "shifted"
    )

    # pack density estimates on original and bootstrap samples into a list
    gn_pred_natural_all <- c(list(gn_pred_natural), gn_pred_natural_boot)
    names(gn_pred_natural_all) <- c("original", boot_samples$id)
    gn_pred_shifted_all <- c(list(gn_pred_shifted), gn_pred_shifted_boot)
    names(gn_pred_shifted_all) <- c("original", boot_samples$id)

    # IPW estimate for regularization sequence for each bootstrap sample
    psi_ipw_lambda_boot <-
      lapply(seq_along(gn_pred_natural_all), function(boot_idx) {
        ip_wts <- gn_pred_shifted_all[[boot_idx]] /
          gn_pred_natural_all[[boot_idx]]
        psi_ipw_lambda <- apply(ip_wts, 2, function(ip_wts_lambda) {
          weighted.mean(Y, ip_wts_lambda)
        })
        return(psi_ipw_lambda)
      })
    psi_ipw_lambda_boot <- do.call(rbind, psi_ipw_lambda_boot)
    rownames(psi_ipw_lambda_boot) <- names(gn_pred_natural_all)

    # compute bootstrap variance and MSE for each lambda in the sequence
    var_ipw_boot <- matrixStats::colVars(psi_ipw_lambda_boot[-1, ])
    mse_ipw_boot <- matrixStats::colMeans2(
      (psi_ipw_lambda_boot[-1, ] - psi_ipw_lambda_boot[1, ])^2
    )

    # set bootstrap-based estimates to use in selectors
    # NOTE: the bootstrap variance is closer to the true variance than its
    #       faster analogs (e.g., the estimating equation variance)
    psi_ipw_selector <- psi_ipw_lambda_boot[, seq_along(lambda_usm_seq)]
    var_ipw_selector <- var_ipw_boot[seq_along(lambda_usm_seq)]
  } else {
    # empirical variances of the estimation equation and EIF
    var_ipw_dipw <- matrixStats::colVars(dipw_est_mat) / n_obs
    var_ipw_eif <- matrixStats::colVars(eif_est_mat) / n_obs

    # set non-bootstrap-based estimates to use in selectors
    # NOTE: use the conservative estimating equation variance (instead of the
    #       EIF variance) since we only need *changes in SE* for selectors
    psi_ipw_selector <- matrix(psi_ipw_lambda[, seq_along(lambda_usm_seq)],
      nrow = 1
    )
    var_ipw_selector <- var_ipw_dipw[seq_along(lambda_usm_seq)]
  }
  se_ipw <- sqrt(var_ipw_selector)

  # NOTE: selector #1 -- Lepski-type plateau of changes in psi vs. SE
  psi_ipw_diff <- matrixStats::diff2(psi_ipw_selector[1, ])
  lepski_crit <- abs(psi_ipw_diff) <= ci_mult * matrixStats::diff2(se_ipw)
  lepski_idx <- which.max(lepski_crit)
  ## select IPW estimator and get SE + D_IPW
  psi_lepski <- psi_ipw_selector[1, lepski_idx]
  se_lepski <- ifelse(bootstrap, sqrt(var_ipw_boot)[lepski_idx],
    sqrt(var_ipw_eif)[lepski_idx]
  )
  lambda_lepski <- lambda_usm_seq[lepski_idx]
  ip_wts_lepski <- ip_wts_mat[, lepski_idx]
  eif_lepski <- eif_est_mat[, lepski_idx]

  # NOTE: selector #2 -- plateau in the point estimate
  #       defn of plateau *hard-coded* above :(
  psi_ipw_diff_rel <- abs(psi_ipw_diff) / cummax(abs(psi_ipw_diff))
  plateau_idx <- do.call(c, lapply(plateau_cutoff, function(cutoff) {
    which.max(psi_ipw_diff_rel <= cutoff)
  }))
  psi_plateau <- psi_ipw_selector[1, plateau_idx]
  if (bootstrap) {
    se_plateau <- sqrt(var_ipw_boot)[plateau_idx]
  } else {
    se_plateau <- sqrt(var_ipw_eif)[plateau_idx]
  }
  lambda_plateau <- lambda_usm_seq[plateau_idx]
  ip_wts_plateau <- ip_wts_mat[, plateau_idx]
  eif_plateau <- eif_est_mat[, plateau_idx]

  # bundle each selector's choice together for return object
  est_mat <- list(
    psi = c(psi_lepski, psi_plateau),
    se_est = c(se_lepski, se_plateau),
    lambda_idx = c(lepski_idx, plateau_idx),
    type = c("lepski_plateau", paste("psi_plateau", plateau_cutoff, sep = "_"))
  ) %>% tibble::as_tibble()

  eif_mat <- cbind(eif_lepski, eif_plateau)
  colnames(eif_mat) <- c(
    "lepski_plateau",
    paste("psi_plateau", plateau_cutoff, sep = "_")
  )
  eif_mat <- eif_mat %>% tibble::as_tibble()

  # output
  out <- list(est = est_mat, eif = eif_mat)
  return(out)
}
