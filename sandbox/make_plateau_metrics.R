# packages and settings
library(here)
library(tidyverse)
set.seed(72439)

# preset number of example experimental datasets
n_exp <- 4

for (iter in seq_len(n_exp)) {
  # load example analysis dataset
  file_name <- paste0("dgp1a_haldensify_bootstrap_testing123_",
                      "exp", iter, ".rds")
  if (str_detect(Sys.info()[["nodename"]], "brc") |
      str_detect(Sys.info()[["nodename"]], "savio3")) {
    load(paste0("/global/scratch/nhejazi/haldensify_sims/data/", file_name))
  } else {
    load(here("data", file_name))
  }

  # build IPW estimate for regularization sequence for each bootstrap sample
  psi_ipw_lambda <- lapply(seq_along(gn_pred_natural_all), function(boot_idx) {
    ip_wts <- gn_pred_shifted_all[[boot_idx]] / gn_pred_natural_all[[boot_idx]]
    psi_ipw_lambda <- apply(ip_wts, 2, function(ip_wts_lambda) {
      stats::weighted.mean(Y, ip_wts_lambda)
    })
    return(psi_ipw_lambda)
  })
  psi_ipw_lambda_boot <- do.call(rbind, psi_ipw_lambda)
  rownames(psi_ipw_lambda_boot) <- names(gn_pred_natural_all)

  # compute bootstrap MSE for each lambda in the sequence
  var_ipw_boot <- matrixStats::colVars(psi_ipw_lambda_boot[-1, ])
  mse_ipw_boot <- matrixStats::colMeans2(
    (psi_ipw_lambda_boot[-1, ] - psi_ipw_lambda_boot[1, ])^2
  )

  # create plotting helpers
  cv_lambda_idx <- gn_fit_haldensify$cv_tuning_results$lambda_loss_min_idx
  lambda_seq <- gn_fit_haldensify$cv_tuning_results$lambda_seq
  lambda_usm_seq <- lambda_seq[cv_lambda_idx:length(lambda_seq)]
  boot_metrics <- list(
    lambda = lambda_usm_seq,
    mse_boot = mse_ipw_boot,
    var_boot = var_ipw_boot
  ) %>%
  as_tibble()
  plot_data_out <- list(
    psi = psi_ipw_lambda_boot,
    boot = boot_metrics
  )

  # save summary to repo regardless of host
  file_name <- paste0(str_remove(file_name, ".rds"), "_summary", ".rds")
  saveRDS(plot_data_out, here("data", file_name))
}
