# read in command line arguments
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE,
                             defaults = list(n_sim = 200,
                                             dgp = "1a",
                                             delta = 1,
                                             param = "tsm",  #"pie",
                                             estim = "tmle",
                                             Qn_type = "mean"))

# reference for logging
print(args)

# packages, housekeeping
library(here)
library(data.table)
library(tidyverse)
library(haldensify)
library(future)
devtools::load_all(here("ipwshift"))

# load helper package and scripts
source(here("R", "01_dgp.R"))
source(here("R", "02_estimate.R"))
set.seed(72439)
plan(multicore, workers = as.integer(availableCores() - 2))

# sample sizes
n_samp <- cumsum(rep(sqrt(100), 3))^2

# get truth for reference
truth <- get_truth(delta = args$delta, dgp_type = args$dgp,
                   param_type = args$param)

# just pick a sample size
n_obs <- n_samp[2]

# set the number of example experimental datasets to make
n_exp <- 4

for (iter in seq_len(n_exp)) {
  # generate data for this simulation
  dgp <- make_sim_data(n_obs, dgp_type = args$dgp)
  data_obs <- dgp$data_obs
  w_names <- str_subset(colnames(data_obs), "W")

  # set variables for estimation
  W <- data_obs[, w_names]
  A <- data_obs$A
  Y <- data_obs$Y
  delta <- args$delta
  Qn_type <- args$Qn_type
  eff_est <- args$estim
  param <- args$param

  # outcome family for selectors that require Qn
  outcome_family <- ifelse(length(unique(Y)) > 2, "gaussian", "binomial")

  # set simulation parameters
  n_boot <- 500
  n_bins <- c(3, 5)
  cv_folds <- 5
  lambda_seq <- exp(seq(-0.01, -50, length = 2e3))

  # fit haldensify for gPS (generalized propensity score)
  gn_fit_haldensify <- haldensify::haldensify(
    A = A, W = W,
    n_bins = n_bins,
    cv_folds = cv_folds,
    lambda_seq = lambda_seq
  )

  # extract gPS predictions for "natural" A and counterfactual (down)shifted A
  gn_pred_natural <- stats::predict(
    gn_fit_haldensify,
    new_A = A,
    new_W = W,
    lambda_select = "undersmooth"
  )

  gn_pred_shifted <- stats::predict(
    gn_fit_haldensify,
    new_A = (A - delta),
    new_W = W,
    lambda_select = "undersmooth"
  )

  # useful constants and IP weights
  n_obs <- length(Y)
  ip_wts <- gn_pred_shifted / gn_pred_natural

  # construct bootstrap samples
  data_obs <- tibble::tibble(W, A, Y)
  boot_samples <- rsample::bootstraps(data = data_obs, times = n_boot)

  # fit haldensify on each of the bootstrap re-samples
  # NOTE: pass in CV-selected tuning parameters and basis list for HAL fits
  # TODO: consider switching lapply loop to future_lapply, nb, this complicates
  #       the simulation process since future topologies must be properly set.
  haldensify_pred_boot <-
    future.apply::future_lapply(seq_along(boot_samples$splits),
                                function(boot_idx) {
    # get bootstrap sample
    boot_samp <- boot_samples$splits[[boot_idx]]
    data_boot <- tibble::as_tibble(boot_samp)

    # fit HAL model on bootstrap sample
    haldensify_boot <- haldensify::haldensify(
      A = data_boot$A,
      W = data_boot %>% dplyr::select(dplyr::contains("W")),
      grid_type = gn_fit_haldensify$grid_type_cvselect,
      n_bins = gn_fit_haldensify$n_bins_cvselect,
      cv_folds = cv_folds,
      lambda_seq = gn_fit_haldensify$hal_fit$lambda_star,
      hal_basis_list = gn_fit_haldensify$hal_fit$basis_list
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
    print(paste("Completed estimation on bootstrap sample", boot_idx))
    gn_preds_out <- list(gn_natural = gn_pred_natural_boot,
                         gn_shifted = gn_pred_shifted_boot)
    return(gn_preds_out)
  }, future.seed = TRUE)

  # reduce predicted CDE from HAL models fit on bootstrap samples to only those
  # lambdas selected as part of the undersmoothed sequence (CV -> more relaxed)
  trim_boot_preds <- function(hal_boot_preds, hal_orig_preds,
                              type = c("natural", "shifted")) {
    # extract predicted CDE
    hal_boot_preds_typed <- hal_boot_preds[[paste("gn", type, sep = "_")]]

    # get subset of lambdas by matching column names
    lambda_col_idx <- which(colnames(hal_boot_preds_typed) %in%
                            colnames(hal_orig_preds))
    return(hal_boot_preds_typed[, lambda_col_idx])
  }
  gn_pred_natural_boot <- future.apply::future_lapply(
    haldensify_pred_boot, trim_boot_preds, gn_pred_natural, "natural"
  )
  gn_pred_shifted_boot <- future.apply::future_lapply(
    haldensify_pred_boot, trim_boot_preds, gn_pred_shifted, "shifted"
  )

  # pack density estimates on original and bootstrap samples into a list
  gn_pred_natural_all <- c(list(gn_pred_natural), gn_pred_natural_boot)
  names(gn_pred_natural_all) <- c("original", boot_samples$id)

  gn_pred_shifted_all <- c(list(gn_pred_shifted), gn_pred_shifted_boot)
  names(gn_pred_shifted_all) <- c("original", boot_samples$id)

  # NOTE: stopping point for example datasets for methods formulation
  file_name <- paste0("dgp", args$dgp, "_haldensify_bootstrap_testing123_",
                      "exp", iter, ".rds")
  if (str_detect(Sys.info()[["nodename"]], "brc")) {
    save.image(paste0("/global/scratch/nhejazi/haldensify_sims/data/",
                      file_name))
  } else {
    save.image(here("data", file_name))
  }
}
