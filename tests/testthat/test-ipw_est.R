library(data.table)
library(stringr)
set.seed(76924)

# load data simulation helpers
source("helpers_dgp.R")

# simulate test data
dgp <- "2a"
delta <- 1
n_samp <- 100
data_with_dgp <- make_sim_data(n_samp = n_samp, dgp_type = dgp)
data_obs <- data_with_dgp$data_obs

# estimate TSM for all types of undersmoothing selectors
est_ipw <- ipw_shift(
  W = data_obs[, c("W1", "W2", "W3")],
  A = data_obs$A, Y = data_obs$Y,
  delta = delta,
  n_bins = 5L,
  cv_folds = 3L,
  lambda_seq = exp(seq(-1, -8, length = 100L)),
  ## arguments passed to hal9001::fit_hal()
  max_degree = 3,
  smoothness_orders = 0,
  # reduce_basis = 1 / n_samp,
  ## ...more non-hal9001 args
  bin_type = "equal_range",
  undersmooth_type = "all"
)

# generate confidence intervals for IPW estimates
ci_ipw <- confint(est_ipw)
ci_ipw$type <- est_ipw$est$type

# get approximately true TSM value for this DGP
tsm_true <- get_truth(
  n_samp = 1e4, delta = delta, dgp_type = dgp, param_type = "tsm"
)

# extract IPW estimators of each type from output
est_ipw_tbl <- as.data.table(merge(est_ipw$est, ci_ipw, by = "type"))
est_ipw_tbl[, bias := (psi - tsm_true$psi)]
est_ipw_tbl[, mse := bias^2 + se_est^2]
est_ipw_tbl[, nmse := (mse * n_samp) / tsm_true$eff_bound]

# test: bias not too high (within 5% of truth)
est_ipw_tbl[, bias_rel := abs(bias) / tsm_true$psi]
test_that("Bias of IPW based on global CV is within 5% of the truth", {
  expect_lt(est_ipw_tbl[str_detect(type, "gcv"), bias_rel], 0.05)
})
test_that("Bias of plateau-based IPW is within 5% of the truth", {
  expect_true(all(est_ipw_tbl[str_detect(type, "plateau"), bias_rel] < 0.05))
})
test_that("Bias of D_CAR-based IPW is within 5% of the truth", {
  expect_true(all(est_ipw_tbl[str_detect(type, "dcar"), bias_rel] < 0.05))
})

# test: MSE is very small (well-behaved estimators)
test_that("MSE of IPW based on global CV is very low", {
  expect_lt(est_ipw_tbl[str_detect(type, "gcv"), mse], 0.01)
})
test_that("MSE of plateau-based IPW is very low", {
  expect_true(all(est_ipw_tbl[str_detect(type, "plateau"), mse] < 0.01))
})
test_that("MSE of D_CAR-based IPW is very low", {
  expect_true(all(est_ipw_tbl[str_detect(type, "dcar"), mse] < 0.01))
})

# test: scaled MSE is larger than efficiency bound (no super-efficiency)
# test_that("Scaled MSE of IPW based on global CV exceeds efficiency bound", {
# expect_gt(est_ipw_tbl[str_detect(type, "gcv"), nmse], 1)
# })
# test_that("Scaled MSE of plateau-based IPW exceeds efficiency bound", {
# expect_true(all(est_ipw_tbl[str_detect(type, "plateau"), nmse] > 1))
# })
# test_that("Scaled MSE of D_CAR-based IPW exceeds efficiency bound", {
# expect_true(all(est_ipw_tbl[str_detect(type, "dcar"), nmse] > 1))
# })

# test: estimators further in lambda sequence (than GCV) have lower bias
test_that("Bias of undersmoothed IPW no worse than of cross-validated IPW", {
  bias_ipw_gcv <- est_ipw_tbl[type == "gcv", unique(abs(bias))]
  bias_ipw_usm <- est_ipw_tbl[type != "gcv", abs(bias)]
  expect_true(all(bias_ipw_usm <= bias_ipw_gcv))
})

# test: estimators further in lambda sequence (than GCV) have better efficiency
test_that("MSE of undersmoothed IPW is better than for cross-validated IPW", {
  # TODO: why do the plateau selectors at lambda_idx = 1 have better variance??
  mse_ipw_gcv <- est_ipw_tbl[type == "gcv", mse]
  mse_ipw_usm <- est_ipw_tbl[type != "gcv", mse]
  lapply(mse_ipw_usm, function(x) expect_equal(x - mse_ipw_gcv, 0, tol = 1e-2))
})

# test: confidence intervals cover truth for all IPW estimators
test_that("Wald-style CIs cover truth for all IPW estimators", {
  covers <- between(tsm_true$psi, est_ipw_tbl$lwr_ci, est_ipw_tbl$upr_ci)
  expect_true(all(covers))
})
