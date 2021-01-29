library(tidymodels)
library(haldensify)
library(abind)

# set up simple DGP
g0 <- function(W1, W2, W3) {
  mu <- 2 * W1 - W2 - W1 * W2
  sigma2 <- 4
  return(list(mu = mu, sigma2 = sigma2))
}
# define outcome mechanism
Q0 <- function(A, W1, W2, W3) {
  plogis(3 * A + W1 + W2 - 2 * W3 - W1 * W3)
}

# simulate data
n_samp <- 500
W1 <- rbinom(n_samp, 1, 0.6)
W2 <- rbinom(n_samp, 1, 0.2)
W3 <- rpois(n_samp, 3)
g_obs <- g0(W1, W2, W3)
A <- rnorm(n_samp, g_obs$mu, sqrt(g_obs$sigma2))
Y <- rbinom(n_samp, 1, Q0(A, W1, W2, W3))
data_obs <- as_tibble(list(W1 = W1, W2 = W2, W3 = W3, A = A, Y = Y))

# fit haldensify for conditional density estimation
v_folds <- 5
haldensify_fit <- haldensify(
  A = data_obs$A,
  W = data_obs %>% select(contains("W")),
  grid_type = "equal_range",
  n_bins = c(3, 5),
  cv_folds = v_folds,
  lambda_seq = exp(seq(-1, -10, length = 50))
)
haldensify_pred <- predict(
  haldensify_fit,
  new_A = data_obs$A,
  new_W = data_obs %>% select(contains("W")),
  lambda_select = "undersmooth"
)

# construct bootstrap samples and fit haldensify on each of these resamples
# NOTE: pass in CV-selected tuning parameters and basis list for HAL fits
n_boot <- 5
boot_samples <- bootstraps(data_obs, times = n_boot)
haldensify_pred_boot <- lapply(boot_samples$splits, function(data_split) {
  # get bootstrap sample
  data_boot <- as_tibble(data_split)

  # fit HAL model on bootstrap sample
  haldensify_boot <- haldensify(
    A = data_boot$A,
    W = data_boot %>% select(contains("W")),
    grid_type = haldensify_fit$grid_type_cvselect,
    n_bins = haldensify_fit$n_bins_cvselect,
    cv_folds = v_folds,
    lambda_seq = haldensify_fit$hal_fit$lambda_star,
    hal_basis_list = haldensify_fit$hal_fit$basis_list
   )

  # extract conditional density estimates on bootstrap sample
  haldensify_pred_boot <- predict(
    haldensify_boot,
    new_A = data_boot$A,
    new_W = data_boot %>% select(contains("W")),
    lambda_select = "all"
  )

  # return predicted CDE from HAL model on given bootstrap resample
  return(haldensify_pred_boot)
})

# reduce predicted CDE from HAL models fit on bootstrap samples to only those
# lambdas selected as part of the undersmoothed sequence
haldensify_pred_boot_usm <- lapply(haldensify_pred_boot, function(hal_pred) {
  # get subset of lambdas by matching column names
  lambda_col_idx <- which(colnames(hal_pred) %in% colnames(haldensify_pred))
  return(hal_pred[, lambda_col_idx])
})

# pack density estimates on original and bootstrap samples into array
test <- abind(haldensify_pred,
              unlist(haldensify_pred_boot_usm, recursive = FALSE),
              along = 0)

