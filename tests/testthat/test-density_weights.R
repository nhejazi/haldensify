library(data.table)
library(ggplot2)
library(dplyr)
library(hal9001)
library(future)
plan(sequential)
set.seed(76924)

# data simulation
sim_data_set <- function(n_obs = 1000, w_prob = 0.5, shift_delta = 0.5) {
  w <- rbinom(n = n_obs, size = 1, prob = w_prob)
  ipc_delta <- rbinom(n = n_obs, size = 1, prob = plogis(w))
  a <- rnorm(n = n_obs, mean = 2 * w, sd = 1)
  y <- a + w + rnorm(n_obs, mean = 0, sd = 1)
  data_in <- as.data.frame(cbind(y, a, ipc_delta, w, 1 / plogis(w))) %>%
    dplyr::filter(ipc_delta == 1) %>%
    dplyr::select(-ipc_delta) %>%
    as.data.table()
  setnames(data_in, c("Y", "A", "W", "Weights"))
  return(data_in)
}
data_in <- sim_data_set()

# learn relationship A|W using HAL-based density estimation procedure
dens_lrn <- with(
  data_in,
  haldensify(
    A = A, W = W,
    wts = Weights,
    lambda_seq = exp(seq(-1, -13, length = 100))
  )
)

# predictions to recover conditional density of A, given W = 0 or W = 1
new_a <- seq(-2, 2, by = 0.01)
new_w_neg <- rep(-1, length(new_a))
new_w_pos <- rep(1, length(new_a))
new_dat <- as.data.table(list(a = new_a, w_neg = new_w_neg, w_pos = new_w_pos))
new_dat$pred_w_neg <- predict(mod_haldensify,
  new_A = new_dat$a, new_W = new_dat$w_neg
)
new_dat$pred_w_pos <- predict(mod_haldensify,
  new_A = new_dat$a, new_W = new_dat$w_pos
)

# test that maximum value of prediction happens at appropriate mean of the
# conditional density N(mu = \pm 1, sd = 0.5)
test_that("Maximum predicted probability of p(A|W = -1) matches N(-1, 0.5)", {
  obs_a_max_prob_w_neg <- new_dat[which.max(new_dat$pred_w_neg), ]$a
  expect_equal(round(obs_a_max_prob_w_neg), unique(new_w_neg))
})

test_that("Maximum predicted probability of p(A|W = +1) matches N(+1, 0.5)", {
  obs_a_max_prob_w_pos <- new_dat[which.max(new_dat$pred_w_pos), ]$a
  expect_equal(round(obs_a_max_prob_w_pos), unique(new_w_pos))
})
