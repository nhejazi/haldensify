library(data.table)
library(ggplot2)
library(dplyr)
library(hal9001)
library(future)
plan(transparent)
set.seed(76921)

# data simulation
sim_data_set <- function(n_obs = 1000, w_prob = 0.5, shift_delta = 0.5) {
  w <- rbinom(n = n_obs, size = 1, prob = w_prob)
  w[w == 0] <- -1
  ipc_delta <- rbinom(n = n_obs, size = 1, prob = plogis(w))
  a <- rnorm(n = n_obs, mean = 2 * w, sd = 0.5)
  y <- a + w + rnorm(n_obs, mean = 0, sd = 1)
  data_in <- as.data.frame(cbind(y, a, ipc_delta, w, 1 / plogis(w))) %>%
    dplyr::filter(ipc_delta == 1) %>%
    dplyr::select(-ipc_delta) %>%
    as.data.table()
  setnames(data_in, c("Y", "A", "W", "Weights"))
  return(data_in)
}
data_in <- sim_data_set(n_obs = 100)

# learn relationship A|W using HAL-based density estimation procedure
dens_lrn <- with(
  data_in,
  haldensify(
    A = A, W = W,
    wts = Weights,
    n_bins = c(5, 10, 15),
    lambda_seq = exp(seq(-1, -13, length = 200)),
    use_future = FALSE
  )
)

# predictions to recover conditional density of A, given W = 0 or W = 1
new_a <- seq(-4, 4, by = 0.05)
new_w_neg <- rep(-1, length(new_a))
new_w_pos <- rep(1, length(new_a))
new_dat <- as.data.table(list(a = new_a, w_neg = new_w_neg, w_pos = new_w_pos))
new_dat$pred_w_neg <- predict(dens_lrn,
  new_A = new_dat$a, new_W = new_dat$w_neg
)
new_dat$pred_w_pos <- predict(dens_lrn,
  new_A = new_dat$a, new_W = new_dat$w_pos
)

# test that maximum value of prediction happens at appropriate mean of the
# conditional density N(mu = \pm 2, sd = 0.5)
test_that("Maximum predicted probability of p(A|W = -1) matches N(-2, 0.5)", {
  obs_a_max_prob_w_neg <- new_dat[which.max(new_dat$pred_w_neg), ]$a
  expect_equal(
    round(obs_a_max_prob_w_neg),
    round(mean(data_in$A[data_in$W == -1]))
  )
})

test_that("Maximum predicted probability of p(A|W = +1) matches N(+2, 0.5)", {
  obs_a_max_prob_w_pos <- new_dat[which.max(new_dat$pred_w_pos), ]$a
  expect_equal(
    round(obs_a_max_prob_w_pos),
    round(mean(data_in$A[data_in$W == 1]))
  )
})
