library(data.table)
set.seed(76921)

# data simulation
sim_data_set <- function(n_obs = 1000) {
  w <- rbinom(n_obs, 1, 0.5)
  w[w == 0] <- -1
  prob_delta <- runif(n_obs, 0.05, 1)
  delta <- rbinom(n = n_obs, size = 1, prob = prob_delta)
  a <- rnorm(n = n_obs, mean = 2 * w, sd = 0.5)
  data_in <- as.data.table(cbind(a, delta, w, 1 / prob_delta))
  data_in <- data_in[delta == 1, ]
  data_in[, delta := NULL]
  setnames(data_in, c("a", "w", "wts"))
  return(data_in)
}
data_in <- sim_data_set(n_obs = 100)

# learn relationship A|W using HAL-based density estimation procedure
dens_lrn <- with(
  data_in,
  haldensify(
    A = a, W = w,
    wts = wts,
    n_bins = c(3, 5),
    lambda_seq = exp(seq(-1, -8, length = 100)),
    max_degree = 3, smoothness_orders = 0
  )
)

# predictions to recover conditional density of A, given W = 0 or W = 1
new_a <- seq(-4, 4, by = 0.1)
new_w_neg <- rep(-2, length(new_a))
new_w_pos <- rep(2, length(new_a))
new_dat <- as.data.table(
  list(a = new_a, w_neg = new_w_neg, w_pos = new_w_pos)
)
new_dat$pred_w_neg <- predict(dens_lrn,
  new_A = new_dat$a, new_W = new_dat$w_neg,
  trim = FALSE
)
new_dat$pred_w_pos <- predict(dens_lrn,
  new_A = new_dat$a, new_W = new_dat$w_pos,
  trim = FALSE
)

# test that maximum value of prediction happens at appropriate mean of the
# conditional density N(mu = \pm 2, sd = 0.5)
# test_that("Maximum predicted probability of p(A|W = -1) matches N(-2, 0.5)", {
# expect_equal(
# new_dat[which.max(pred_w_neg), a],
# data_in[w == -1, mean(a)]
# )
# })

# test_that("Maximum predicted probability of p(A|W = +1) matches N(+2, 0.5)", {
# expect_equal(
# new_dat[which.max(pred_w_pos), a],
# data_in[w == 1, mean(a)]
# )
# })
