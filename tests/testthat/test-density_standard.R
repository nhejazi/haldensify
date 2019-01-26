library(data.table)
library(ggplot2)
library(dplyr)
library(hal9001)
set.seed(76924)

# simulate data: W ~ Rademacher and A|W ~ N(mu = \pm 1, sd = 0.5)
n_train <- 1000
w <- rbinom(n_train, 1, 0.5)
w[w == 0] <- -1
a <- rnorm(n_train, w, 0.5)

# learn relationship A|W using HAL-based density estimation procedure
mod_haldensify <- haldensify(
  A = a, W = w,
  grid_type = "equal_range",
  n_bins = 10,
  lambda_seq = exp(seq(-1, -13, length = 1000))
)

# predictions to recover conditional density of A|W
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
