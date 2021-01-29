# author: David Benkeser

library(data.table)
library(ggplot2)
library(dplyr)
library(hal9001)
devtools::load_all()
set.seed(76924)

# simulate data: W ~ Rademacher and A|W ~ N(mu = \pm 1, sd = 0.5)
n_train <- 100
w <- runif(n_train, -4, 4)
a <- rnorm(n_train, w, 0.5)

# learn relationship A|W using HAL-based density estimation procedure
# tune over different choices of n_bins and grid_type
haldensify_fit <- haldensify(
  A = a, W = w,
  grid_type = c("equal_range","equal_mass"),
  n_bins = c(5, 10),
  lambda_seq = exp(seq(-1, -13, length = 250))
)

# predictions to recover conditional density of A|W
new_a <- seq(-8, 8, by = 0.01)
w_val <- c(-3, -1, 1, 3)
add_line <- function(a_val = seq(-5, 5, by = 0.01),
                     w_val = 0, new_plot = TRUE, 
                     col_true = 1, col_est = 1, ...) {
  pred <- predict(haldensify_fit,
  new_A = a_val, new_W = rep(w_val, length(a_val)))
  if (new_plot) {
    plot(0, 0, pch = "", xlim = range(a_val), ylim = c(0, max(pred) * 1.10),
         xlab = "a", ylab = "Density")
  }
  # add true density
  lines(x = a_val, y = dnorm(a_val, w_val, 0.5), lty = 2, col = col_true)

  # add predicted density
  lines(x = a_val, y = pred, lty = 1, col = col_est)
}

add_line(col_true = 1, col_est = 1)
add_line(w_val = -2, new_plot = FALSE, col_true = 2, col_est = 2)
add_line(w_val = 2, new_plot = FALSE, col_true = 4, col_est = 4)
