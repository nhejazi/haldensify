library(data.table)
library(ggplot2)
library(dplyr)
library(hal9001)

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
pred_haldensify <- predict(mod_haldensify, new_A = a, new_W = w)

