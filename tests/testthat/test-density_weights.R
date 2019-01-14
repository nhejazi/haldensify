library(data.table)
library(ggplot2)
library(dplyr)
library(hal9001)

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
dens_lrn <- haldensify(
  A = data_in$A, W = data_in$W,
  # wts = data_in$Weights,
  grid_type = "equal_range",
  n_bins = 10,
  lambda_seq = exp(seq(-1, -13, length = 1000))
)

# predictions to recover conditional density of A, given W = 0 or W = 1
n_samp <- 5000
A_supp <- seq(-5, 5, length = n_samp)
W0 <- rep(0, n_samp)
predictions_W0 <- predict(dens_lrn, new_A = A_supp, new_W = W0)

hist_W0 <- as.data.table(list(A = A_supp, likelihood = predictions_W0)) %>%
  ggplot(aes(x = A, y = likelihood)) +
  geom_smooth(method = "loess", se = FALSE) +
  ggtitle("Conditional density p(A | W = 0)") +
  theme_bw()
hist_W0

W1 <- rep(1, n_samp)
predictions_W1 <- predict(dens_lrn, new_A = A_supp, new_W = W1)

hist_W1 <- as.data.table(list(A = A_supp, likelihood = predictions_W1)) %>%
  ggplot(aes(x = A, y = likelihood)) +
  geom_smooth(method = "loess", se = FALSE) +
  ggtitle("Conditional density p(A | W = 1)") +
  theme_bw()
hist_W1
