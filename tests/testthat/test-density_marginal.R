set.seed(76924)

# simulate data: A ~ N(mu = 0, sd = 2)
n_train <- 500
a <- rnorm(n_train, 0, 2)

# learn marginal density of A using HAL regression
mod_haldensify <- haldensify(
  A = a, W = NULL,
  n_bins = c(5, 10),
  lambda_seq = exp(seq(-1, -13, length = 200))
)

# estimate density via Gaussian kernel density
gauss_dens <- density(new_a)
gauss_emprisk <- mean(-log(gauss_dens$y))

# HAL predictions of density using support from the fitted density function
hal_dens <- predict(mod_haldensify,
  new_A = gauss_dens$x, new_W = rep(0, length(new_a))
)
hal_emprisk <- mean(-log(hal_dens))

# HAL-predicted density matches Gaussian kernel density closely in loss?
test_that("Empirical risk of HAL density less than that of kernel density", {
  expect_lt(hal_emprisk, gauss_emprisk)
})
