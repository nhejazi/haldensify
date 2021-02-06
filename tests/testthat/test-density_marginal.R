set.seed(76924)

# simulate data: A ~ N(mu = 0, sd = 2)
n_train <- 500
a <- rnorm(n_train, 0, 2)

# learn marginal density of A using HAL regression
mod_haldensify <- haldensify(
  A = a, W = NULL,
  n_bins = c(3, 5),
  lambda_seq = exp(seq(-1, -13, length = 200))
)

# estimate density via Gaussian kernel density
gauss_dens <- density(a)

# HAL predictions of density using support from the fitted density function
hal_dens <- predict(mod_haldensify,
  new_A = gauss_dens$x, new_W = rep(0, length(a))
)

# compute empirical risk of each estimator based on -log(P) loss
hal_emprisk <- mean(-log(hal_dens[hal_dens > 0]))
gauss_emprisk <- mean(-log(gauss_dens$y)[hal_dens > 0])

# HAL-predicted density matches Gaussian kernel density closely in loss?
test_that("Empirical risk of HAL density less than that of kernel density", {
  expect_lt(hal_emprisk, gauss_emprisk)
})
