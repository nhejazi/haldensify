set.seed(76924)

# simulate data: A ~ N(mu = 0, sd = 2)
n_train <- 50
a <- rnorm(n_train, 0, 2)

# learn marginal density of A using HAL regression
haldensify_fit <- haldensify(
  A = a, W = NULL,
  n_bins = 5,
  lambda_seq = exp(seq(-1, -13, length = 100)),
  max_degree = 3
)

# estimate density via Gaussian kernel density
gauss_dens <- density(a)

# HAL predictions of density using support from the fitted density function
hal_dens <- predict(haldensify_fit,
  new_A = gauss_dens$x, new_W = rep(1, length(a)),
  trim = FALSE
)

# compute empirical risk of each estimator based on -log(P) loss
gauss_emprisk <- mean(-log(gauss_dens$y))
hal_emprisk <- mean(-log(hal_dens))

# HAL-predicted density matches Gaussian kernel density closely in loss?
test_that("Empirical risk of HAL density less than that of kernel density", {
  expect_lt(hal_emprisk, gauss_emprisk)
})
