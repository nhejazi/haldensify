set.seed(11249)

# simulate data: W ~ Rademacher and A|W ~ N(mu = ±1, sd = 0.5)
n_train <- 100
w <- rbinom(n_train, 1, 0.5)
w[w == 0] <- -1
a <- rnorm(n_train, 2 * w, 0.5)

# specify cut-points in A manually
break_min <- floor(min(a))
break_max <- ceiling(max(a))
breaks <- seq(break_min, break_max, by = 1)

# create pooled hazard data frame with manual cut-points
hazard_df <- format_long_hazards(
  A = a, W = w, breaks = breaks
)

test_that("Re-formatting data for hazard estimation with specified breaks", {
  expect_equal(breaks, hazard_df$breaks)
})

test_that("Density estimation with specified breaks", {
  # try to learn the conditional density with pre-set break points
  haldensify_fit <- haldensify(
    A = a, W = w,
    n_bins = NULL,
    breaks = breaks,
    lambda_seq = exp(seq(-1, -13, length = 100)),
    max_degree = 2
  )

  expect_equal(haldensify_fit$breaks, breaks)
})

test_that("Predict from fitted density model with specified breaks", {
  # try to learn the conditional density with pre-set break points
  haldensify_fit2 <- haldensify(
    A = a, W = w,
    lambda_seq = exp(seq(-1, -13, length = 100)),
    max_degree = 2
  )

  # new data for prediction
  new_a <- seq(round(min(a)), round(max(a)), by = 0.1)
  new_w <- rep(1, length(new_a))
  haldensify_pred <- predict(haldensify_fit2, new_A = new_a, new_W = new_w)
})

# try again with uneven binning
breaks_neg <- seq(break_min, 0, by = 0.5)
breaks_pos <- seq(0, break_max, by = 0.25)
breaks <- unique(c(breaks_neg, breaks_pos))

# create pooled hazard data frame with manual cut-points
hazard_df <- format_long_hazards(
  A = a, W = w, breaks = breaks
)
