library(data.table)
set.seed(76924)

# simulate data: W ~ Rademacher and A|W ~ N(mu = \pm 1, sd = 0.5)
n_train <- 100
w <- rbinom(n_train, 1, 0.5)
w[w == 0] <- -1
a <- rnorm(n_train, 2 * w, 0.5)

# learn relationship A|W using HAL-based density estimation procedure
haldensify_fit <- haldensify(
  A = a, W = w,
  n_bins = c(3, 5),
  lambda_seq = exp(seq(-1, -13, length = 100)),
  max_degree = 2
)

# predictions to recover conditional density of A|W
new_a <- seq(-1, 1, by = 0.01)
new_w_neg <- rep(-1, length(new_a))
new_w_pos <- rep(1, length(new_a))
new_dat <- as.data.table(list(a = new_a, w_neg = new_w_neg, w_pos = new_w_pos))
new_dat$pred_w_neg <- predict(haldensify_fit,
  new_A = new_dat$a, new_W = new_dat$w_neg
)
new_dat$pred_w_pos <- predict(haldensify_fit,
  new_A = new_dat$a, new_W = new_dat$w_pos
)

# NOTE: these tests are poorly thought out, so temporarily removing
# test that maximum value of prediction happens at appropriate mean of the
# conditional density N(mu = \pm 1, sd = 0.5)
#test_that("Maximum predicted probability of p(A|W = -1) matches N(-1, 0.5)", {
  #obs_a_max_prob_w_neg <- new_dat[which.max(new_dat$pred_w_neg), ]$a
  #expect_equal(round(obs_a_max_prob_w_neg), unique(new_w_neg))
#})

#test_that("Maximum predicted probability of p(A|W = +1) matches N(+1, 0.5)", {
  #obs_a_max_prob_w_pos <- new_dat[which.max(new_dat$pred_w_pos), ]$a
  #expect_equal(round(obs_a_max_prob_w_pos), unique(new_w_pos))
#})

# supply fit_control additional arguments
n_lambda <- 100L
haldensify_fit_cntrl <- haldensify(
  A = a, W = w,
  n_bins = c(3, 5),
  lambda_seq = exp(seq(-1, -13, length = n_lambda)),
  max_degree = 2,
  fit_control = list(cv_select = TRUE, n_folds = 3L, use_min = TRUE)
)
cv_lambda_idx <- haldensify_fit_cntrl$cv_tuning_results$lambda_loss_min_idx

# prediction with lambda_selected by cross-validation
pred_w_cv <- predict(haldensify_fit_cntrl,
  new_A = new_dat$a, new_W = new_dat$w_neg, lambda_select = "cv"
)

# prediction with lambda_select undersmooth
pred_w_undersmooth <- predict(haldensify_fit_cntrl,
  new_A = new_dat$a, new_W = new_dat$w_neg, lambda_select = "undersmooth"
)
test_that("Prediction for undersmoothed lambda is of correct dimensions", {
  # in case the CV-chosen lambda is the last in the sequence
  if (cv_lambda_idx == n_lambda) {
    # number of rows should match input data nrows
    expect_equal(length(pred_w_undersmooth), nrow(new_dat))

    # first lambda in sequence should be the cross-validation selector's choice
    expect_equal(pred_w_undersmooth, pred_w_cv)
  } else {
    # number of rows should match input data nrows
    expect_equal(nrow(pred_w_undersmooth), nrow(new_dat))

    # number of columns should be less than the full sequence of lambda
    expect_lt(ncol(pred_w_undersmooth), n_lambda)

    # first lambda in sequence should be the cross-validation selector's choice
    expect_equal(pred_w_undersmooth[, 1], pred_w_cv)
  }
})

# prediction with lambda_select all
pred_w_all <- predict(haldensify_fit_cntrl,
  new_A = new_dat$a, new_W = new_dat$w_neg, lambda_select = "all"
)
test_that("Prediction for all lambda is of correct dimensions", {
  # number of rows should match input data nrows
  expect_equal(nrow(pred_w_all), nrow(new_dat))

  # number of columns should match the full sequence of lambda
  expect_equal(ncol(pred_w_all), n_lambda)
})

# print a fit
suppressWarnings(
  print(haldensify_fit_cntrl)
)
