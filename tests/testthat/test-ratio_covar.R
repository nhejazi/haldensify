library(data.table)
set.seed(76924)

# HAL for A ~ W1/W2, considering W2 as continuous and discrete
n_train <- 100
w1_norm <- rnorm(n_train, 0, 0.5)
w2_cats <- sample(c(1, 2, 3), n_train, replace = TRUE, prob = c(0.1, 0.6, 0.3))
w2_cont <- runif(n_train, 0.25, 0.75)
a_w2_cont <- rnorm(n_train, w1_norm / w2_cont, 0.5)
a_w2_cats <- rnorm(n_train, w1_norm / w2_cats, 0.5)

# learn density A|W based on categorical W2
fit_cats <- haldensify(
  A = a_w2_cats, W = cbind(w1_norm, w2_cats),
  lambda_seq = exp(seq(-1, -20, length = 500)),
  # arguments passed to hal9001::fit_hal()
  max_degree = NULL, num_knots = NULL, smoothness_orders = 0
)
pred_cats_haz <- predict(
  fit_cats, new_A = a_w2_cats, new_W = cbind(w1_norm, w2_cats)
)
emprisk_haz_cats <- mean(-log(pred_cats_haz))
if (FALSE) {
  library(sl3)
  cats_data <- as.data.table(list(W1 = w1_norm, W2 = w2_cats, A = a_w2_cats))
  cats_task <- sl3_Task$new(
    data = cats_data,
    covariates = c("W1", "W2"),
    outcome = "A"
  )

  # compare against GLM-based semiparametric strategy with constant variance
  hose_glm_lrnr <- Lrnr_density_semiparametric$new(
    mean_learner = Lrnr_glm_fast$new(family = gaussian())
  )
  fit_cats_hose <- hose_glm_lrnr$train(cats_task)
  pred_cats_hose <- fit_cats_hose$predict()
  emprisk_hose_cats <- mean(-log(pred_cats_hose))
  expect_lt(emprisk_haz_cats, emprisk_hose_cats)

  # compare against GLM-based semiparametric strategy with fitted variance
  hese_glm_lrnr <- Lrnr_density_semiparametric$new(
    mean_learner = Lrnr_glm_fast$new(family = gaussian()),
    var_learner = Lrnr_glm_fast$new(family = gaussian())
  )
  fit_cats_hese <- hese_glm_lrnr$train(cats_task)
  pred_cats_hese <- fit_cats_hese$predict()
  emprisk_hese_cats <- mean(-log(pred_cats_hese))
  expect_lt(emprisk_haz_cats, emprisk_hese_cats)
}

# learn density A|W based on continuous-valued W2
fit_cont <- haldensify(
  A = a_w2_cont, W = cbind(w1_norm, w2_cont),
  lambda_seq = exp(seq(-1, -20, length = 500)),
  # arguments passed to hal9001::fit_hal()
  max_degree = NULL, num_knots = NULL, smoothness_orders = 0
)
pred_cont_haz <- predict(
  fit_cont, new_A = a_w2_cont, new_W = cbind(w1_norm, w2_cont)
)
emprisk_haz_cont <- mean(-log(pred_cont_haz))
if (FALSE) {
  library(sl3)
  cont_data <- as.data.table(list(W1 = w1_norm, W2 = w2_cont, A = a_w2_cont))
  cont_task <- sl3_Task$new(
    data = cont_data,
    covariates = c("W1", "W2"),
    outcome = "A"
  )

  # compare against GLM-based semiparametric strategy with constant variance
  hose_glm_lrnr <- Lrnr_density_semiparametric$new(
    mean_learner = Lrnr_glm_fast$new(family = gaussian())
  )
  fit_cont_hose <- hose_glm_lrnr$train(cont_task)
  pred_cont_hose <- fit_cont_hose$predict()
  emprisk_hose_cont <- mean(-log(pred_cont_hose))
  expect_lt(emprisk_haz_cont, emprisk_hose_cont)

  # compare against GLM-based semiparametric strategy with fitted variance
  hese_glm_lrnr <- Lrnr_density_semiparametric$new(
    mean_learner = Lrnr_glm_fast$new(family = gaussian()),
    var_learner = Lrnr_glm_fast$new(family = gaussian())
  )
  fit_cont_hese <- hese_glm_lrnr$train(cont_task)
  pred_cont_hese <- fit_cont_hese$predict()
  emprisk_hese_cont <- mean(-log(pred_cont_hese))
  expect_lt(emprisk_haz_cont, emprisk_hese_cont)
}
