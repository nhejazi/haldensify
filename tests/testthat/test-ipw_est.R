library(here)
library(data.table)
set.seed(76924)

# load data simulation helpers
source(here("tests", "testthat", "helpers_dgp.R"))

# simulate test data + get truth for DGP


# estimate counterfactual mean
est_ipw_all <- ipw_shift(
  W = W, A = A, Y = Y, delta = delta, cv_folds = cv_folds,
  lambda_seq = exp(seq(-1, -17, length = 3000L)),
  ## arguments passed to hal9001::fit_hal()
  max_degree = NULL,
  smoothness_orders = 0,
  num_knots = NULL,
  reduce_basis = NULL, #1 / length(A),
  ## ...more undersmoothing args
  bin_type = bin_type,
  undersmooth_type = "all"
)
