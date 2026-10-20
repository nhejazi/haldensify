# IPW Estimator of the Causal Effects of Additive Modified Treatment Policies

IPW Estimator of the Causal Effects of Additive Modified Treatment
Policies

## Usage

``` r
ipw_shift(
  W,
  A,
  Y,
  delta = 0,
  n_bins = make_bins(A, "hist"),
  cv_folds = 10L,
  lambda_seq,
  ...,
  bin_type = c("equal_range", "equal_mass"),
  selector_type = c("dcar", "plateau", "gcv", "all")
)
```

## Arguments

- W:

  A `matrix`, `data.frame`, or similar containing a set of baseline
  covariates.

- A:

  A `numeric` vector corresponding to a exposure variable. The parameter
  of interest is defined as a location shift of this quantity.

- Y:

  A `numeric` vector of the observed outcomes.

- delta:

  A `numeric` value indicating the shift in the exposure to be used in
  defining the target parameter. This is defined with respect to the
  scale of the exposure (A).

- n_bins:

  A `numeric`, scalar or vector, indicating the number of bins into
  which the support of A is to be partitioned for constructing
  conditional density estimates.

- cv_folds:

  A `numeric` giving the number of folds to be used for
  cross-validation. Note that this form of sample splitting is used for
  the selection of tuning parameters by empirical risk minimization, not
  for the estimation of nuisance parameters (i.e., to relax regularity
  conditions).

- lambda_seq:

  A `numeric` sequence of the regularization parameter (L1 norm of HAL
  coefficients) to be used in fitting HAL models.

- ...:

  Additional arguments for model fitting to be passed directly to
  [`haldensify`](https://codex.nimahejazi.org/haldensify/reference/haldensify.md).

- bin_type:

  A `character` indicating the strategy to be used in creating bins
  along the observed support of `A`. For bins of equal range, use
  `"equal_range"`; to ensure each bin has the same number of
  observations, use instead `"equal_mass"`. For more information, see
  documentation of `grid_type` in
  [`haldensify`](https://codex.nimahejazi.org/haldensify/reference/haldensify.md).

- selector_type:

  A `character` indicating the selection strategy for identifying an
  efficent IPW estimator. The choices include `"gcv"` for global
  cross-validation, `"dcar"` for solving the EIF equation, and
  `"plateau"` for agnostic approaches (1) balancing changes in the IPW
  estimate and its standard error (adapting Lepski's method) and (2) a
  plateau detector for inflection points in the IPW estimator's
  trajectory. The option `"all"` runs all three selection strategies
  while sharing redundant computation between each.

## Examples

``` r
# simulate data
set.seed(11249)
n_obs <- 50
W1 <- rbinom(n_obs, 1, 0.6)
W2 <- rbinom(n_obs, 1, 0.2)
W3 <- rpois(n_obs, 3)
A <- rpois(n_obs, 3 * W1 - W2 + 2 * W1 * W2 + 4)
Y <- rbinom(n_obs, 1, plogis(A + W1 + W2 - W3 - W1 * W3))

# fit the IPW estimator
est_ipw <- ipw_shift(
  W = cbind(W1, W2, W3), A = A, Y = Y,
  delta = 0.5, cv_folds = 3L,
  n_bins = 4L, bin_type = "equal_range",
  lambda_seq = exp(seq(-1, -10, length = 100L)),
  # arguments passed to hal9001::fit_hal()
  max_degree = 1L,
  smoothness_orders = 0,
  reduce_basis = 1 / sqrt(n_obs)
)
#> 8% of observations outside training support...predictions trimmed.
```
