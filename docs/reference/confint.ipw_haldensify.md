# Confidence Intervals for IPW Estimates of the Causal Effects of Stochatic Shift Interventions

Confidence Intervals for IPW Estimates of the Causal Effects of
Stochatic Shift Interventions

## Usage

``` r
# S3 method for class 'ipw_haldensify'
confint(object, parm = seq_len(object$psi), level = 0.95, ...)
```

## Arguments

- object:

  An object of class `ipw_haldensify`, produced by invoking the function
  [`ipw_shift`](https://codex.nimahejazi.org/haldensify/reference/ipw_shift.md),
  for which a confidence interval is to be computed.

- parm:

  A `numeric` vector indicating indices of `object$est` for which to
  return confidence intervals.

- level:

  A `numeric` indicating the nominal level of the confidence interval to
  be computed.

- ...:

  Other arguments. Not currently used.

## Value

A named `numeric` vector containing the parameter estimate from a
`ipw_haldensify` object, alongside lower/upper Wald-style confidence
intervals at a specified coverage level.

## Details

Compute confidence intervals for estimates produced by
[`ipw_shift`](https://codex.nimahejazi.org/haldensify/reference/ipw_shift.md).

## Examples

``` r
# simulate data
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
  n_bins = 5L, bin_type = "equal_range",
  lambda_seq = exp(seq(-1, -10, length = 100L)),
  # arguments passed to hal9001::fit_hal()
  max_degree = 2,
  smoothness_orders = 0,
  reduce_basis = 1 / sqrt(n_obs)
)
#> 2% of observations outside training support...predictions trimmed.
confint(est_ipw)
#> # A tibble: 2 × 8
#>   lwr_ci   psi upr_ci se_est type     l1_norm lambda_idx gn_nbins
#>    <dbl> <dbl>  <dbl>  <dbl> <chr>      <dbl>      <dbl>    <int>
#> 1  0.601  0.74  0.843 0.0627 dcar_tol    6.84         39        5
#> 2  0.601  0.74  0.843 0.0627 dcar_min    6.84         39        5
```
