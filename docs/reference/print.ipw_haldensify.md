# Print: IPW Estimates of the Causal Effects of Stochatic Shift Interventions

Print: IPW Estimates of the Causal Effects of Stochatic Shift
Interventions

## Usage

``` r
# S3 method for class 'ipw_haldensify'
print(x, ..., ci_level = 0.95)
```

## Arguments

- x:

  An object of class `ipw_haldensify`.

- ...:

  Other options (not currently used).

- ci_level:

  A `numeric` indicating the level of the confidence interval to be
  computed.

## Value

None. Called for the side effect of printing an informative summary of
slots of objects of class `ipw_haldensify`.

## Details

The `print` method for objects of class `ipw_haldensify`

## Examples

``` r
# simulate data
set.seed(11249)
n_obs <- 50
W1 <- rbinom(n_obs, 1, 0.6)
W2 <- rbinom(n_obs, 1, 0.2)
A <- rnorm(n_obs, (2 * W1 - W2 - W1 * W2), 2)
Y <- rbinom(n_obs, 1, plogis(3 * A + W1 + W2 - W1 * W2))

# fit the IPW estimator
est_ipw_shift <- ipw_shift(
  W = cbind(W1, W2), A = A, Y = Y,
  delta = 0.5, n_bins = 3L, cv_folds = 3L,
  lambda_seq = exp(seq(-1, -10, length = 100L)),
  # arguments passed to hal9001::fit_hal()
  max_degree = 1,
  # ...continue arguments for IPW
  selector_type = "gcv"
)
#> 2% of observations outside training support...predictions trimmed.
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
print(est_ipw_shift)
#> Counterfactual Mean of Shifted Treatment
#> Intervention: Treatment + 0.5
#> IPW Estimator Criterion: Global CV
#> Estimate: 0.7018
#> Std. Error: 0.0749
#> 95% CI: [0.5386, 0.8259]
#> EIF Mean: -0.0497
```
