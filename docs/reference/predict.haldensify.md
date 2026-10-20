# Prediction Method for HAL Conditional Density Estimation

Prediction Method for HAL Conditional Density Estimation

## Usage

``` r
# S3 method for class 'haldensify'
predict(
  object,
  ...,
  new_A,
  new_W,
  trim = TRUE,
  trim_min = NULL,
  lambda_select = c("cv", "undersmooth", "all")
)
```

## Arguments

- object:

  An object of class
  [`haldensify`](https://codex.nimahejazi.org/txshift/reference/haldensify.md),
  containing the results of fitting the highly adaptive lasso for
  conditional density estimation, as produced by a call to
  [`haldensify`](https://codex.nimahejazi.org/txshift/reference/haldensify.md).

- ...:

  Additional arguments passed to `predict` as necessary.

- new_A:

  The `numeric` vector or similar of the observed values for which a
  conditional density estimate is to be generated.

- new_W:

  A `data.frame`, `matrix`, or similar giving the values of baseline
  covariates (potential confounders) for the conditioning set of the
  observed values `A`.

- trim:

  A `logical` indicating whether estimates of the conditional density
  below the value indicated in `trim_min` should be truncated. The
  default value of `TRUE` enforces truncation of any values below the
  cutoff specified in `trim_min` and similarly truncates predictions for
  any of `new_A` falling outside of the training support.

- trim_min:

  A `numeric` indicating the minimum allowed value of the resultant
  density predictions. Any predicted density values below this tolerance
  threshold are set to the indicated minimum. The default is to use a
  scaled inverse square root of the sample size of the prediction set,
  i.e., 5/sqrt(n)/log(n) (another notable choice is 1/sqrt(n)). If there
  are observations in the prediction set with values of `new_A` outside
  of the support of the training set (i.e., provided in the argument `A`
  to
  [`haldensify`](https://codex.nimahejazi.org/txshift/reference/haldensify.md)),
  their predictions are similarly truncated.

- lambda_select:

  A `character` indicating whether to return the predicted density for
  the value of the regularization parameter chosen by the global
  cross-validation selector or whether to return an undersmoothed
  sequence (which starts with the cross-validation selector's choice but
  also includes all values in the sequence that are less restrictive).
  The default is `"cv"` for the global cross-validation selector.
  Setting the choice to `"undersmooth"` returns a matrix of predicted
  densities, with each column corresponding to a value of the
  regularization parameter less than or equal to the choice made by the
  global cross-validation selector. When `"all"` is set, predictions are
  returned for the full sequence of the regularization parameter on
  which the HAL model `object` was fitted.

## Value

A `numeric` vector of predicted conditional density values from a fitted
`haldensify` object.

## Details

Method for computing and extracting predictions of the conditional
density estimates based on the highly adaptive lasso estimator, returned
as an S3 object of class `haldensify` from
[`haldensify`](https://codex.nimahejazi.org/txshift/reference/haldensify.md).

## Examples

``` r
# simulate data: W ~ U[-4, 4] and A|W ~ N(mu = W, sd = 0.5)
n_train <- 50
w <- runif(n_train, -4, 4)
a <- rnorm(n_train, w, 0.5)
# HAL-based density estimator of A|W
haldensify_fit <- haldensify(
  A = a, W = w, n_bins = 10L, lambda_seq = exp(seq(-1, -10, length = 100)),
  # the following arguments are passed to hal9001::fit_hal()
  max_degree = 2, smoothness_orders = 0L, reduce_basis = 1 / sqrt(length(a))
)
# predictions to recover conditional density of A|W
new_a <- seq(-4, 4, by = 0.1)
new_w <- rep(0, length(new_a))
pred_dens <- predict(haldensify_fit, new_A = new_a, new_W = new_w)
```
