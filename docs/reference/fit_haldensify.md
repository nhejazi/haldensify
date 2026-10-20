# Fit Conditional Density Estimation over a Sequence of HAL Models

Fit Conditional Density Estimation over a Sequence of HAL Models

## Usage

``` r
fit_haldensify(
  A,
  W,
  wts = rep(1, length(A)),
  grid_type = "equal_range",
  n_bins = round(c(0.5, 1, 1.5, 2) * sqrt(length(A))),
  cv_folds = 5L,
  lambda_seq = exp(seq(-1, -13, length = 1000L)),
  smoothness_orders = 0L,
  ...
)
```

## Arguments

- A:

  The `numeric` vector of observed values.

- W:

  A `data.frame`, `matrix`, or similar giving the values of baseline
  covariates (potential confounders) for the observed units. These make
  up the conditioning set for the conditional density estimate.

- wts:

  A `numeric` vector of observation-level weights. The default is to
  weight all observations equally.

- grid_type:

  A `character` indicating the strategy to be used in creating bins
  along the observed support of `A`. For bins of equal range, use
  `"equal_range"`; consult the documentation of
  [`cut_interval`](https://ggplot2.tidyverse.org/reference/cut_interval.html)
  for more information. To ensure each bin has the same number of
  observations, use `"equal_mass"`; consult the documentation of
  [`cut_number`](https://ggplot2.tidyverse.org/reference/cut_interval.html)
  for details.

- n_bins:

  This `numeric` value indicates the number(s) of bins into which the
  support of `A` is to be divided. As with `grid_type`, multiple values
  may be specified, in which case cross-validation will be used to
  choose the optimal number of bins. The default sets the candidate
  choices of the number of bins based on heuristics tested in
  simulation.

- cv_folds:

  A `numeric` indicating the number of cross-validation folds to be used
  in fitting the sequence of HAL conditional density models.

- lambda_seq:

  A `numeric` sequence of values of the regularization parameter of
  Lasso regression; passed to
  [`fit_hal`](https://rdrr.io/pkg/hal9001/man/fit_hal.html).

- smoothness_orders:

  A `integer` indicating the smoothness of the HAL basis functions;
  passed to [`fit_hal`](https://rdrr.io/pkg/hal9001/man/fit_hal.html).
  The default is set to zero, for indicator basis functions.

- ...:

  Additional (optional) arguments of
  [`fit_hal`](https://rdrr.io/pkg/hal9001/man/fit_hal.html) that may be
  used to control fitting of the HAL regression model. Possible choices
  include `use_min`, `reduce_basis`, `return_lasso`, and
  `return_x_basis`, but this list is not exhaustive. Consult the
  documentation of
  [`fit_hal`](https://rdrr.io/pkg/hal9001/man/fit_hal.html) for complete
  details.

## Value

A `list`, containing density predictions for the sequence of fitted HAL
models; the index and value of the L1 regularization parameter
minimizing the density loss; and the sequence of empirical risks for the
sequence of fitted HAL models.

## Details

Estimation of the conditional density of A\|W via a cross-validated
highly adaptive lasso, used to estimate the conditional hazard of
failure in a given bin over the support of A.

## Examples

``` r
# simulate data: W ~ U[-4, 4] and A|W ~ N(mu = W, sd = 0.5)
set.seed(11249)
n_train <- 50
w <- runif(n_train, -4, 4)
a <- rnorm(n_train, w, 0.5)
# fit cross-validated HAL-based density estimator of A|W
haldensify_cvfit <- fit_haldensify(
  A = a, W = w, n_bins = 10L, lambda_seq = exp(seq(-1, -10, length = 100)),
  # the following arguments are passed to hal9001::fit_hal()
  max_degree = 3, reduce_basis = 1 / sqrt(length(a))
)
```
