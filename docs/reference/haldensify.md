# Cross-validated HAL Conditional Density Estimation

Cross-validated HAL Conditional Density Estimation

## Usage

``` r
haldensify(
  A,
  W,
  wts = rep(1, length(A)),
  grid_type = "equal_range",
  n_bins = round(c(0.5, 1, 1.5, 2) * sqrt(length(A))),
  cv_folds = 5L,
  lambda_seq = exp(seq(-1, -13, length = 1000L)),
  smoothness_orders = 0L,
  hal_basis_list = NULL,
  ...
)
```

## Arguments

- A:

  The `numeric` vector observed values.

- W:

  A `data.frame`, `matrix`, or similar giving the values of baseline
  covariates (potential confounders) for the observed units. These make
  up the conditioning set for the density estimate. For estimation of a
  marginal density, specify a constant `numeric` vector or `NULL`.

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
  for details. The default is `"equal_range"` since this has been found
  to provide better performance in simulation experiments; however, both
  types may be specified (i.e., `c("equal_range", "equal_mass")`)
  together, in which case cross-validation will be used to select the
  optimal binning strategy.

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
  [`fit_hal`](https://rdrr.io/pkg/hal9001/man/fit_hal.html) via its
  argument `lambda`.

- smoothness_orders:

  A `integer` indicating the smoothness of the HAL basis functions;
  passed to [`fit_hal`](https://rdrr.io/pkg/hal9001/man/fit_hal.html).
  The default is set to zero, for indicator basis functions.

- hal_basis_list:

  A `list` consisting of a preconstructed set of HAL basis functions, as
  produced by [`fit_hal`](https://rdrr.io/pkg/hal9001/man/fit_hal.html).
  The default of `NULL` results in creating such a set of basis
  functions. When specified, this is passed directly to the HAL model
  fitted upon the augmented (repeated measures) data structure,
  resulting in a much lowered computational cost. This is useful, for
  example, in fitting HAL conditional density estimates with external
  cross-validation or bootstrap samples.

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

Object of class `haldensify`, containing a fitted `hal9001` object; a
vector of break points used in binning `A` over its support `W`; sizes
of the bins used in each fit; the tuning parameters selected by
cross-validation; the full sequence (in lambda) of HAL models for the
CV-selected number of bins and binning strategy; and the range of `A`.

## Details

Estimation of the conditional density A\|W through using the highly
adaptive lasso to estimate the conditional hazard of failure in a given
bin over the support of A. Cross-validation is used to select the
optimal value of the penalization parameters, based on minimization of
the weighted log-likelihood loss for a density.

## Note

Parallel evaluation of the cross-validation procedure to select tuning
parameters for density estimation may be invoked via the framework
exposed in the future ecosystem. Specifically, set
[`plan`](https://future.futureverse.org/reference/plan.html) for
[`future_mapply`](https://future.apply.futureverse.org/reference/future_mapply.html)
to be used internally.

## Examples

``` r
# simulate data: W ~ U[-4, 4] and A|W ~ N(mu = W, sd = 0.5)
set.seed(11249)
n_train <- 50
w <- runif(n_train, -4, 4)
a <- rnorm(n_train, w, 0.5)
# learn relationship A|W using HAL-based density estimation procedure
haldensify_fit <- haldensify(
  A = a, W = w, n_bins = 10L, lambda_seq = exp(seq(-1, -10, length = 100)),
  # the following arguments are passed to hal9001::fit_hal()
  max_degree = 3, reduce_basis = 1 / sqrt(length(a))
)
```
