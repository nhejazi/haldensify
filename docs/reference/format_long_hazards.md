# Generate Augmented Repeated Measures Data for Pooled Hazards Regression

Generate Augmented Repeated Measures Data for Pooled Hazards Regression

## Usage

``` r
format_long_hazards(
  A,
  W,
  wts = rep(1, length(A)),
  grid_type = c("equal_range", "equal_mass"),
  n_bins = NULL,
  breaks = NULL
)
```

## Arguments

- A:

  The `numeric` vector or similar of the observed values of an
  intervention for a group of observational units of interest.

- W:

  A `data.frame`, `matrix`, or similar giving the values of baseline
  covariates (potential confounders) for the observed units whose
  observed intervention values are provided in the previous argument.

- wts:

  A `numeric` vector of observation-level weights. The default is to
  weight all observations equally.

- grid_type:

  A `character` indicating the strategy (or strategies) to be used in
  creating bins along the observed support of the intervention `A`. For
  bins of equal range, use "equal_range"; consult documentation of
  [`cut_interval`](https://ggplot2.tidyverse.org/reference/cut_interval.html)
  for more information. To ensure each bin has the same number of
  points, use "equal_mass"; consult documentation of
  [`cut_number`](https://ggplot2.tidyverse.org/reference/cut_interval.html)
  for details.

- n_bins:

  Only used if `grid_type` is set to `"equal_range"` or `"equal_mass"`.
  This `numeric` value indicates the number(s) of bins into which the
  support of `A` is to be divided.

- breaks:

  A `numeric` vector of break points to be used in dividing up the
  support of `A`. This is passed through the `...` argument to
  [`cut.default`](https://rdrr.io/r/base/cut.html) by
  [`cut_interval`](https://ggplot2.tidyverse.org/reference/cut_interval.html)
  or
  [`cut_number`](https://ggplot2.tidyverse.org/reference/cut_interval.html).

## Value

A `list` containing the break points used in dividing the support of `A`
into discrete bins, the length of each bin, and the reformatted,
"repeated measures" dataset. The reformatted dataset is a
[`data.table`](https://rdrr.io/pkg/data.table/man/data.table.html) of
repeated entries for observations up until the bin in which their `A`
falls, including an indicator for which bin an observation falls in, the
bin ID, observation ID, values of `W` for each observation, and,
possibly, observation-level weights.

## Details

Generates an augmented (long format, or repeated measures) dataset that
includes multiple records for each observation, a single record for each
discretized bin up to and including the bin in which a given observed
value of A falls. Such bins are derived from selecting break points over
the support of A. This repeated measures dataset is suitable for
estimating the hazard of failing in a particular bin over A using a
highly adaptive lasso (or other) classification model.
