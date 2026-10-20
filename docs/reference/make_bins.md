# Histogram Binning Procedures for Pooled Hazards Regression

Histogram Binning Procedures for Pooled Hazards Regression

## Usage

``` r
make_bins(grid_var, grid_type = c("hist", "scaled"), max_bins = 30L)
```

## Arguments

- grid_var:

  The `numeric` vector over which histogram-based binning is to be
  performed.

- grid_type:

  A `character` indicating the choice of binning rule, with `"hist"`
  corresponding to the use of several rules proposed for optimal
  histogram construction and `"scaled"` corresponding to the use of
  various pre-set multiples of the square root of the sample size.

- max_bins:

  A `numeric` indicating the maximum number of bins that are allowed in
  the grid for building the histogram based discretization.
