# Map Predicted Hazard to Predicted Density for a Single Observation

Map Predicted Hazard to Predicted Density for a Single Observation

## Usage

``` r
map_hazard_to_density(hazard_pred_single_obs)
```

## Arguments

- hazard_pred_single_obs:

  A `numeric` vector of predicted hazard of failure in a given bin
  (under a given partitioning of the support) for a single observational
  unit based on a long format data structure (from
  [`format_long_hazards`](https://codex.nimahejazi.org/txshift/reference/format_long_hazards.md)).
  This is the probability that a given value falls in a corresponding
  bin, given that it has not yet failed (fallen in a preceding bin), as
  per Díaz I, van der Laan MJ (2011). “Super learner based conditional
  density estimation with application to marginal structural models.”
  *International Journal of Biostatistics*, **7**(1), 1–20.
  [doi:10.2202/1557-4679.1356](https://doi.org/10.2202/1557-4679.1356) .
  .

## Value

A `matrix` composed of a single row and a number of columns specified by
the grid of penalization parameters used in fitting of the highly
adaptive lasso. This is the predicted conditional density for a given
observation, re-mapped from the hazard scale.

## Details

For a single observation, map a predicted hazard of failure (as an
occurrence in a particular bin, under a given partitioning of the
support) to a density.
