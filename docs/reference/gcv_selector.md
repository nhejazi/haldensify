# IPW Estimator Selector via Global Cross-Validation

IPW Estimator Selector via Global Cross-Validation

## Usage

``` r
gcv_selector(
  W,
  A,
  Y,
  delta = 0,
  gn_pred_natural,
  gn_pred_shifted,
  Qn_pred_natural,
  Qn_pred_shifted
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

- gn_pred_natural:

  A `matrix` of conditional density estimates of the exposure mechanism
  g(A\|W) along a grid of the regularization parameter, at the natural
  (observed, actual) values of the exposure.

- gn_pred_shifted:

  A `matrix` of conditional density estimates of the exposure mechanism
  g(A+delta\|W) along a grid of the regularization parameter, at the
  shifted (counterfactual) values of the exposure.

- Qn_pred_natural:

  A `numeric` of the outcome mechanism estimate at the natural (i.e.,
  observed) values of the exposure. HAL regression is used for the
  estimate, with the regularization term chosen by cross-validation.

- Qn_pred_shifted:

  A `numeric` of the outcome mechanism estimate at the shifted (i.e.,
  counterfactual) values of the exposure. HAL regression is used for the
  estimate, with the regularization term chosen by cross-validation.
