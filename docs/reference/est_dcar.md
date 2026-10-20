# DCAR Projection of the Efficient Influence Function

DCAR Projection of the Efficient Influence Function

## Usage

``` r
est_dcar(
  psi_ipw,
  gn_pred_natural,
  gn_pred_shifted,
  Qn_pred_natural,
  Qn_pred_shifted
)
```

## Arguments

- psi_ipw:

  A `numeric` vector of IPW point estimates corresponding to the
  conditional density estimates of the exposure mechanism provided in
  the arguments `gn_pred_natural` and `gn_pred_shifted`.

- gn_pred_natural:

  A `matrix` of conditional density estimates of the exposure mechanism
  g(A\|W) along a grid of the regularization parameter, at the natural
  (i.e., observed) values of the exposure.

- gn_pred_shifted:

  A `matrix` of conditional density estimates of the exposure mechanism
  g(A+delta\|W) along a grid of the regularization parameter, at the
  shifted (i.e., counterfactual) values of the exposure.

- Qn_pred_natural:

  A `numeric` of the outcome mechanism estimate at the natural (i.e.,
  observed) values of the exposure. HAL regression is used for the
  estimate, with the regularization term chosen by cross-validation.

- Qn_pred_shifted:

  A `numeric` of the outcome mechanism estimate at the shifted (i.e.,
  counterfactual) values of the exposure. HAL regression is used for the
  estimate, with the regularization term chosen by cross-validation.
