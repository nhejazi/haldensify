# Agnostic IPW Estimator Selector via Lepski's and Variance-Blind Methods

Agnostic IPW Estimator Selector via Lepski's and Variance-Blind Methods

## Usage

``` r
plateau_selector(
  W,
  A,
  Y,
  delta = 0,
  gn_pred_natural,
  gn_pred_shifted,
  gn_fit_haldensify,
  Qn_pred_natural,
  Qn_pred_shifted,
  cv_folds = 10L,
  ci_level = 0.95,
  l1norm_mult = 10L
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

- gn_fit_haldensify:

  An object of class `haldensify` of the fitted conditional density
  model for the natural exposure mechanism. This should be the fit
  object returned by
  [`haldensify`](https://codex.nimahejazi.org/txshift/reference/haldensify.md)`[haldensify]`
  as part of a call to
  [`ipw_shift`](https://codex.nimahejazi.org/txshift/reference/ipw_shift.md).

- Qn_pred_natural:

  A `numeric` of the outcome mechanism estimate at the natural (i.e.,
  observed) values of the exposure. HAL regression is used for the
  estimate, with the regularization term chosen by cross-validation.

- Qn_pred_shifted:

  A `numeric` of the outcome mechanism estimate at the shifted (i.e.,
  counterfactual) values of the exposure. HAL regression is used for the
  estimate, with the regularization term chosen by cross-validation.

- cv_folds:

  A `numeric` giving the number of folds to be used for
  cross-validation. Note that this form of sample splitting is used for
  the selection of tuning parameters by empirical risk minimization, not
  for the estimation of nuisance parameters (i.e., to relax regularity
  conditions).

- ci_level:

  A `numeric` indicating the confidence level to be used in determining
  the cutoff used by the Lepski-type selector. This is only exposed for
  the sake of accommodating experimentation.

- l1norm_mult:

  A `numeric` indicating the multipler to be used by the plateau-based
  selector in reducing the candidate set of L1 norms relative to the
  choice made by the cross-validation selector.
