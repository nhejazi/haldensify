# haldensify development

* Add a `plot()` method to more easily visualize how empirical risk changes
  across the sequence of regularization parameter values.
* Add example of working with bootstrap samples, e.g., using `rsample`.
* Add a `ipw_est()` function for constructing IPW estimators of the mean
  counterfactual outocme of a stochastic shift intervention via `haldensify`.
* Remove the `use_future` argument to `haldensify()`, instead reducing to just
  using `future_mapply()`, which should allow sequential evaluation under
  `plan(sequential)`.
