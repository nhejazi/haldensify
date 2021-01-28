# `haldensify` development

* Add example of working with bootstrap samples, e.g., using `rsample`.
* Add a `plot()` method to more easily visualize how empirical risk changes
  across the sequence of regularization parameter values.
* Remove `use_future` argument to `haldensify()`, instead reducing to calling
  `future_mapply()`, with sequential evaluation under `plan(sequential)`.
* Add a new argument to `haldensify()` to allow normalization of the density
  estimates (summing to 1) to improve estimation stability.
* Add an `ipw_est()` function for constructing IPW estimators of the mean
  counterfactual outcome of a stochastic shift intervention via `haldensify`.
