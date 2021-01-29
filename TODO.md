# `haldensify` development

- [x] Add `...` args to `haldensify()` to allow arbitrary arguments to be
  passed directly to `fit_hal()`.
- [x] Remove `use_future` argument to `haldensify()`, instead reducing to
  calling `future_mapply()`, with sequential evaluation via `plan(sequential)`.
- [x] Add example of working with bootstrap samples, e.g., using
  [`rsample`](https://rsample.tidymodels.org/reference/bootstraps.html).
- [ ] Add a `plot()` method to more easily visualize how empirical risk changes
  across the sequence of explored regularization parameter values.
- [ ] Add an `ipw_est()` function for constructing IPW estimators of the mean
  counterfactual outcome of a stochastic shift intervention via `haldensify`.
- [ ] Add an argument to `haldensify()` to allow normalization of the density
  estimates to improve estimation stability. Note that this normalized density
  is actually g(A|W)/g(A), instead of the currently estimated g(A|W).
