# haldensify development notes

- [x] Add `...` args to
  [`haldensify()`](https://codex.nimahejazi.org/haldensify/reference/haldensify.md)
  to allow arbitrary arguments to be passed directly to `fit_hal()`.
- [x] Remove `use_future` argument to
  [`haldensify()`](https://codex.nimahejazi.org/haldensify/reference/haldensify.md),
  instead reducing to calling `future_mapply()`, with sequential
  evaluation via `plan(sequential)`.
- [x] Add example of working with bootstrap samples, e.g., using
  [`rsample`](https://rsample.tidymodels.org/reference/bootstraps.html).
- [x] Add a [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method to more easily visualize how empirical risk changes across the
  sequence of explored regularization parameter values.
- [x] Add an
  [`ipw_shift()`](https://codex.nimahejazi.org/haldensify/reference/ipw_shift.md)
  function for constructing IPW estimators of the mean counterfactual
  outcome of a stochastic shift intervention via `haldensify`.
- [ ] Add an argument to
  [`haldensify()`](https://codex.nimahejazi.org/haldensify/reference/haldensify.md)
  to allow normalization of the density estimates to improve estimation
  stability. Note that this normalized density is actually g(A\|W)/g(A),
  instead of the currently estimated g(A\|W).
