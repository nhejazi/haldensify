# haldensify 0.2.0

As of September 2021:
* Refinements of internal calls to `hal9001::fit_hal()` in keeping with updates
  to that package, for compatibility with its v0.4.0 CRAN release.
* The `smoothness_orders` argument of `hal9001::fit_hal()` previously was set
  through the `...` argument of `haldensify`; however, it has now been made a
  named argument to both `haldensify` and the internal `cv_haldensify` and
  `fit_haldensify` functions. The default is set to zero, for indicator basis
  functions, which differs from the default of `hal9001::fit_hal()` as of its
  v0.4.0 release.

# haldensify 0.1.5

As of April 2021:
* Changes to internal calls of `hal9001::fit_hal()` in order to correctly use
  the pared-down interface introduced in v0.4.0, contributed by @rachaelvp.
* The default for the grid of bins used for discretization of the variable `A`
  has been altered to be multiples of `sqrt(length(A))`.

# haldensify 0.1.0

As of April 2021:
* Updates to `haldensify` arguments (removal of `hal_max_degree` as a named
  argument) to simplify and better match use of `fit_hal` in `hal9001` v0.3.0+.
  This overhaul also included the addition of `...` arguments, now passed
  through `haldensify` and `fit_haldensify` to `cv_haldensify`, allowing all
  internal calls to `hal9001::fit_hal()` to specify the same arguments be
  passed for the fitting of HAL models.
* Changes to the default values of the argument `n_bins`, now setting this to
  (much) larger values that are themselves based on the sample size. This is in
  accordance with evidence from simulation experiments indicating that higher
  values of `n_bins` lead to significantly improved density estimates.
* Addition of argument `trim` and `trim_dens` to `predict.haldensify` to support
  the use of truncation more transparently. While the default was to set
  predictions for values of `new_A` outside the training support to zero, this
  has been changed to avoid trimming and, when the choice is made to trim the
  predictions, to set this value to `1/sqrt(length(new_A))`.
* Addition of a new method `print.haldensify` for a more user-friendly display
   of the prediction procedure's output, including the selected number of bins,
   the CV-selected choice of the regularization parameter, and the `summary`
   of the fitted HAL model.

# haldensify 0.0.9

As of February 2021:
* Addition of a method `plot.haldensify` to simplify visualizing the empirical
  risks of the sequence of HAL-based conditional density estimators across the
  grid of the regularization parameter, and necessary changes to the vignette.
* Preparation to add an option to visualize the conditional density estimates
  (of the estimator selected by cross-validation) via a `type` argument in the
  `plot.haldensify` method. Not yet implemented.
* Simplification of unit tests to remove unnecessary reliance on `dplyr`.
* Limit re-fitting of HAL model (after CV-selection of tuning parameters) in
  `haldensify()` to full-data fit by explicitly passing `n_folds = 1`.
* Avoid cross-validation procedure conditionally when the arguments `n_bin` and
  `grid_type` are fixed; add related assertion check in `predict()` when
  `haldensify()` skips cross-validation (since lambda selection skipped).
* Change how long-format repeated measures data is passed around in both
  `haldensify()` and `predict()` to clarify variable passing.
* Correct `predict()` method to truncate small conditional density estimates to
  a minimum value of [1 / sqrt(n)], based on the prediction set sample size.

# haldensify 0.0.8

As of January 2021:
* Addition of argument `hal_basis_list` to `haldensify()`, allowing for a HAL
  basis produced by `fit_hal()` to be passed into the HAL regression used for
  density estimation. This facilitates reduced computational overhead when
  requiring external cross-validation of nuisance functions (e.g., CV-TMLE) as
  well as working with bootstrap samples.
* Addition of argument `hal_max_degree` to `haldensify()`, allowing for control
  of the highest degree of interactions considered in the HAL model for density
  estimation. Like the above, this can reduce computational overhead.
* Fix a minor bug in `haldensify()` by passing `cv_folds` to the `n_folds`
  argument of `fit_hal()` when fitting HAL regression for density estimation.
  Previously, `cv_folds` was only used in constructing cross-validation (CV)
  folds for choosing tuning parameters, but the subsequent HAL regression was
  fiex to use the default number of folds specified in `fit_hal()` to choose
  the regularization parameter of the HAL regression for density estimation.
  Now, both CV to choose density estimation tuning parameters and CV to choose
  the lasso tuning parameter use the same number of folds.
* Addition of argument `...` to `haldensify()` so that arbitrary arguments can
  be passed to `fit_hal()` for density estimation, when not already specified
  as other arguments of the `haldensify()` constructor.
* Remove the unnecessary argument `use_future`, specifying parallel evaluation
  in a note instead.
* Add an option `"all"` to the `lambda_select` argument of the `predict()`
  method, allowing for predictions on the full (non-truncated) sequence of
  lambdas fitted on to be returned.
* Change truncation option in `predict()` method to 1/n instead of zero.

# haldensify 0.0.7

As of January 2021:
* Adds support to facilitate convenient marginal density estimation by creating
  automatically a constant vector when `W = NULL` is set in `haldensify()`.
* The `hal9001` dependency has been upgraded to v0.2.8 of that package, which
  introduced breaking changes in the names of slots in fitted model objects.
* The sequence of HAL models re-fit after identification of the regularization
  parameter selected by cross-validation has been padded with more aggressive
  choices of the parameter to ameliorate convergence issues in model fitting.
* Re-fitting of the HAL model with cross-validated choices of the number of
   bins, binning procedure, and regularization sequence has been altered to
   reuse the regularization sequence provided as input rather than subsetting
   the sequence to start with the cross-validation selector's choice of the
   parameter. Though convenient for undersmoothing `haldensify` estimates, this
   subsetting proved problematic for convergence of `glmnet()`.
* The `predict()` method's `cv_select` argument has been replaced in order to
  better facilitate undersmoothing. The new argument `lambda_select` defaults
  to the cross-validation selector but now easily allows access to the sequence
  of undersmoothed density estimates (less restrictive regularization values).
* The names of three slots in the `haldensify` S3 output class have been changed
  * `grid_type_tune_opt` is now `grid_type_cvselect`,
  * `n_bins_tune_opt` is now `n_bins_cvselect`, and
  * `cv_hal_fits_tune_opt` is now `cv_tuning_results`.

# haldensify 0.0.6

As of December 2020:
* Use of `plan(transparent)` has been changed to `plan(sequential)` based on
  ongoing development in the `future` package ecosystem.

As of June 2020:
* A short software paper for inclusion in JOSS has been added.

As of May 2020:
* The core cross-validation routine in `haldensify` for fitting HAL models has
  been slightly abstracted and moved to the new function `fit_haldensify`.
* The `haldensify()` wrapper function serves to cross-validate over choices of
  the histogram binning strategy and the number of bins.
* The defaults of `haldensify()` have been changed based on results of
  simulation experiments.
* The unnecessary argument `seed_int` in `haldensify()` has been removed.
* Fixes a bug introduced by returning predicted hazards as a vector instead of
  a matrix.
* An argument `cv_select`, defaulting to `TRUE`, has been added to the
  `predict` method, to make undersmoothing more accessible.
* A simple vignette has been added.

# haldensify 0.0.5

* First CRAN release.
* Documentation updates.

# haldensify 0.0.3

* First stable version of the package.
