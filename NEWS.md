# haldensify 0.0.8

As of January 2021:
* Addition of a `basis_list` argument to `haldensify()`, allowing for a HAL
  basis produced by `fit_hal()` to be passed in to the HAL regression for
  density estimation. This facilitates reduced computational overhead when
  requiring external cross-validation of nuisance functions (e.g., CV-TMLE) as
  well as working with bootstrap samples.
* [TO FILL IN]
* [TO FILL IN]
* [TO FILL IN]

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
