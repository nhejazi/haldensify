# haldensify 0.0.6

As of June 2020:
* A short software paper for inclusion in JOSS has been added.

As of May 2020:
* The core cross-validation routine in `haldensify` for fitting HAL models has
  been slightly abstracted and moved to the new function `fit_haldensify`.
* The `haldensify` wrapper function serves to cross-validate over choices of
  the histogram binning strategy and the number of bins.
* The defaults of `haldensify` have been changed based on results of simulation
  experiments; the unnecessary argument `seed_int` has also been removed.
* An argument `cv_select`, defaulting to `TRUE`, has been added to the
  `predict` method, to make undersmoothing more accessible.
* A simple vignette has been added.

# haldensify 0.0.5

* First CRAN release.
* Documentation updates.

# haldensify 0.0.3

* First stable version of package.
