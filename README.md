
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`haldensify`

[![Travis-CI Build
Status](https://travis-ci.com/nhejazi/haldensify.svg?branch=master)](https://travis-ci.com/nhejazi/haldensify)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/nhejazi/haldensify?branch=master&svg=true)](https://ci.appveyor.com/project/nhejazi/haldensify)
[![Coverage
Status](https://img.shields.io/codecov/c/github/nhejazi/haldensify/master.svg)](https://codecov.io/github/nhejazi/haldensify?branch=master)
[![CRAN](http://www.r-pkg.org/badges/version/haldensify)](https://www.r-pkg.org/pkg/haldensify)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/haldensify)](https://CRAN.R-project.org/package=haldensify)
[![CRAN total
downloads](http://cranlogs.r-pkg.org/badges/grand-total/haldensify)](https://CRAN.R-project.org/package=haldensify)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3698329.svg)](https://doi.org/10.5281/zenodo.3698329)

> Highly Adaptive Lasso Conditional Density Estimation

**Authors:** [Nima Hejazi](https://nimahejazi.org), [David
Benkeser](https://www.sph.emory.edu/faculty/profile/#!dbenkes), and
[Mark van der Laan](https://vanderlaan-lab.org/about/)

-----

## What’s `haldensify`?

The `haldensify` R package is designed to provide facilities for
nonparametric conditional density estimation based on a flexible
procedure proposed initially by Dı́az and van der Laan (2011). The core
of the implemented methodology involves recovering conditional density
estimates by performing pooled hazards regressions so as to assess the
conditional hazard that an observed value falls in a given bin over the
(conditional) support of the variable of interest. Such conditional
density estimates are useful, for example, in causal inference problems
in which the *generalized propensity score* (for continuous-valued
exposures) must be estimated (Dı́az and van der Laan 2012, 2018; Dı́az
and Hejazi 2020). `haldensify` implements this condtional density
estimation strategy specifically for use only with the highly adaptive
lasso (Benkeser and van der Laan 2016; van der Laan 2017; van der Laan
and Benkeser 2018; Coyle, Hejazi, and van der Laan 2020; Hejazi, Coyle,
and van der Laan 2020).

-----

## Installation

For standard use, we recommend installing the package from
[CRAN](https://CRAN.R-project.org/package=haldensify) via

``` r
install.packages("haldensify")
```

To contribute, install the *development version* of `haldensify` from
GitHub via [`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("nhejazi/haldensify")
```

-----

## Example

A simple example illustrates how `haldensify` may be used to train a
highly adaptive lasso model to obtain conditional density estimates:

``` r
library(haldensify)
#> haldensify v0.1.5: Highly Adaptive Lasso Conditional Density Estimation
set.seed(76924)

# simulate data: W ~ U[-4, 4] and A|W ~ N(mu = W, sd = 0.25)
n_train <- 100
w <- runif(n_train, -4, 4)
a <- rnorm(n_train, w, 0.25)

# HAL-based density estimate of A|W
haldensify_fit <- haldensify(
  A = a, W = w,
  n_bins = 10, grid_type = "equal_range",
  lambda_seq = exp(seq(-1, -10, length = 100)),
  # arguments passed to hal9001::fit_hal()
  max_degree = 5, smoothness_orders = 0,
  num_knots = NULL, reduce_basis = 0.05
)
haldensify_fit
#> HAL Conditional Density Estimation
#> Number of bins over support of A: 10
#> CV-selected lambda: 0.0016
#> Summary of fitted HAL:
#> Warning in summary.hal9001(x$hal_fit): Coefficients for many lamdba exist --
#> Summarizing coefficients corresponding to minimum lambda.
#>          coef               term
#>  1:  5.977489        (Intercept)
#>  2: 10.481606 [ I(bin_id >= 1) ]
#>  3: 10.440778 [ I(W >= -2.371) ]
#>  4: -9.663779 [ I(W >= -3.546) ]
#>  5:  8.965791 [ I(bin_id >= 5) ]
#>  6:  8.621709 [ I(bin_id >= 6) ]
#>  7:  8.621184 [ I(bin_id >= 4) ]
#>  8:  8.299777 [ I(bin_id >= 8) ]
#>  9: -8.253294  [ I(W >= -3.12) ]
#> 10:  8.091661 [ I(bin_id >= 3) ]
```

We can also visualize the empirical risk (with respect to density loss)
in terms of the solution path of the lasso regularization parameter:

``` r
# just use the built-in plot method
plot(haldensify_fit)
```

<img src="README-example-plot-1.png" width="80%" />

Finally, we can obtain conditional density estimates from the trained
model on the training (or on new) data:

``` r
# use the built-in predict method to get predictions
pred_haldensify <- predict(haldensify_fit, new_A = a, new_W = w)
head(pred_haldensify)
#> [1] 0.8677402 0.4276165 0.4430710 0.5334161 0.8721339 0.6149775
```

For more details, check out the [package
vignette](https://code.nimahejazi.org/haldensify/articles/intro_haldensify)
on the corresponding `pkgdown` site.

-----

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/nhejazi/haldensify/issues).

-----

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/nhejazi/haldensify/blob/master/CONTRIBUTING.md)
prior to submitting a pull request.

-----

## Citation

After using the `haldensify` R package, please cite the following:

``` 
    @software{hejazi2020haldensify,
      author = {Hejazi, Nima S and Benkeser, David C and {van der Laan},
        Mark J},
      title = {{haldensify}: Highly adaptive lasso conditional density
        estimation },
      year  = {2020},
      howpublished = {\url{https://github.com/nhejazi/haldensify}},
      note = {{R} package version 0.0.5},
      url = {https://doi.org/10.5281/zenodo.3698329},
      doi = {10.5281/zenodo.3698329}
    }
```

-----

## Related

  - [R/`hal9001`](https://github.com/tlverse/hal9001) – The highly
    adaptive lasso estimator used internally to constructed conditional
    density estimates.

-----

## License

© 2019-2021 [Nima S. Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    
    Copyright (c) 2019-2021 Nima S. Hejazi
    
    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

-----

## References

<div id="refs" class="references">

<div id="ref-benkeser2016highly">

Benkeser, David, and Mark J van der Laan. 2016. “The Highly Adaptive
Lasso Estimator.” In *Proceedings of the International Conference on
Data Science and Advanced Analytics. IEEE International Conference on
Data Science and Advanced Analytics*, 2016:689. NIH Public Access.

</div>

<div id="ref-coyle2020hal9001-rpkg">

Coyle, Jeremy R, Nima S Hejazi, and Mark J van der Laan. 2020. *hal9001:
The Scalable Highly Adaptive Lasso*.
<https://doi.org/10.5281/zenodo.3558313>.

</div>

<div id="ref-diaz2020causal">

Dı́az, Iván, and Nima S Hejazi. 2020. “Causal Mediation Analysis for
Stochastic Interventions.” *Journal of the Royal Statistical Society:
Series B (Statistical Methodology)*.
<https://doi.org/10.1111/rssb.12362>.

</div>

<div id="ref-diaz2011super">

Dı́az, Iván, and Mark J van der Laan. 2011. “Super Learner Based
Conditional Density Estimation with Application to Marginal Structural
Models.” *The International Journal of Biostatistics* 7 (1): 1–20.

</div>

<div id="ref-diaz2012population">

———. 2012. “Population Intervention Causal Effects Based on Stochastic
Interventions.” *Biometrics* 68 (2): 541–49.

</div>

<div id="ref-diaz2018stochastic">

———. 2018. “Stochastic Treatment Regimes.” In *Targeted Learning in Data
Science: Causal Inference for Complex Longitudinal Studies*, 167–80.
Springer Science & Business Media.

</div>

<div id="ref-hejazi2020hal9001-joss">

Hejazi, Nima S, Jeremy R Coyle, and Mark J van der Laan. 2020. “hal9001:
Scalable Highly Adaptive Lasso Regression in R.” *Journal of Open Source
Software*. <https://doi.org/10.21105/joss.02526>.

</div>

<div id="ref-vdl2017generally">

van der Laan, Mark J. 2017. “A Generally Efficient Targeted Minimum Loss
Based Estimator Based on the Highly Adaptive Lasso.” *The International
Journal of Biostatistics* 13 (2).

</div>

<div id="ref-vdl2018highly">

van der Laan, Mark J, and David Benkeser. 2018. “Highly Adaptive Lasso
(HAL).” In *Targeted Learning in Data Science*, 77–94. Springer.

</div>

</div>
