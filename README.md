
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`haldensify`

<!-- badges: start -->

[![R-CMD-check](https://github.com/nhejazi/haldensify/workflows/R-CMD-check/badge.svg)](https://github.com/nhejazi/haldensify/actions)
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
<!-- badges: end -->

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
and Hejazi 2020). `haldensify` implements this conditional density
estimation strategy for use only with the highly adaptive lasso (HAL)
(Benkeser and van der Laan 2016; van der Laan 2017; van der Laan and
Benkeser 2018; Coyle et al. 2022; Hejazi, Coyle, and van der Laan 2020).
As the (generalized) propensity score is the primary ingredient in
inverse probability weighted (IPW) methods, `haldensify` builds loosely
on the advances of Ertefaie, Hejazi, and van der Laan (2022) to provide
nonparametric IPW estimators of the causal effects of continuous
treatments (Hejazi et al. 2022), which can be made to achieve the
non/semi-parametric efficiency bound by undersmoothing (lowering
regularization) over a family of the HAL conditional density estimators.

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
#> haldensify v0.2.2: Highly Adaptive Lasso Conditional Density Estimation
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
  max_degree = 3,
  reduce_basis = 1 / sqrt(n_train)
)
haldensify_fit
#> HAL Conditional Density Estimation
#> Number of bins over support of A: 10
#> CV-selected lambda: 0.0016
#> Summary of fitted HAL:
#> Warning in summary.hal9001(x$hal_fit): Coefficients for many lambda exist --
#> Summarizing coefficients corresponding to minimum lambda.
#>          coef                                    term
#>  1:  5.989688                             (Intercept)
#>  2: 10.498800                      [ I(bin_id >= 2) ]
#>  3: -9.673620                      [ I(W >= -3.353) ]
#>  4:  8.659440                      [ I(bin_id >= 6) ]
#>  5: -8.272041 [ I(bin_id >= 2) ] * [ I(W >= -2.371) ]
#>  6: -8.261273                      [ I(W >= -3.109) ]
#>  7:  8.054827                      [ I(bin_id >= 7) ]
#>  8:  8.013383                      [ I(bin_id >= 4) ]
#>  9:  8.001995                      [ I(bin_id >= 5) ]
#> 10: -7.649731                      [ I(W >= -2.157) ]
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
#> [1] 0.2818730 0.5513780 0.4449961 0.5329549 0.8722028 0.6150810
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
    @software{hejazi2021haldensify,
      author = {Hejazi, Nima S and Benkeser, David C and {van der Laan},
        Mark J},
      title = {{haldensify}: Highly adaptive lasso conditional density
        estimation},
      year  = {2021},
      howpublished = {\url{https://github.com/nhejazi/haldensify}},
      note = {{R} package version 0.2.0},
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

## Funding

The development of this software was supported in part through grants
from the National Library of Medicine (award number [T32
LM012417](https://reporter.nih.gov/project-details/9248418)), the
National Institute of Allergy and Infectious Diseases (award number [R01
AI074345](https://reporter.nih.gov/project-details/9926564)) of the
National Institutes of Health, and the National Science Foundation
(award number
[DMS 2102840](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2102840)).

-----

## License

© 2019-2022 [Nima S. Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    
    Copyright (c) 2019-2022 Nima S. Hejazi
    
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
Data Science and Advanced Analytics*, 2016:689. NIH Public Access.

</div>

<div id="ref-coyle2022hal9001-rpkg">

Coyle, Jeremy R, Nima S Hejazi, Rachael V Phillips, Lars WP van der
Laan, and Mark J van der Laan. 2022. *hal9001: The Scalable Highly
Adaptive Lasso*. <https://doi.org/10.5281/zenodo.3558313>.

</div>

<div id="ref-diaz2020causal">

Dı́az, Iván, and Nima S Hejazi. 2020. “Causal Mediation Analysis for
Stochastic Interventions.” *Journal of the Royal Statistical Society:
Series B (Statistical Methodology)* 82 (3): 661–83.
<https://doi.org/10.1111/rssb.12362>.

</div>

<div id="ref-diaz2011super">

Dı́az, Iván, and Mark J van der Laan. 2011. “Super Learner Based
Conditional Density Estimation with Application to Marginal Structural
Models.” *International Journal of Biostatistics* 7 (1): 1–20.

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

<div id="ref-ertefaie2022nonparametric">

Ertefaie, Ashkan, Nima S Hejazi, and Mark J van der Laan. 2022.
“Nonparametric Inverse Probability Weighted Estimators Based on the
Highly Adaptive Lasso.” *<span class="csl-no-emph">Revision Invited
at</span> Biometrics*. <https://arxiv.org/abs/2005.11303>.

</div>

<div id="ref-hejazi2022efficient">

Hejazi, Nima S, David C Benkeser, Iván Dı́az, and Mark J van der Laan.
2022. “Efficient Estimation of Modified Treatment Policy Effects Based
on the Generalized Propensity Score.”
*<span class="csl-no-emph">Forthcoming</span>*.

</div>

<div id="ref-hejazi2020hal9001-joss">

Hejazi, Nima S, Jeremy R Coyle, and Mark J van der Laan. 2020. “hal9001:
Scalable Highly Adaptive Lasso Regression in R.” *Journal of Open Source
Software*. <https://doi.org/10.21105/joss.02526>.

</div>

<div id="ref-vdl2017generally">

van der Laan, Mark J. 2017. “A Generally Efficient Targeted Minimum Loss
Based Estimator Based on the Highly Adaptive Lasso.” *International
Journal of Biostatistics* 13 (2).

</div>

<div id="ref-vdl2018highly">

van der Laan, Mark J, and David Benkeser. 2018. “Highly Adaptive Lasso
(HAL).” In *Targeted Learning in Data Science: Causal Inference for
Complex Longitudinal Studies*, 77–94. Springer.

</div>

</div>
