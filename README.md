
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`haldensify`

[![Travis-CI Build
Status](https://travis-ci.org/nhejazi/haldensify.svg?branch=master)](https://travis-ci.org/nhejazi/haldensify)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/nhejazi/haldensify?branch=master&svg=true)](https://ci.appveyor.com/project/nhejazi/haldensify)
[![Coverage
Status](https://img.shields.io/codecov/c/github/nhejazi/haldensify/master.svg)](https://codecov.io/github/nhejazi/haldensify?branch=master)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

> Nonparametric Conditional Density Estimation with the Highly Adaptive
> Lasso

**Author:** [Nima Hejazi](https://nimahejazi.org)

-----

## What’s `haldensify`?

The `haldensify` R package is designed to provide facilities for
nonparametric conditional density estimation based on the procedure
proposed by Díaz and van der Laan (2011). The core procedure involves
recovering conditional density estimates by performing pooled hazards
regressions on a long format data set built from the input data. For the
time being, `haldensify` is a minimal implementation of this density
estimation strategy, though future generalization of the core routines
may be possible.

-----

## Installation

Install the most recent *stable release* from GitHub via
[`devtools`](https://www.rstudio.com/products/rpackages/devtools/):

``` r
devtools::install_github("nhejazi/haldensify", build_vignettes = FALSE)
```

-----

## Example

A simple example illustrates how `haldensify` may be used to construct
conditional density estimates:

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
    @manual{hejazi2019haldensify,
      author = {Hejazi, Nima S},
      title = {haldensify: Nonparametric conditional density estimation
        with the highly adaptive lasso in {R}},
      year  = {2019},
      url = {https://github.com/nhejazi/haldensify},
      note = {R package version 0.0.1}
    }
```

-----

## Related

  - [R/`condensier`](https://github.com/osofr/condensier) - An R package
    providing an independent implementation of the same core
    methodology, though `condensier` allows for arbitrary selection of
    regression functions and a greater variety of hazard regression
    strategies, thus making it more general.

-----

## License

© 2019 [Nima S. Hejazi](https://nimahejazi.org)

The contents of this repository are distributed under the MIT license.
See below for details:

    MIT License
    
    Copyright (c) 2019 Nima S. Hejazi
    
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

<div id="ref-diaz2011super">

Díaz, Iván, and Mark J van der Laan. 2011. “Super Learner Based
Conditional Density Estimation with Application to Marginal Structural
Models.” *The International Journal of Biostatistics* 7 (1). De Gruyter:
1–20.

</div>

</div>
