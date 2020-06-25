---
title: "`haldensify`: Highly adaptive lasso conditional density estimation in `R`"
tags:
  - machine learning
  - density estimation
  - causal inference
  - propensity score
  - R
authors:
  - name: Nima S. Hejazi
    orcid: 0000-0002-7127-2789
    affiliation: 1, 3
  - name: David C. Benkeser
    orcid: 0000-0002-1019-8343
    affiliation: 5
  - name: Mark J. van der Laan
    orcid: 0000-0002-1019-8343
    affiliation: 2, 3, 4
affiliations:
  - name: Graduate Group in Biostatistics, University of California, Berkeley
    index: 1
  - name: Division of Epidemiology & Biostatistics, School of Public Health University of California, Berkeley
    index: 2
  - name: Center for Computational Biology, University of California, Berkeley
    index: 3
  - name: Department of Statistics, University of California, Berkeley
    index: 4
  - name: Department of Biostatistics and Bioinformatics, Rollins School of Public Health, Emory University
    index: 5
date: 24 June 2020
bibliography: refs.bib
---

# Summary

The `haldensify` `R` package provides a toolbox for nonparametric conditional
density estimation based on the highly adaptive lasso, a flexible nonparametric
regression algorithm capable of reliably estimating functions falling in a large
class under only mild assumptions. Building on an earlier proposal by
@diaz2011super, the approach to conditional density estimation implemented in
the `haldensify` package leverages a histogram-like formulation that partitions
the support of the variable of interest into bins, recovering the conditional
density through a rescaling of estimated discretized conditional hazards. Such
conditional density estimates are useful, for example, in causal inference
problems in which the _generalized propensity score_ (for continuous-valued
exposures) [@hirano2004propensity; @imai2004causal; @zhu2015boosting], including
in evaluating causal effects defined by stochastic treatment regimes
[@diaz2012population; @diaz2018stochastic; @diaz2020causal;
@hejazi2020efficient].

# Background

Estimation of conditional density functions is a challenging problem in
statistical learning theory, for, unlike a regression problem in which one
estimates a conditional mean (e.g., $\mathbb{E}\{Y \mid X\}$), the goal is
instead to estimate the much more complex density function $f(Y \mid X)$.
A range of methods have been proposed, including kernel-based methods approaches
[e.g., @bashtannyk2001bandwidth; @takeuchi2009nonparametric], density ratio
estimation [e.g., @sugiyama2010conditional], particular neural network
architectures [e.g., @neuneier1994estimation], support vector machines
[@vapnik2000svm]; however, many such methods either place restrictive
assumptions on the form of the density function (e.g., location-scale families)
or are computationally expensive (e.g., neural network models). In order to
build conditional density estimators that could accommodate arbitrary ensemble
learning, @diaz2011super proposed an approach that proceeds in three simple
steps. First, the support of the variable of interest is partitioned into
a user-specified number of bins, with the dataset cast to a "long" format to
include repeated observations and an indicator of bin membership. Then, the
probability, conditional on covariates, of falling in a given bin (i.e., the
hazard) may be estimated by an arbitrary machine learning algorithm. Finally,
estimates of the conditional hazard could be rescaled to conditional density
estimates by dividing by the respective bin widths. Such an approach proved
a convenient and flexible strategy for conditional density estimation.

Recently, the highly adaptive lasso (HAL) regression algorithm was proposed as
a nonparametric approach for estimating functions that are càdlàg (right-hand
continuous with left-hand limits) and have bounded sectional variation norm
[@vdl2017generally; @vdl2017uniform]. HAL regression has been proven capable of
reliably estimating a wide variety of functions at a near-parametric
($n^{-1/3}$) rate [@bibaut2019fast]. A loss-based HAL estimator may be
constructed by using $\ell_1$-penalized (i.e., lasso) regression
[@tibshirani1996regression] to select
indicator basis functions that minimize the loss-specific empirical risk under
an appropriate loss function. Selection of the $\ell_1$ penalization parameter
may utilize a cross-validation selector [@vdl2003unified; @vdv2006oracle] or
alternative criteria that allow undersmoothing when HAL estimand is only
a nuisance function of the target parameter of interest [e.g.,
@vdl2019efficient; @ertefaie2020nonparametric]. The `hal9001` package
[@coyle2020hal9001], for the `R` language and environment for statistical
computing [@R], implements zeroth order HAL regression, which relies upon
indicator basis functions to construct a representation of the target function.

# `haldensify`

The `haldensify` `R` package combines the flexible binning and hazard estimation
strategy of @diaz2011super with highly adaptive lasso regression. In order to
estimate the conditional hazard of falling in a given bin over the support of
the variable of interest, HAL regression is peformed using the `hal9001` `R`
package [@coyle2020hal9001]. Cross-validation, using the `origami` `R` package
[@coyle2018origami], is implemented in order to select among several
hyperparameters: (1) the number of bins into which the variable of interest
should be discretized; (2) the strategy for bin creation (e.g., an equal number
of bins, or bins of equal mass); and (3) the $\ell_1$ penalization parameter for
HAL regression. Building upon the proposal of @diaz2011super, `haldensify`
additionally (1) adjusts the proposed algorithm to incorporate sample-level
weights; (2) replaces their use of an arbitrary classification model with the
HAL regression function; and (3) alters the HAL regression algorithm to use a
loss function tailored for hazard estimation, invoking $\ell_1$-penalization in
a manner consistent with this loss.

To provide a convenient and accessible interface, the `haldensify` package
requires the use of only the eponymous `haldensify()` function, to construct an
estimated conditional density function, and a `predict()` method, to evaluate
the estimated conditional density function at new values of the target variable
and conditioning set. When calling `haldensify()`, the arguments `n_bins` and
`cv_folds` can be adjusted to invoke a cross-validation selector to choose, by
empirical risk minimization, the optimal number of bins into which to partition
the support; moreover, the argument `lambda_seq` can be used to control the
range of $\ell_1$-penalization parameters defining the sequence of HAL models
explored for pooled estimation of the conditional hazard. Several internal
utility functions are also provided that may be of independent interest:
`fit_haldensify()` fits a sequence of HAL conditional density models with
cross-validation, `cv_haldensify()` fits HAL models for the conditional hazard
(mapping these to the density scale via `map_hazard_to_density()`),
`format_long_hazards()` creates an artifical repeated measures dataset by
partitioning the the target variable into a specified number of bins over its
support.

Future software development efforts are intended to focus primarily upon
improving computational aspects of the conditional density estimation procedure,
including random sampling rows from the artificial repeated measures dataset,
combining subset-specific HAL fits [@sapp2014subsemble], and performance
adjustments in tandem with the updates to the `hal9001` package
[@coyle2020hal9001]. Currently, stable releases of the `haldensify` package are
made available on the Comprehensive `R` Archive Network at
https://CRAN.R-project.org/package=haldensify, while both stable (branch
`master`) and development (branch `devel`) versions of the package are hosted at
https://github.com/nhejazi/haldensify.

# References

