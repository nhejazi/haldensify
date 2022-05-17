---
title: "`haldensify`: Highly adaptive lasso conditional density estimation in `R`"
tags:
  - machine learning
  - causal inference
  - semiparametric estimation
  - conditional density estimation
  - generalized propensity score
  - inverse probability weighting
  - R
authors:
  - name: Nima S. Hejazi
    orcid: 0000-0002-7127-2789
    affiliation: 1
  - name: David Benkeser
    orcid: 0000-0002-1019-8343
    affiliation: 4
  - name: Mark J. van der Laan
    orcid: 0000-0002-1019-8343
    affiliation: 2, 3
affiliations:
  - name: Division of Biostatistics, Department of Population Health Sciences, Weill Cornell Medicine, USA
    index: 1
  - name: Division of Biostatistics, School of Public Health, University of California, Berkeley, USA
    index: 2
  - name: Department of Statistics, University of California, Berkeley, USA
    index: 3
  - name: Department of Biostatistics and Bioinformatics, Rollins School of Public Health, Emory University, USA
    index: 4
date: 19 May 2022
bibliography: ../inst/REFERENCES.bib
---

# Summary

The `haldensify` `R` package serves as a toolbox for nonparametric conditional
density estimation based on the highly adaptive lasso, a flexible nonparametric
algorithm for the estimation of functional statistical parameters (e.g.,
conditional mean, hazard, density). Building upon an earlier proposal
[@diaz2011super], `haldensify` leverages the relationship between the hazard and
density functions to estimate the latter by applying pooled hazard regression to
a synthetic repeated measures dataset created from the input data, relying upon
the framework of cross-validated loss-based estimation to yield an optimal
estimator [@vdl2004asymptotic; @dudoit2005asymptotics]. While conditional
density estimation is a fundamental problem in statistics, arising naturally in
a variety of applications (including machine learning), it plays a critical role
in estimating the causal effects continuous- or ordinal-valued treatments. In
such settings this covariate-conditional treatment density has been termed the
_generalized propensity score_ [@hirano2004propensity; @imai2004causal], and,
like its analog for binary treatments [@rosenbaum1983central], serves as a key
ingredient in developing inverse probability weighted or doubly robust
estimators of causal effects [@diaz2012population; @haneuse2013estimation;
@diaz2018stochastic; @hejazi2020efficient].

# Statement of Need

Conditional density estimation is an important fundamental problem in the
computational sciences and statistics, having garnered (independent) attention
in machine learning [@takeuchi2009nonparametric; @sugiyama2010conditional;
@sugiyama2012density], semiparametric estimation [@qin1998inferences;
@cheng2004semiparametric], and causal inference [@hirano2004propensity;
@vdl2010targeted; @diaz2011super; @zhu2015boosting]. Techniques for the
nonparametric estimation of this quantity, complete with asymptotic optimality
guarantees, have received comparatively limited attention. Similarly, despite
the critical role of the generalized propensity score in the estimation of the
causal effects of continuous treatments, this nuisance parameter is usually
estimated with restrictive parametric modeling strategies, ultimately
compromising the quality of downstream estimates and corresponding inferences.
Approaches for flexibly estimating the generalized propensity score have
received limited attention [@diaz2011super; @zhu2015boosting], and software
implementations of these techniques are, to the best of our knowledge,
exceedingly rare. `haldensify` aims to resolve this need by implementing
a flexible, nonparametric estimator of a conditional (or marginal) density,
appropriate for estimation of the generalized propensity score and inverse
probability weighted (IPW) estimators of a class of causal effect parameters
applicable to continuous and ordinal treatments.

# Flexible Estimation of Conditional Density Functions

Conditional density estimation is a challenging problem in statistical learning
theory, significantly more complex than standard regression problems (i.e.,
conditional mean estimation). Given the ubiquity of this functional parameter,
it is unsurprising, then, that a range of techniques have been proposed,
including kernel-based approaches [e.g., @bashtannyk2001bandwidth;
@takeuchi2009nonparametric], conditional density ratio estimation [e.g.,
@sugiyama2010conditional; @sugiyama2012density], and specialized neural network
architectures [e.g., @neuneier1994estimation]. Yet, many of these approaches
place restrictive assumptions on the density function (e.g., $k$-times
differentiability), are computationally expensive or prohibitive (e.g., neural
network models), avoid estimating the quantity altogether (opting instead to
estimate a density ratio), or fail to achieve convergence rates necessary for
semiparametric estimation and inference.

@vdl2010targeted and @diaz2011super provided solutions to a subset of these
issues by proposing a conditional density estimation algorithm that proceeds in
three simple steps. First, the support of the variable of interest (e.g., the
treatment) is partitioned into a user-specified number of bins, casting the
dataset into a "long format" repeated measures structure, in which each unit is
represented by a number of records (the number of records matches the position
of the bin into which the observed value variable of interest falls). Then, the
probability, conditional on covariates, of the variable of interest falling in
a given bin along its discretized support (i.e., the hazard) is estimated by
applying any machine learning algorithm for binary regression to this repeated
measures data structure. This is merely an application of pooled hazards
regression. Finally, the conditional hazard estimates are translated to the
conditional density scale by dividing by the estimated hazard probabilities by
the respective bin widths. This approach is made particularly flexible by the
fact that it may be readily applied to arbitrary data structures and
straightforwardly incorporates machine learning at the hazard estimation step.
@diaz2011super advocated for the use of the Super Learner algorithm
[@vdl2007super], a theoretically principled version of ensemble machine learning
or model stacking [@wolpert1992stacked; @breiman1996stacked].

Although Super Learning is both highly flexible and compatible with the
framework of loss-based estimation, it is generally incapable of yielding
estimators with desirable rate-convergence properties. The highly adaptive lasso
(HAL) algorithm was developed to a nonparametric approach for estimating
functions that are càdlàg (right-hand continuous with left-hand limits) and have
bounded sectional variation norm [@vdl2017generally; @vdl2017uniform].
Importantly, the HAL estimation procedure has been proven to reliably estimate
target functionals at a suitable rate ($\approx n^{-1/3}$, per @bibaut2019fast)
for standard semiparametric theory to be applied to resultant estimators.
A loss-based HAL estimator may be constructed by using $\ell_1$-penalized (i.e.,
lasso) regression [@tibshirani1996regression] to select indicator basis
functions that minimize the empirical risk with respect to an appropriately
chosen loss function. Selection of the $\ell_1$ penalization term be based on
a cross-validation selector [@vdl2003unified; @vdv2006oracle] or alternative
selection criteria compatible with sieve estimation (via undersmoothing). This
latter class of selection criteria are appropriate when HAL is used only to
estimate nuisance function of a (usually low-dimensional) target parameter. In
such cases, undersmoothing has been shown to allow for the formulation of
regular, asymptotically linear and efficient IPW estimators of causal effects
[@ertefaie2020nonparametric; @hejazi2022efficient].

# Nonparametric IPW Estimation with the Generalized Propensity Score

TODO

# `haldensify`'s Scope

The `haldensify` `R` package combines the binning and hazard estimation strategy
of @diaz2011super with HAL regression [@benkeser2016highly], resulting in
a flexible, nonparametric conditional density estimator. This procedure --
accessible via the eponymous `haldensify()` function -- relies upon the
`hal9001` `R` package [@coyle2022hal9001-rpkg; @hejazi2020hal9001-joss] for the
HAL regression step and upon the `origami` `R` package [@coyle2018origami] for
cross-validated selection of tuning parameters (e.g., number of bins, $\ell_1$
regularization) so as to empirically minimize the negative log-density loss
[@dudoit2005asymptotics]. `haldensify` additionally adjusts the proposal of
@diaz2011super to (1) incorporate sample-level weights and (2) apply HAL
regression to the repeated measures data structure in a manner tailored for
estimation on the hazard scale.

The `haldensify` package exposes only a limited set of functions in order to
ensure a simple API: (1) the `haldensify()` function, which facilitates the
estimation of conditional (or marginal) densities as described above, and (2)
the `ipw_shift()` function, which implements nonparametric IPW estimators of the
causal effect of an additive modified treatment policy. As IPW estimation
requires estimation of the generalized propensity score as an intermediate step,
this latter function internally calls the former; moreover, the `ipw_shift()`
function and the various selectors to which it provides access (e.g.,
`selector_gcv()` for estimator selection based on "global" cross-validation)
have been studied from a methodological perspective in @hejazi2022efficient. The
`haldensify()` function is complemented by an appropriate `predict()` method, to
allow for the estimated conditional density to be evaluated at new values of the
variable of interest and its conditioning set, while the `ipw_shift()` function
is accompanied by a corresponding `confint()` method to generate confidence
intervals for the IPW estimates; the custom S3 classes returned by each function
have custom `print()` methods to allow for the results to be easily evaluated.
Several internal utility functions, including, for example, `cv_haldensify()`,
`map_hazard_to_density()`, and `selector_dcar()` implement core aspects of the
conditional density estimation and nonparametric IPW estimation methodology.

# Availability

Future software development efforts will be focused primarily along two avenues:
(1) improving the computational aspects of the conditional density estimation
procedure, possibly to include random sampling from the synthetic repeated
measures dataset created internally, and (2) further adjustments to the various
estimator selection strategies implemented for IPW estimation, to be based on
future methodological progress. Software maintenance efforts will focus on
ensuring that the package remains compatible with future versions of the
`hal9001` package [@coyle2022hal9001-rpkg; @hejazi2020hal9001-joss]. Currently,
stable releases of the `haldensify` package are made available via the
Comprehensive `R` Archive Network [CRAN, @R] at
https://CRAN.R-project.org/package=haldensify, while development efforts are
carried out on the package's version-controlled repository, publicly hosted at
https://github.com/nhejazi/haldensify. To date (mid-May 2022), CRAN records
indicate that `haldensify` has been downloaded over 13,000 times.

# Acknowledgments

NSH's contributions to this work were supported in part by a grant from the
National Science Foundation (award number [DMS
2102840](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2102840)).

# References

