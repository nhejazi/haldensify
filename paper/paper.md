---
title: "`haldensify`: Highly adaptive lasso conditional density estimation in `R`"
tags:
  - machine learning
  - causal inference
  - conditional density estimation
  - generalized propensity score
  - inverse probability weighting
  - semiparametric inference
  - R
authors:
  - name: Nima S. Hejazi
    orcid: 0000-0002-7127-2789
    affiliation: 1
  - name: Mark J. van der Laan
    orcid: 0000-0002-1019-8343
    affiliation: 2, 3
  - name: David Benkeser
    orcid: 0000-0002-1019-8343
    affiliation: 4
affiliations:
  - name: Department of Biostatistics, T.H. Chan School of Public Health, Harvard University
    index: 1
  - name: Division of Biostatistics, School of Public Health, University of California, Berkeley
    index: 2
  - name: Department of Statistics, University of California, Berkeley
    index: 3
  - name: Department of Biostatistics and Bioinformatics, Rollins School of Public Health, Emory University
    index: 4
date: 07 September 2022
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
in estimating the causal effects of continuous- or ordinal-valued treatments. In
such settings this covariate-conditional treatment density has been termed the
_generalized propensity score_ [@hirano2004propensity; @imai2004causal], and,
like its analog for binary treatments [@rosenbaum1983central], serves as a key
ingredient in developing both inverse probability weighted and doubly robust
estimators of causal effects [@diaz2012population; @haneuse2013estimation;
@diaz2018stochastic; @hejazi2022efficient].

# Statement of Need

Conditional density estimation is an important fundamental problem in the
computational sciences and statistics, having garnered (independent) attention
in machine learning [@takeuchi2009nonparametric; @sugiyama2012density],
semiparametric estimation [@qin1998inferences; @cheng2004semiparametric], and
causal inference [@hirano2004propensity; @diaz2011super; @zhu2015boosting].
Techniques for the nonparametric estimation of this quantity, complete with
asymptotic optimality guarantees, have received comparatively limited attention.
Similarly, despite the critical role of the generalized propensity score in the
estimation of the causal effects of continuous treatments, this nuisance
parameter is usually estimated with restrictive parametric modeling strategies,
ultimately sharply limiting the quality of downstream point estimates and
corresponding statistical inference (e.g., hypothesis tests, confidence
intervals). Approaches for flexibly estimating the generalized propensity score
have received limited attention [@diaz2011super; @zhu2015boosting], and software
implementations of these techniques are, to the best of our knowledge,
exceedingly rare, compared to, for example, regression algorithms for estimating
conditional means. `haldensify` aims to partially fill this gap by implementing
a flexible, nonparametric estimator of a conditional (or marginal) density,
appropriate for estimation of the generalized propensity score and useful for
the construction of inverse probability weighted or doubly robust estimators of
a class of causal effect parameters tailored to continuous treatments.

# Conditional Density Estimation and Modern Causal Inference

Conditional density estimation is a challenging and fundamental problem in
statistical learning theory. Owing to the high frequency with which conditional
density estimation arises in statistics and machine learning, a wide range of
techniques have been proposed -- under a correspondingly wide range of
assumptions. Some techniques are based in kernel smoothing [e.g.,
@takeuchi2009nonparametric], others in specialized neural network architectures
[e.g., @neuneier1994estimation], and others still in the direct estimation of
ratios of conditional densities [e.g., @sugiyama2012density]. Most approaches
make restrictive (parametric) assumptions about the form of the underlying
density functional or fail to achieve convergence rates (of the estimator to the
true, underlying conditional density) necessary for semiparametric inference. As
such, analysts must often negotiate a difficult tradeoff between tractability,
ease of implementation, and optimality properties of the chosen estimator. To
partially resolve this open challenge, `haldensify` implements a nonparametric
conditional density estimation procedure, making few assumptions regarding the
underlying form of the density functional, with convergence-rate guarantees
suitable for use in modern semiparametric inference and causal machine learning
applications.

The algorithm implemented in `haldensify` is an improved and tailored version of
the proposal of @diaz2011super, who formulated a nonparametric conditional
density estimator based on the relationship between the density and hazard
functions. This algorithm proceeds by, first, partitioning the support of the
dependent variable into a user-specified number of bins and recasting the input
dataset into a repeated measures structure, in which each observational unit is
represented by a variable number of records (with the last record corresponding
to the position of the bin over the discretized support into which the observed
value of the dependent variable falls). Next, the hazard probability,
conditional on any covariates, of the dependent variable falling in a given bin
along the discretized support is estimated by applying the highly adaptive lasso
(HAL) algorithm [@vdl2015generally; @benkeser2016highly; @vdl2017generally] (in
this case, for binary regression), via the `hal9001` package
[@hejazi2020hal9001-joss; @coyle2022hal9001-rpkg]; this step is often labeled
"pooled hazards" regression. Under plausible assumptions on the global variation
of the target functional, HAL has been shown to converge at a suitable rate
($\approx n^{-1/3}$ per @bibaut2019fast) for standard semiparametric efficiency
theory to apply to any estimators incorporating this conditional density
estimator; however, in this application, the $\ell_1$ (i.e., lasso) penalty of
the HAL estimator is updated to utilize a loss function suitable for density
estimation [@vdl2004asymptotic; @dudoit2005asymptotics]. In a final step, the
conditional hazard estimates are rescaled to conditional density estimates by
dividing the estimated hazard probabilities by the respective widths of the bins
along the support.

The advantages derived from the flexibility and rate-convergence properties of
this algorithm are especially apparent in causal inference problems with
continuous-valued treatments. In such problems, a key nuisance parameter is the
generalized propensity score (GPS), the conditional density of the treatment,
given covariates. This nuisance parameter is required to be well-estimated (in
a rate-convergence sense) for the construction of asymptotically efficient
estimators (e.g., of treatment effects), which attain the minimal possible
variance in a given regularity class. Such estimators are desirable since,
theoretically speaking, they are admit the tightest confidence intervals and
most sensitive hypothesis tests, making inference based upon these more
informative for downstream decision making. For example, the GPS is a nuisance
parameter required for the estimation of the counterfactual mean of a modified
treatment policy (MTP) [@haneuse2013estimation; @diaz2018stochastic], a type of
intervention that perturbs the natural (or observed) value of the treatment.
Doubly robust estimators of this causal effect are implemented in the `txshift`
`R` package [@hejazi2020txshift-joss; @hejazi2022txshift-rpkg], which relies
upon `haldensify` for estimation of the GPS and has been used in estimating
counterfactual vaccine efficacy based on MTPs interpretable as corresponding to
hypothetical (next-generation) vaccines that modulate the activity of target
immunologic biomarkers in vaccine efficacy clinical trials
[@hejazi2020efficient]. Alternative, asymptotically efficient and nonparametric
inverse probability weighted (IPW) estimators [@ertefaie2020nonparametric] of
such a counterfactual mean parameter are implemented in `haldensify`'s
`ipw_shift()` function, which constructs these IPW estimators by combining
`haldensify`'s GPS estimator (implemented in the eponymous `haldensify()`
function) with the sieve estimation framework to select an asymptotically
optimal IPW estimator with respect to a criterion rooted in semiparametric
efficiency theory; a formal description of these novel IPW estimators is given
in @hejazi2022efficient.

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
density estimation on the hazard scale. The nonparametric IPW estimators of
@hejazi2022efficient have been implemented in the `ipw_shift()` function.

In order to ensure a simplified and minimal API, the `haldensify` package
exposes only a limited set of functions: (1) the `haldensify()` function, which
facilitates the estimation of conditional (or marginal) densities as described
above, and (2) the `ipw_shift()` function, which implements nonparametric IPW
estimators of the causal effect of an additive modified treatment policy. As IPW
estimators require estimation of the generalized propensity score as an
intermediate (nuisance) step, this latter function internally calls the former;
moreover, the `ipw_shift()` function and the various selectors to which it
provides access (e.g., `selector_gcv()` for estimator selection based on
"global" cross-validation) have been studied from a theoretical-methodological
perspective in @hejazi2022efficient. The `haldensify()` function is complemented
by appropriate `predict()` and `plot()` methods, the former to allow for the
estimated conditional density to be evaluated at new values of the variable of
interest and its conditioning set and the latter to visualize the resultant
estimators along the regularization trajectory. The `ipw_shift()` function is
accompanied by a corresponding `confint()` method to easily generate confidence
intervals around the IPW point estimates. The S3 classes returned by both of
these functions have custom `print()` methods to allow for their results to be
easily inspected. Several internal utility functions, including, for example,
`cv_haldensify()`, `map_hazard_to_density()`, and `selector_dcar()`, implement
core aspects of the conditional density estimation and nonparametric IPW
estimation methodology.

# Availability

Future software development efforts will be focused primarily along two avenues:
(1) improving the computational aspects of the conditional density estimation
procedure, possibly to include random sampling from the internally generated
repeated measures dataset, and (2) further adjustments to the undersmoothing
estimator selection strategies implemented for nonparametric IPW estimation, to
be based on future methodological progress. Software maintenance efforts will
focus on ensuring that the package remains compatible with future versions of
the `hal9001` package [@coyle2022hal9001-rpkg; @hejazi2020hal9001-joss].
Currently, stable releases of the `haldensify` package are made available via
the Comprehensive `R` Archive Network [CRAN, @R] at
https://CRAN.R-project.org/package=haldensify, while development efforts are
carried out on the package's version-controlled repository, publicly hosted at
https://github.com/nhejazi/haldensify. To date (mid-September 2022), CRAN
records indicate that `haldensify` has been downloaded over 14,400 times.

# Acknowledgments

NSH's contributions to this work were supported in part by a grant from the
National Science Foundation (award number [DMS
2102840](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2102840)).

# References

