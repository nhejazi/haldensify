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
date: 18 May 2022
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
in machine learning [@takeuchi2009nonparametric; @sugiyama2010conditional;
@sugiyama2012density], semiparametric estimation [@qin1998inferences;
@cheng2004semiparametric], and causal inference [@hirano2004propensity;
@vdl2010targeted; @diaz2011super; @zhu2015boosting]. Techniques for the
nonparametric estimation of this quantity, complete with asymptotic optimality
guarantees, have received comparatively limited attention. Similarly, despite
the critical role of the generalized propensity score in the estimation of the
causal effects of continuous treatments, this nuisance parameter is usually
estimated with restrictive parametric modeling strategies, ultimately sharply
limiting the quality of downstream point estimates and corresponding statistical
inference (e.g., hypothesis tests, confidence intervals). Approaches for
flexibly estimating the generalized propensity score have received limited
attention [@diaz2011super; @zhu2015boosting], and software implementations of
these techniques are, to the best of our knowledge, exceedingly rare, compared
to, for example, regression algorithms for estimating conditional means.
`haldensify` aims to partially fill this gap by implementing a flexible,
nonparametric estimator of a conditional (or marginal) density, appropriate for
estimation of the generalized propensity score and useful for the construction
of inverse probability weighted or doubly robust estimators of a class of causal
effect parameters tailored to continuous treatments.

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
network models), avoid estimating the target quantity altogether (opting instead
to estimate a density ratio), or fail to achieve convergence rates necessary for
semiparametric estimation and inference.

Both @vdl2010targeted and @diaz2011super provided solutions to a subset of these
issues by proposing a conditional density estimation algorithm that proceeds in
three simple steps. First, the support of the variable of interest (e.g., the
treatment) is partitioned into a user-specified number of bins, casting the
dataset into a "long format" repeated measures structure, in which each unit is
represented by a number of records (with the number of records matching the
position of the bin into which the observed value of the variable of interest
falls). Then, the hazard probability, conditional on any covariates, of the
variable of interest falling in a given bin along its discretized support is
estimated by applying any machine learning algorithm for binary regression to
this repeated measures data structure; this is merely an application of pooled
hazard regression. Finally, the conditional hazard estimates are translated to
the conditional density scale by dividing the estimated hazard probabilities by
the respective bin widths. This approach is made particularly flexible by the
fact that it may be readily applied to arbitrary data structures and
straightforwardly incorporates machine learning at the hazard estimation step
(via the binary regression procedure). @diaz2011super advocated for the use of
the Super Learner algorithm [@vdl2007super], a theoretically principled variant
of ensemble machine learning or model stacking [@wolpert1992stacked;
@breiman1996stacked] for this regression step.

Although Super Learning is both highly flexible and compatible with the
framework of loss-based estimation, it is generally incapable of yielding
nuisance estimators with desirable rate-convergence properties for
semiparametric estimation of low-dimensional statistical parameters. The highly
adaptive lasso (HAL) algorithm was developed as a nonparametric approach for
estimating functions that are càdlàg (right-hand continuous with left-hand
limits) and have bounded sectional variation norm [@vdl2017generally;
@vdl2017uniform]. Importantly, the HAL estimation procedure has been proven to
reliably estimate target functionals at a suitable rate (per @bibaut2019fast,
$\approx n^{-1/3}$) for standard semiparametric theory to remain applicable to
any downstream estimators. A loss-based HAL estimator may be constructed by
using $\ell_1$-penalized (i.e., lasso) regression [@tibshirani1996regression] to
select basis functions that minimize the estimator's empirical risk with respect
to an appropriately chosen loss function. Selection of the $\ell_1$ penalization
term can be based on cross-validation principles [@vdl2004asymptotic;
@dudoit2005asymptotics] or alternative selection criteria compatible with sieve
estimation techniques (e.g., undersmoothing). This latter class of selection
criteria are appropriate when HAL is leveraged to estimate the nuisance
functions of a (usually low-dimensional) target parameter, common in the
semiparametric estimation of quantities motivated by causal inference
[@bang2005doubly; @vdl2006targeted; @vdl2011targeted]. In such cases,
undersmoothing methods have been developed to allow for the formulation of
regular, asymptotically linear and efficient inverse probability weighted (IPW)
estimators of causal effects [e.g., @ertefaie2020nonparametric;
@hejazi2022efficient].

# Nonparametric IPW Estimation with the Generalized Propensity Score

A popular framework for defining and evaluating the causal effects of continuous
treatments is that of modified treatment policies [@haneuse2013estimation;
@diaz2018stochastic; @hejazi2022efficient], which define interventions that
shift (or modify) the natural value of the treatment. For example, in a setting
with a continuous treatment $A$, in which we additionally collect baseline
covariates $W$ and measure an outcome $Y$, we could consider a hypothetical
intervention setting the value of $A$ via $d(A,W; \delta) = A + \delta(W)$, for
a user-defined function $d(A,W;\delta)$ indexed by a scalar $\delta(W)$. This
intervention regime is a simple example of a modified treatment policy (MTP); it
can be thought of as mapping the observed $A$ to a counterfactual $A_{\delta}$
that is itself an additive shift of the natural value of $A$. The counterfactual
mean of such an intervention may be expressed $\mathbb{E}[Y(A_{\delta})]$, where
$Y(A_{\delta})$ is the potential outcome that would have been observed had the
treatment taken the value $A_{\delta}$. Both @haneuse2013estimation and
@diaz2018stochastic proposed classical substitution (G-computation), inverse
probability weighted, and doubly robust estimators of a statistical functional
$\psi$ that identifies this counterfactual mean under standard assumptions.
Doubly robust estimators of $\psi$ have been implemented in the `txshift` `R`
package [@hejazi2020txshift-joss; @hejazi2022txshift-rpkg], used in the
estimation of the causal effects of MTPs in vaccine efficacy trials
[@hejazi2020efficient]; such estimation frameworks are usually necessary in
order to take advantage of flexible estimators of nuisance parameters.

Despite the popularity of doubly robust estimation procedures, IPW estimators
can be modified to accommodate data adaptive estimation of the (generalized)
propensity score. Such nonparametric IPW estimators, based on HAL, have been
studied by @ertefaie2020nonparametric in the context of binary treatments and by
@hejazi2022efficient for continuous treatments. In the case of MTPs, the IPW
estimator of $\psi$ is $\psi_{n}^{\text{IPW}} = n^{-1} \sum_{i=1}^n
\{\tilde{g}_{n,A}(A_i \mid W_i) / g_{n,A}(A_i \mid W_i)\} Y_i$, where $g_{n,A}$
is a generalized propensity score estimator (e.g., produced by `haldensify()`)
and $\tilde{g}_{n,A}$ is this same quantity evaluated at the post-intervention
value of the treatment $A_{\delta}$. Usually, $g_{n,A}$ must be estimated via
parametric modeling strategies in order for $\psi_{n}^{\text{IPW}}$ to exhibit
desirable asymptotic properties (consistency, efficiency); in such cases, the
IPW estimator is only unbiased if the parametric conditional density estimator
is _correctly specified_. To limit chances for model misspecification, flexible,
data adaptive strategies may be used in estimating $g_{n,A}$; however, IPW
estimators are incompatible with such estimators "out of the box." Instead,
sieve estimation strategies (undersmoothing) must be used to select an estimator
$g_{n,A}$, from among an appropriate family, that allows for optimal estimation
of $\psi$. This issue arises in part because strategies for optimal selection of
$g_{n,A}$ (e.g., cross-validation) optimize for estimation of the conditional
density function itself, ignoring the fact that it is only a nuisance parameter
in the IPW estimation procedure. When `haldensify()` is used for this purpose,
a family of conditional density estimators $g_{n,A,\lambda}$, indexed by the
$\ell_1$ regularization term $\lambda$, are generated, with cross-validation
usually used to select an optimal conditional density estimator along the
trajectory in $\lambda$. When selected in this way, the estimator
$g_{n,A,\lambda_n}$ will fail to yield an IPW estimator with desirable
asymptotic properties, but undersmoothing may be used to select a more
appropriate estimator. The `haldensify` package implements nonparametric IPW
estimators that incorporate several undersmoothing selectors in the
`ipw_shift()` function. For a formal description of these selectors and
numerical experiments examining the performance of the resultant IPW estimators,
see @hejazi2022efficient.

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
estimation requires estimation of the generalized propensity score as an
intermediate step, this latter function internally calls the former; moreover,
the `ipw_shift()` function and the various selectors to which it provides access
(e.g., `selector_gcv()` for estimator selection based on "global"
cross-validation) have been studied from a methodological perspective in
@hejazi2022efficient. The `haldensify()` function is complemented by appropriate
`predict()` and `plot()` methods, the former to allow for the estimated
conditional density to be evaluated at new values of the variable of interest
and its conditioning set and the latter to visualize the resultant estimators
along the regularization trajectory. The `ipw_shift()` function is accompanied
by a corresponding `confint()` method to easily generate confidence intervals
around the IPW point estimates. The S3 classes returned by both of these
functions have custom `print()` methods to allow for their results to be easily
inspected. Several internal utility functions, including, for example,
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
https://github.com/nhejazi/haldensify. To date (mid-May 2022), CRAN records
indicate that `haldensify` has been downloaded over 13,000 times.

# Acknowledgments

NSH's contributions to this work were supported in part by a grant from the
National Science Foundation (award number [DMS
2102840](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2102840)).

# References

