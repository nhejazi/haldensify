% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{format_long_hazards}
\alias{format_long_hazards}
\title{Generate Augmented Repeated Measures Data for Pooled Hazards Regression}
\usage{
format_long_hazards(
  A,
  W,
  wts = rep(1, length(A)),
  grid_type = c("equal_range", "equal_mass"),
  n_bins = NULL,
  breaks = NULL
)
}
\arguments{
\item{A}{The \code{numeric} vector or similar of the observed values of an
intervention for a group of observational units of interest.}

\item{W}{A \code{data.frame}, \code{matrix}, or similar giving the values of
baseline covariates (potential confounders) for the observed units whose
observed intervention values are provided in the previous argument.}

\item{wts}{A \code{numeric} vector of observation-level weights. The default
is to weight all observations equally.}

\item{grid_type}{A \code{character} indicating the strategy (or strategies)
to be used in creating bins along the observed support of the intervention
\code{A}. For bins of equal range, use "equal_range"; consult documentation
of \code{\link[ggplot2]{cut_interval}} for more information. To ensure each
bin has the same number of points, use "equal_mass"; consult documentation
of \code{\link[ggplot2]{cut_number}} for details.}

\item{n_bins}{Only used if \code{grid_type} is set to \code{"equal_range"}
or \code{"equal_mass"}. This \code{numeric} value indicates the number(s)
of bins into which the support of \code{A} is to be divided.}

\item{breaks}{A \code{numeric} vector of break points to be used in dividing
up the support of \code{A}. This is passed through the \code{...} argument
to \code{\link[base]{cut.default}} by \code{\link[ggplot2]{cut_interval}}
or \code{\link[ggplot2]{cut_number}}.}
}
\value{
A \code{list} containing the break points used in dividing the
 support of \code{A} into discrete bins, the length of each bin, and the
 reformatted data. The reformatted data is a \code{\link{data.table}} of
 repeated measures data, with an indicator for which bin an observation
 fails in, the bin ID, observation ID, values of \code{W} for each given
 observation, and observation-level weights.
}
\description{
Generate Augmented Repeated Measures Data for Pooled Hazards Regression
}
\details{
Generates an augmented (long format, or repeated measures) dataset
 that includes multiple records for each observation, a single record for
 each discretized bin up to and including the bin in which a given observed
 value of A falls. Such bins are derived from selecting break points over
 the support of A. This repeated measures dataset is suitable for estimating
 the hazard of failing in a particular bin over A using a highly adaptive
 lasso (or other) classification model.
}
