#' Generate Augmented Repeated Measures Data for Pooled Hazards Regression
#'
#' @details Generates an augmented (long format, or repeated measures) dataset
#'  that includes multiple records for each observation, a single record for
#'  each discretized bin up to and including the bin in which a given observed
#'  value of A falls. Such bins are derived from selecting break points over
#'  the support of A. This repeated measures dataset is suitable for estimating
#'  the hazard of failing in a particular bin over A using a highly adaptive
#'  lasso (or other) classification model.
#'
#' @param A The \code{numeric} vector or similar of the observed values of an
#'  intervention for a group of observational units of interest.
#' @param W A \code{data.frame}, \code{matrix}, or similar giving the values of
#'  baseline covariates (potential confounders) for the observed units whose
#'  observed intervention values are provided in the previous argument.
#' @param wts A \code{numeric} vector of observation-level weights. The default
#'  is to weight all observations equally.
#' @param grid_type A \code{character} indicating the strategy (or strategies)
#'  to be used in creating bins along the observed support of the intervention
#'  \code{A}. For bins of equal range, use "equal_range"; consult documentation
#'  of \code{\link[ggplot2]{cut_interval}} for more information. To ensure each
#'  bin has the same number of points, use "equal_mass"; consult documentation
#'  of \code{\link[ggplot2]{cut_number}} for details.
#' @param n_bins Only used if \code{grid_type} is set to \code{"equal_range"}
#'  or \code{"equal_mass"}. This \code{numeric} value indicates the number(s)
#'  of bins into which the support of \code{A} is to be divided.
#' @param breaks A \code{numeric} vector of break points to be used in dividing
#'  up the support of \code{A}. This is passed through the \code{...} argument
#'  to \code{\link[base]{cut.default}} by \code{\link[ggplot2]{cut_interval}}
#'  or \code{\link[ggplot2]{cut_number}}.
#'
#' @importFrom data.table as.data.table setnames
#' @importFrom ggplot2 cut_interval cut_number
#' @importFrom future.apply future_lapply
#' @importFrom assertthat assert_that
#'
#' @return A \code{list} containing the break points used in dividing the
#'  support of \code{A} into discrete bins, the length of each bin, and the
#'  reformatted data. The reformatted data is a \code{\link{data.table}} of
#'  repeated measures data, with an indicator for which bin an observation
#'  fails in, the bin ID, observation ID, values of \code{W} for each given
#'  observation, and observation-level weights.
format_long_hazards <- function(A, W, wts = rep(1, length(A)),
                                grid_type = c(
                                  "equal_range", "equal_mass"
                                ),
                                n_bins = NULL, breaks = NULL) {
  # clean up arguments
  grid_type <- match.arg(grid_type)

  # set grid along A and find interval membership of observations along grid
  if (is.null(breaks) & !is.null(n_bins)) {
    if (grid_type == "equal_range") {
      bins <- ggplot2::cut_interval(A, n_bins,
        right = FALSE,
        ordered_result = TRUE, dig.lab = 12
      )
    } else if (grid_type == "equal_mass") {
      bins <- ggplot2::cut_number(A, n_bins,
        right = FALSE,
        ordered_result = TRUE, dig.lab = 12
      )
    }
    # https://stackoverflow.com/questions/36581075/extract-the-breakpoints-from-cut
    breaks_left <- as.numeric(sub(".(.+),.+", "\\1", levels(bins)))
    breaks_right <- as.numeric(sub(".+,(.+).", "\\1", levels(bins)))
    bin_length <- round(breaks_right - breaks_left, 3)
    bin_id <- as.numeric(bins)
    all_bins <- matrix(seq_len(max(bin_id)), ncol = 1)
    # for predict method, only need to assign observations to existing intervals
  } else if (!is.null(breaks)) {
    # NOTE: findInterval() and cut() might return slightly different results...
    bin_id <- findInterval(A, breaks, all.inside = TRUE)
    all_bins <- matrix(seq_along(breaks), ncol = 1)
  } else {
    stop("Combination of arguments `breaks`, `n_bins` incorrectly specified.")
  }

  # loop over observations to create expanded set of records for each
  reformat_each_obs <- future.apply::future_lapply(seq_along(A), function(i) {
    # create indicator and "turn on" indicator for interval membership
    bin_indicator <- rep(0, nrow(all_bins))
    bin_indicator[bin_id[i]] <- 1
    id <- rep(i, nrow(all_bins))

    # get correct value of baseline variables and repeat along intervals
    if (is.null(dim(W))) {
      # assume vector
      obs_w <- rep(W[i], nrow(all_bins))
      names_w <- "W"
    } else {
      # assume two-dimensional array
      obs_w <- rep(as.numeric(W[i, ]), nrow(all_bins))
      obs_w <- matrix(obs_w, ncol = ncol(W), byrow = TRUE)

      # use names from array if present
      if (is.null(names(W))) {
        names_w <- paste("W", seq_len(ncol(W)), sep = "_")
      } else {
        names_w <- names(W)
      }
    }

    # get correct value of weights and repeat along intervals
    # NOTE: the weights are always a vector
    obs_wts <- rep(wts[i], nrow(all_bins))

    # create data table with membership indicator and interval limits
    suppressWarnings(
      hazards_df <- data.table::as.data.table(cbind(
        id, bin_indicator,
        all_bins, obs_w,
        obs_wts
      ))
    )

    # trim records to simply end at the failure time for a given observation
    hazards_df_reduced <- hazards_df[seq_len(bin_id[i]), ]

    # give explicit names and add to appropriate position in list
    hazards_df <-
      data.table::setnames(
        hazards_df_reduced,
        c("obs_id", "in_bin", "bin_id", names_w, "wts")
      )
    return(hazards_df)
  })

  # combine observation-level hazards data into larger structure
  reformatted_data <- do.call(rbind, reformat_each_obs)
  out <- list(
    data = reformatted_data,
    breaks =
      if (exists("breaks_left")) {
        breaks_left
      } else {
        NULL
      },
    bin_length =
      if (exists("bin_length")) {
        bin_length
      } else {
        NULL
      }
  )
  return(out)
}

###############################################################################

#' Map Predicted Hazard to Predicted Density for a Single Observation
#'
#' @details For a single observation, map a predicted hazard of failure (as an
#'  occurrence in a particular bin, under a given partitioning of the support)
#'  to a density.
#'
#' @param hazard_pred_single_obs A \code{numeric} vector of predicted hazard of
#'  failure in a given bin (under a given partitioning of the support) for a
#'  single observational unit based on a long format data structure (from
#'  \code{\link{format_long_hazards}}). This is the probability that a given
#'  value falls in a corresponding bin, given that it has not yet failed
#'  (fallen in a preceding bin), as per \insertRef{diaz2011super}{haldensify}.
#'
#' @importFrom assertthat assert_that
#'
#' @return A \code{matrix} composed of a single row and a number of columns
#'  specified by the grid of penalization parameters used in fitting of the
#'  highly adaptive lasso. This is the predicted conditional density for a
#'  given observation, re-mapped from the hazard scale.
map_hazard_to_density <- function(hazard_pred_single_obs) {
  # number of records for the given observation
  n_records <- nrow(hazard_pred_single_obs)

  # NOTE: pred_hazard = (1 - pred) if 0 in this bin * pred if 1 in this bin
  if (n_records > 1) {
    hazard_prefailure <- matrix(1 - hazard_pred_single_obs[-n_records, ],
      nrow = (n_records - 1)
    )
    hazard_at_failure <- hazard_pred_single_obs[n_records, ]
    hazard_predicted <- rbind(hazard_prefailure, hazard_at_failure)
    rownames(hazard_predicted) <- NULL
  } else {
    hazard_predicted <- hazard_pred_single_obs
  }

  # sanity check of dimensions
  assertthat::assert_that(all(dim(hazard_pred_single_obs) ==
    dim(hazard_predicted)))

  # multiply hazards across rows to construct the individual-level density
  density_pred_from_hazards <- matrix(apply(hazard_predicted, 2, prod),
    nrow = 1
  )
  return(density_pred_from_hazards)
}

###############################################################################

#' Histogram Binning Procedures for Pooled Hazards Regression
#'
#' @param grid_var The \code{numeric} vector over which histogram-based binning
#'  is to be performed.
#' @param grid_type A \code{character} indicating the choice of binning rule,
#'  with \code{"hist"} corresponding to the use of several rules proposed for
#'  optimal histogram construction and \code{"scaled"} corresponding to the use
#'  of various pre-set multiples of the square root of the sample size.
#' @param max_bins A \code{numeric} indicating the maximum number of bins that
#'  are allowed in the grid for building the histogram based discretization.
#'
#' @importFrom stats IQR
#' @importFrom dplyr case_when
#'
#' @keywords internal
make_bins <- function(grid_var,
                      grid_type = c("hist", "scaled"),
                      max_bins = 30L) {
  # set default grid type
  grid_type <- match.arg(grid_type)

  # make bins based on sample size and binning rules
  n_obs <- length(grid_var)
  if (grid_type == "hist") {
    # histogram binning rules for pooled hazards
    ## 1) simplest rule, just based on the root-n
    k_sqrt <- ceiling(sqrt(n_obs))
    ## 2) Rice's rule looking at n^(1/3)
    k_rice <- ceiling(2 * (n_obs^(1 / 3)))
    ## 3) Freedman-Diaconis rule, also based on n^{1/3}
    h_diaconis <- 2 * (stats::IQR(grid_var) / (n_obs^(1 / 3)))
    k_diaconis <- ceiling((max(grid_var) - min(grid_var)) / h_diaconis)

    # take multiples of each of the numbers of bins
    max_k_bins <- round(c(1.5, 2, 2.5) * max(c(k_sqrt, k_rice, k_diaconis)))

    # construct grid
    bin_grid <- sort(unique(c(k_sqrt, k_rice, k_diaconis, max_k_bins)))
  } else if (grid_type == "scaled") {
    # set different multiplers for root-n based on sample size
    bin_mult <- dplyr::case_when(
      n_obs >= 900 ~ c(0.5, 0.75, 1.0, 1.25),
      TRUE ~ c(0.5, 1.0, 1.5, 2.0)
    )

    # construct grid
    bin_grid <- round(sqrt(n_obs) * bin_mult)
  }

  # return grid of bins
  bin_grid <- bin_grid[bin_grid < max_bins]
  return(bin_grid)
}

###############################################################################

#' Print: Highly Adaptive Lasso Conditional Density Estimates
#'
#' @details The \code{print} method for objects of class \code{haldensify}
#'
#' @param x An object of class \code{haldensify}.
#' @param ... Other options (not currently used).
#'
#' @method print haldensify
#'
#' @importFrom utils head
#'
#' @return None. Called for the side effect of printing an informative summary
#'  of slots of objects of class \code{haldensify}.
#'
#' @export
#'
#' @examples
#' # simulate data: W ~ U[-4, 4] and A|W ~ N(mu = W, sd = 0.5)
#' set.seed(429153)
#' n_train <- 50
#' w <- runif(n_train, -4, 4)
#' a <- rnorm(n_train, w, 0.5)
#'
#' # learn relationship A|W using HAL-based density estimation procedure
#' haldensify_fit <- haldensify(
#'   A = a, W = w, n_bins = c(3, 5),
#'   lambda_seq = exp(seq(-1, -15, length = 50L)),
#'   max_degree = 3, reduce_basis = 0.1
#' )
#' print(haldensify_fit)
print.haldensify <- function(x, ...) {
  # construct and print output
  message("HAL Conditional Density Estimation")
  message(
    "Number of bins over support of A: ",
    x$n_bins_cvselect
  )
  message(
    "CV-selected lambda: ",
    round(x$cv_tuning_results$lambda_loss_min, 4)
  )
  message(
    "Summary of fitted HAL:"
  )
  suppressWarnings(
    print(utils::head(summary(x$hal_fit)$table, 10))
  )
}

###############################################################################

#' Print: IPW Estimates of the Causal Effects of Stochatic Shift Interventions
#'
#' @details The \code{print} method for objects of class \code{ipw_haldensify}
#'
#' @param x An object of class \code{ipw_haldensify}.
#' @param ... Other options (not currently used).
#' @param ci_level A \code{numeric} indicating the level of the confidence
#'  interval to be computed.
#'
#' @method print ipw_haldensify
#'
#' @importFrom stats confint
#' @importFrom scales percent
#' @importFrom dplyr case_when
#' @importFrom matrixStats colVars
#'
#' @return None. Called for the side effect of printing an informative summary
#'  of slots of objects of class \code{ipw_haldensify}.
#'
#' @export
#'
#' @examples
#' # simulate data
#' n_obs <- 50
#' W1 <- rbinom(n_obs, 1, 0.6)
#' W2 <- rbinom(n_obs, 1, 0.2)
#' A <- rnorm(n_obs, (2 * W1 - W2 - W1 * W2), 2)
#' Y <- rbinom(n_obs, 1, plogis(3 * A + W1 + W2 - W1 * W2))
#'
#' # fit the IPW estimator
#' est_ipw_shift <- ipw_shift(
#'   W = cbind(W1, W2), A = A, Y = Y,
#'   delta = 0.5, n_bins = 3L, cv_folds = 2L,
#'   lambda_seq = exp(seq(-1, -10, length = 100L)),
#'   # arguments passed to hal9001::fit_hal()
#'   max_degree = 1,
#'   # ...continue arguments for IPW
#'   selector_type = "gcv"
#' )
#' print(est_ipw_shift)
print.ipw_haldensify <- function(x, ..., ci_level = 0.95) {
  # compute confidence interval
  ci_est <- stats::confint(x, level = ci_level)

  # dictionary of human-readable names for estimator variants
  x_est <- x$est %>%
    dplyr::mutate(
      name = dplyr::case_when(
        x$est$type == "gcv" ~ "Global CV",
        x$est$type == "dcar_tol" ~ "D_CAR Minimizer (Tolerance)",
        x$est$type == "dcar_min" ~ "D_CAR Minimizer (Absolute)",
        x$est$type == "lepski_plateau" ~ "Plateau (Lepski)",
        x$est$type == "smooth_plateau" ~ "Plateau (Smoothed)",
        x$est$type == "hybrid_plateau" ~ "Plateau (Hybrid)"
      )
    )

  # display only the _most efficient_ estimator (minimal variance)
  if (nrow(x$est) > 1L) {
    idx_eff <- which.min(matrixStats::colVars(as.matrix(x$eif)))
    x_eif <- x$eif[, idx_eff]
  } else {
    idx_eff <- 1L
    x_eif <- x$eif
  }
  x_est <- x_est[idx_eff, ]
  ci_est <- ci_est[idx_eff, ]

  # construct and print output
  message("Counterfactual Mean of Shifted Treatment")
  message("Intervention: ", "Treatment + ", x$.delta)
  message("IPW Estimator Criterion: ", x_est$name)
  message("Estimate: ", round(x_est$psi, 4L))
  message("Std. Error: ", round(x_est$se_est, 4L))
  message(paste0(
    scales::percent(ci_level), " CI: [",
    round(ci_est$lwr_ci, 4L), ", ", round(ci_est$upr_ci, 4L), "]"
  ))
  message("EIF Mean: ", round(mean(x_eif), 4L))
}

###############################################################################

is.haldensify <- function(x) {
  class(x) == "haldensify"
}

is.ipw_haldensify <- function(x) {
  class(x) == "ipw_haldensify"
}
