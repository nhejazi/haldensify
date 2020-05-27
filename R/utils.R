#' Generate long format hazards data for pooled hazards estimation
#'
#' @details Generates a long-form dataset that represents each observation in
#'  terms of repeated measures across discretized bins derived from selecting
#'  break points over the support of A. This repeated measures dataset is
#'  suitable for estimating the hazard of failing in a particular bin over A
#'  using a highly adaptive lasso classification model.
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

#' Map a predicted hazard to a predicted density for a single observation
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
