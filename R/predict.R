utils::globalVariables(c("wts"))

#' Prediction method for HAL-based conditional density estimation
#'
#' @param object An object of class \code{\link{haldensify}}, containing the
#'  results of fitting the highly adaptive lasso for conditional density
#'  estimation, as produced by a call to \code{\link{haldensify}}.
#' @param ... Additional arguments passed to \code{predict} as necessary.
#' @param new_A The \code{numeric} vector or similar of the observed values of
#'  an intervention for a group of observational units of interest.
#' @param new_W A \code{data.frame}, \code{matrix}, or similar giving the values
#'  of baseline covariates (potential confounders) for the observed units whose
#'  observed intervention values are provided in the previous argument.
#'
#' @importFrom stats predict
#' @importFrom future.apply future_lapply
#'
#' @export
#
predict.haldensify <- function(object, ..., new_A, new_W) {
  # make long format data structure with new input data
  long_format_args <- list(
    A = new_A,
    W = new_W,
    type = object$call$grid_type,
    breaks = object$breaks
  )
  reformatted_output <- do.call(format_long_hazards, long_format_args)
  long_data_pred <- reformatted_output$data
  long_data_pred[, wts := NULL]

  # predict conditional density estimate from HAL fit on new long format data
  hazard_pred <-
    stats::predict(object$hal_fit,
      new_data =
        long_data_pred[, 3:ncol(long_data_pred)]
    )

  # compute hazard for a given observation by looping over individuals
  density_pred_each_obs <-
    future.apply::future_lapply(unique(long_data_pred$obs_id), function(id) {
      # get predictions for the current observation only
      hazard_pred_this_obs <- matrix(hazard_pred[long_data_pred$obs_id == id])

      # map hazard to density for a single observation and return
      density_pred_this_obs <-
        map_hazard_to_density(hazard_pred_single_obs = hazard_pred_this_obs)

      # return density for a single observation
      return(as.numeric(density_pred_this_obs))
    })

  # aggregate predicted density at the level of observations
  density_pred_unscaled <- do.call(c, density_pred_each_obs)
  density_pred_scaled <- density_pred_unscaled /
    object$bin_sizes[long_data_pred[in_bin == 1, bin_id]]

  # truncate predictions outside range of observed A
  density_pred_scaled[new_A < object$range_a[1] | new_A > object$range_a[2]] <- 0
  
  # output
  return(density_pred_scaled)
}
