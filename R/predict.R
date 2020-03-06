utils::globalVariables(c("wts"))

#' Prediction method for HAL-based conditional density estimation
#'
#' @details Method for computing and extracting predictions of the conditional
#'  density estimates based on the highly adaptive lasso estimator, returned as
#'  an S3 object of class \code{haldensify} from \code{\link{haldensify}}.
#'
#' @param object An object of class \code{\link{haldensify}}, containing the
#'  results of fitting the highly adaptive lasso for conditional density
#'  estimation, as produced by a call to \code{\link{haldensify}}.
#' @param ... Additional arguments passed to \code{predict} as necessary.
#' @param new_A The \code{numeric} vector or similar of the observed values of
#'  an intervention for a group of observational units of interest.
#' @param new_W A \code{data.frame}, \code{matrix}, or similar giving the
#'  values of baseline covariates (potential confounders) for the observed
#'  units whose observed intervention values are provided in the previous
#'  argument.
#'
#' @importFrom stats predict
#' @importFrom future.apply future_lapply
#'
#' @return A \code{numeric} vector of predicted conditional density values from
#'  a fitted \code{haldensify} object.
#'
#' @examples
#' # simulate data: W ~ U[-4, 4] and A|W ~ N(mu = W, sd = 0.5)
#' n_train <- 50
#' w <- runif(n_train, -4, 4)
#' a <- rnorm(n_train, w, 0.5)
#' # learn relationship A|W using HAL-based density estimation procedure
#' mod_haldensify <- haldensify(
#'   A = a, W = w, n_bins = 3,
#'   lambda_seq = exp(seq(-1, -10, length = 50))
#' )
#' # predictions to recover conditional density of A|W
#' new_a <- seq(-4, 4, by = 0.1)
#' new_w <- rep(0, length(new_a))
#' pred_dens <- predict(mod_haldensify, new_A = new_a, new_W = new_w)
#'
#' @export
predict.haldensify <- function(object, ..., new_A, new_W) {
  # make long format data structure with new input data
  long_format_args <- list(
    A = new_A,
    W = new_W,
    grid_type = object$call$grid_type,
    breaks = object$breaks
  )
  reformatted_output <- do.call(format_long_hazards, long_format_args)
  long_data_pred <- reformatted_output$data
  long_data_pred[, wts := NULL]

  # NOTE: as of v0.2.5 of the hal9001 package, predict.hal9001 checks to see
  #       whether the coefficients are a matrix when fitting a sequence of HAL
  #       models parameterized by lambdas. In our case, we choose a set of
  #       coefficients (from a single lambda) by CV but do this outside of
  #       hal9001::fit_hal, necessitating the following
  object$hal_fit$coefs <- matrix(object$hal_fit$coefs, ncol = 1)

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
  density_pred_scaled[new_A < object$range_a[1] |
    new_A > object$range_a[2]] <- 0

  # output
  return(density_pred_scaled)
}
