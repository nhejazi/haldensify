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
#' @param new_A The \code{numeric} vector or similar of the observed values for
#'  which a conditional density estimate is to be generated.
#' @param new_W A \code{data.frame}, \code{matrix}, or similar giving the
#'  values of baseline covariates (potential confounders) for the conditioning
#'  set of the observed values \code{A}.
#' @param cv_select A \code{logical} indicating whether to return the predicted
#'  density for the value of the regularization parameter selected by global
#'  cross-validation. The default is \code{TRUE}. When set to \code{FALSE}, a
#'  matrix of predicted densities is returned, with each column corresponding
#'  to a value of the regularization parameter less than or equal to the choice
#'  made by the global cross-validation selector.
#'
#' @importFrom stats predict
#' @importFrom data.table ":="
#'
#' @return A \code{numeric} vector of predicted conditional density values from
#'  a fitted \code{haldensify} object.
#'
#' @export
#'
#' @examples
#' # simulate data: W ~ U[-4, 4] and A|W ~ N(mu = W, sd = 0.5)
#' n_train <- 50
#' w <- runif(n_train, -4, 4)
#' a <- rnorm(n_train, w, 0.5)
#' # HAL-based density estimator of A|W
#' mod_haldensify <- haldensify(
#'   A = a, W = w, n_bins = 3,
#'   lambda_seq = exp(seq(-1, -10, length = 50))
#' )
#' # predictions to recover conditional density of A|W
#' new_a <- seq(-4, 4, by = 0.1)
#' new_w <- rep(0, length(new_a))
#' pred_dens <- predict(mod_haldensify, new_A = new_a, new_W = new_w)
predict.haldensify <- function(object, ..., new_A, new_W,
                               cv_select = TRUE) {
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

  # predict conditional density estimate from HAL fit on new long format data
  # over the sequence of lambda less than or equal to CV-selected lambda
  hazard_pred <- stats::predict(
    object = object$hal_fit,
    new_data =
      long_data_pred[, 3:ncol(long_data_pred)]
  )

  # NOTE: we return hazard predictions for the loss minimizer and all lambda
  #       smaller than it, BUT if there are no such lambda, hazard_pred is only
  #       vector rather than the usually expected matrix
  if (!is.matrix(hazard_pred)) hazard_pred <- as.matrix(hazard_pred, ncol = 1)

  # estimate unscaled density for each observation and each lambda
  density_pred_rescaled <- apply(hazard_pred, 2, function(this_hazard_pred) {
    # coerce single column of predictions back to matrix
    this_hazard_pred <- matrix(this_hazard_pred, ncol = 1)

    # compute hazard for a given observation by looping over individuals
    dens_given_lambda <- lapply(unique(long_data_pred$obs_id), function(id) {
      # get predictions for the current observation only
      hazard_pred_this_obs <-
        matrix(this_hazard_pred[long_data_pred$obs_id == id, ])

      # map hazard to density for a single observation and return
      density_pred_this_obs <-
        map_hazard_to_density(hazard_pred_single_obs = hazard_pred_this_obs)

      # return density for a single observation
      return(as.numeric(density_pred_this_obs))
    })
    # aggregate predicted unscaled density at the level of observations
    density_pred_unscaled <- do.call(c, dens_given_lambda)

    # re-scale to densities by dividing by bin widths
    density_pred_scaled <- density_pred_unscaled /
      object$bin_sizes[long_data_pred[in_bin == 1, bin_id]]

    # return re-scaled densities
    return(density_pred_scaled)
  })

  # truncate predictions outside range of observed A
  outside_range <- new_A < object$range_a[1] | new_A > object$range_a[2]
  density_pred_rescaled[outside_range, ] <- 0

  # return predicted densities only for CV-selected lambda
  # NOTE: turn off for access to density estimates for all lambda >= CV-lambda
  if (cv_select) {
    density_pred_rescaled <- density_pred_rescaled[, 1]
  }

  # output
  return(density_pred_rescaled)
}
