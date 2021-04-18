utils::globalVariables(c("wts"))

#' Prediction Method for HAL Conditional Density Estimation
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
#' @param trim_dens A \code{numeric} indicating the minimum allowed value of
#'  the resultant density predictions. Any predicted density values below this
#'  tolerance threshold are set to the indicated minimum. The default is to use
#'  the inverse of the square root of the sample size of the prediction set,
#'  i.e., 1/sqrt(n); another notable choice is 1/sqrt(n)/log(n). If there are
#'  observations in the prediction set with values of \code{new_A} outside of
#'  the support of the training set (i.e., provided in the argument \code{A} to
#'  \code{\link{haldensify}}), their predictions are similarly truncated.
#' @param lambda_select A \code{character} indicating whether to return the
#'  predicted density for the value of the regularization parameter chosen by
#'  the global cross-validation selector or whether to return an undersmoothed
#'  sequence (which starts with the cross-validation selector's choice but also
#'  includes all values in the sequence that are less restrictive). The default
#'  is \code{"cv"} for the global cross-validation selector. Setting the choice
#'  to \code{"undersmooth"} returns a matrix of predicted densities, with each
#'  column corresponding to a value of the regularization parameter less than
#'  or equal to the choice made by the global cross-validation selector. When
#'  \code{"all"} is set, predictions are returned for the full sequence of the
#'  regularization parameter on which the HAL model \code{object} was fitted.
#'
#' @importFrom assertthat assert_that
#' @importFrom data.table ":="
#' @importFrom stats predict
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
#' haldensify_fit <- haldensify(
#'   A = a, W = w, n_bins = c(3, 5),
#'   lambda_seq = exp(seq(-1, -10, length = 500)),
#'   # the following arguments are passed to hal9001::fit_hal()
#'   max_degree = 3, smoothness_orders = 0, reduce_basis = 0.1
#' )
#' # predictions to recover conditional density of A|W
#' new_a <- seq(-4, 4, by = 0.1)
#' new_w <- rep(0, length(new_a))
#' pred_dens <- predict(haldensify_fit, new_A = new_a, new_W = new_w)
predict.haldensify <- function(object,
                               ...,
                               new_A,
                               new_W,
                               trim_dens = 1 / sqrt(length(new_A)),
                               lambda_select = c("cv", "undersmooth", "all")) {
  # set default selection procedure to the cross-validation selector
  lambda_select <- match.arg(lambda_select)
  if (lambda_select %in% c("cv", "undersmooth")) {
    # check existence of CV-selected lambda and extract from model object slot
    assertthat::assert_that(
      !is.na(object$cv_tuning_results$lambda_loss_min_idx),
      msg = "No CV-selected lambda found in fitted haldensify model"
    )
    cv_lambda_idx <- object$cv_tuning_results$lambda_loss_min_idx
  }

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
    new_data = long_data_pred[, -c("obs_id", "in_bin")]
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

  # truncate conditional density estimates below the specified tolerance
  if (min(density_pred_rescaled) < trim_dens) {
    density_pred_rescaled <- apply(density_pred_rescaled, 2, pmax, trim_dens)
  }

  # trim values outside training support to avoid extrapolation
  outside_support <- new_A < object$range_a[1] | new_A > object$range_a[2]
  if (any(outside_support)) {
    density_pred_rescaled[outside_support, ] <- trim_dens
    prop_trimmed <- sum(outside_support) / length(outside_support)
    message(paste0("Trimmed predictions for ", round(100 * prop_trimmed, 2),
                   "% of observations (outside training support)."))
  }

  # return predicted densities only for CV-selected or undersmoothed lambdas
  if (lambda_select == "cv") {
    density_pred_rescaled <- density_pred_rescaled[, cv_lambda_idx]
  } else if (lambda_select == "undersmooth") {
    usm_lambda_idx <- cv_lambda_idx:length(object$cv_tuning_results$lambda_seq)
    density_pred_rescaled <- density_pred_rescaled[, usm_lambda_idx]
  } else if (lambda_select == "all") {
    # pass -- just return conditional density estimates across all lambda
    TRUE
  }
  return(density_pred_rescaled)
}
