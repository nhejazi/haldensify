utils::globalVariables(c("lambda", "risk"))

#' Plot Method for HAL Conditional Density Estimates
#'
#' @param x Object of class \code{haldensify}, containing conditional density
#'  estimates, as produced by \code{\link{haldensify}}.
#' @param ... Additional arguments to be passed \code{plot}, currently ignored.
#' @param type A \code{character} indicating the type of plot to be produced.
#'  Options include visualizing the empirical risks of the conditional density
#'  estimators across a grid of values of the regularization parameter and a
#'  plot of the estimated conditional density (based on the estimator selected
#'  by cross-validation). The latter has yet to be implemented.
#'
#' @importFrom assertthat assert_that
#' @importFrom data.table as.data.table data.table setnames
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line geom_vline xlab
#'  ylab ggtitle theme_bw
#'
#' @return Object of class \code{ggplot} containing a plot of the desired
#'  \code{type}.
#'
#' @export
#'
#' @examples
#' # simulate data: W ~ U[-4, 4] and A|W ~ N(mu = W, sd = 0.5)
#' n_train <- 50
#' w <- runif(n_train, -4, 4)
#' a <- rnorm(n_train, w, 0.5)
#' # learn relationship A|W using HAL-based density estimation procedure
#' haldensify_fit <- haldensify(
#'   A = a, W = w, n_bins = 3,
#'   lambda_seq = exp(seq(-1, -10, length = 50)),
#'   # the following arguments are passed to hal9001::fit_hal()
#'   max_degree = 3, reduce_basis = 0.1
#' )
#' plot(haldensify_fit)
plot.haldensify <- function(x, ..., type = c("risk", "density")) {
  # set default plot type
  type <- match.arg(type)

  # density plot not yet implemented
  assertthat::assert_that(
    type != "density",
    msg = "Density plot method not yet implemented. Check back later..."
  )

  if (type == "risk") {
    # re-organize object output for plotting empirical risk across lambda
    emp_risk_data <- data.table::data.table(
      lambda = x$cv_tuning_results$lambda_seq,
      risk = x$cv_tuning_results$emp_risks
    )
    data.table::setnames(emp_risk_data, c("lambda", "risk"))
    lambda_cvrisk_min <- x$cv_tuning_results$lambda_loss_min

    # create plot of empirical risks across grid in lambda
    p_emprisk <-
      ggplot2::ggplot(
        emp_risk_data,
        ggplot2::aes(x = -log(lambda), y = risk)
      ) +
      ggplot2::geom_point() +
      ggplot2::geom_line(linetype = "dashed") +
      ggplot2::geom_vline(
        xintercept = -log(lambda_cvrisk_min),
        linetype = "dotted"
      ) +
      ggplot2::labs(
        x = "-log(L1 norm regularization)",
        y = "Empirical risk",
        title = "Empirical risk of HAL conditional density estimators",
        subtitle =
          "(dotted line: L1 norm regularization minimizing CV empirical risk)"
      ) +
      ggplot2::theme_bw()
    return(p_emprisk)
  }
}
