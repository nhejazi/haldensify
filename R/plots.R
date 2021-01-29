#' Plot Method for HAL Conditional Density Estimates
#'
#' @param x [TO FILL IN]
#' @param ... [TO FILL IN]
#' @param type [TO FILL IN]
#'
#' @importFrom tibble as_tibble
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line geom_vline xlab
#'  ylab ggtitle theme_bw
#'
#' @return [TO FILL IN]
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
#'   lambda_seq = exp(seq(-1, -10, length = 50))
#' )
#' plot(haldensify_fit, type = "risk")
plot.haldensify <- function(x, ..., type = c("risk", "density")) {
  # set default plot type
  type <- match.arg(type)

  #
  emp_risk_data <- tibble::as_tibble(list(
    lambda = x$cv_tuning_results$lambda_seq,
    risk = x$cv_tuning_results$emp_risks
  ))

  # create plot of empirical risks across grid in lambda
  p_risk <- ggplot2::ggplot(
      emp_risk_dat,
      ggplot2::aes_string(x = lambda, y = risk)
    ) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = x$cv_tuning_results$lambda_loss_min,
                        linetype = "dotted") +
    ggplot2::xlab("Lambda (L1 regularization parameter)") +
    ggplot2::ylab("Empirical Risk") +
    ggplot2::ggtitle("Empirical risk of HAL conditional density estimators") +
    ggplot2::theme_bw()
  return(p_risk)
}
