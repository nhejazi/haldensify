.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "haldensify v", utils::packageDescription("haldensify")$Version,
    ": Conditional Density Estimation with the Highly Adaptive Lasso"
  ))
}
