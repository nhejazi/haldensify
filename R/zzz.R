.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "haldensify v", utils::packageDescription("haldensify")$Version,
    ": ", utils::packageDescription("haldensify")$Title
  ))
}
