set.seed(11249)

# simulate data: W ~ Rademacher and A|W ~ N(mu = \pm 1, sd = 0.5)
n_train <- 100
w <- rbinom(n_train, 1, 0.5)
w[w == 0] <- -1
a <- rnorm(n_train, 2 * w, 0.5)

# specify cut-points in A manually
break_min <- floor(min(a))
break_max <- ceiling(max(a))
breaks <- seq(break_min, break_max, by = 1)

# create pooled hazard data frame with manual cut-points
hazard_df <- format_long_hazards(
  A = a, W = w, breaks = breaks
)

# try again with uneven binning
breaks_neg <- seq(break_min, 0, by = 0.5)
breaks_pos <- seq(0, break_max, by = 0.25)
breaks <- unique(c(breaks_neg, breaks_pos))

# create pooled hazard data frame with manual cut-points
hazard_df <- format_long_hazards(
  A = a, W = w, breaks = breaks
)
