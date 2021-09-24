make_sim_data <- function(n_samp = 1000, dgp_type) {
  # Simulation helper for generating data under a variety of data-generating
  # processes. In each case, the treatment mechanism is normally distributed
  # with the mean (and, sometimes, variance) being functions of the baseline
  # covariates, while the outcome mechanism is a nontrivial function of the
  # exposure and baseline covariates.

  # Scenario definitions
  # - number indicates _same_ outcome mechanism, similar treatment mechanisms
  # - usually, treatment mechanism complexity differs within scenarios

  # SCENARIO 1a: exposure density is N(mu = f(W), sigma2 = c)
  if (dgp_type == "1a") {
    # define treatment mechanism
    g0 <- function(W1, W2, W3) {
      mu <- 2 * W1 + W2 - 2 * (1 - W1) * W2
      sigma2 <- 2
      return(list(mu = mu, sigma2 = sigma2))
    }

    # define outcome mechanism
    Q0 <- function(A, W1, W2, W3) {
      plogis(3 * (A - 1) + W1 + 1.5 * W2^2 - W3 - 3 * (1 - W1) * W2)
    }

    # simulate data
    W1 <- rbinom(n_samp, 1, 0.6)
    W2 <- runif(n_samp, 0.5, 1.5)
    W3 <- rpois(n_samp, 2)
    g_obs <- g0(W1, W2, W3)
    A <- rnorm(n_samp, g_obs$mu, sqrt(g_obs$sigma2))
    Y <- rbinom(n_samp, 1, Q0(A, W1, W2, W3))
    data_obs <- as_tibble(
      list(W1 = W1, W2 = W2, W3 = W3, A = A, Y = Y)
    )
  }

  # SCENARIO 1b: exposure density is N(mu = f_1(W), sigma2 = f_2(W))
  if (dgp_type == "1b") {
    # define treatment mechanism
    g0 <- function(W1, W2, W3) {
      mu <- W1 + 2 * W2 - 2 * (1 - W1) * W2
      sigma2 <- (W1 * 2) + ((1 - W1) * 0.5) + W2
      return(list(mu = mu, sigma2 = sigma2))
    }

    # define outcome mechanism
    Q0 <- function(A, W1, W2, W3) {
      plogis(5 * (A - 2) + W1 + 3 * W2^4 - 2 * W3 - 4 * (1 - W1) * W2)
    }

    # simulate data
    W1 <- rbinom(n_samp, 1, 0.6)
    W2 <- runif(n_samp, 0.5, 1.5)
    W3 <- rpois(n_samp, 2)
    g_obs <- g0(W1, W2, W3)
    A <- rnorm(n_samp, g_obs$mu, sqrt(g_obs$sigma2))
    Y <- rbinom(n_samp, 1, Q0(A, W1, W2, W3))
    data_obs <- as_tibble(
      list(W1 = W1, W2 = W2, W3 = W3, A = A, Y = Y)
    )
  }

  # SCENARIO 2a: exposure density is Poisson(lambda = f(W))
  if (dgp_type == "2a") {
    # define treatment mechanism
    g0 <- function(W1, W2, W3) {
      lambda <- (1 - W1) + 0.25 * W2^3 + 2 * W1 * W2 + 4
      return(list(rate = lambda))
    }

    # define outcome mechanism
    Q0 <- function(A, W1, W2, W3) {
      plogis(A + 2 * (1 - W1) + 0.5 * W2 + 0.5 * W3 + 2 * W1 * W2 - 7)
    }

    # simulate data
    W1 <- rbinom(n_samp, 1, 0.6)
    W2 <- runif(n_samp, 0.5, 1.5)
    W3 <- rpois(n_samp, 2)
    g_obs <- g0(W1, W2, W3)
    A <- rpois(n_samp, lambda = g_obs$rate)
    Y <- rbinom(n_samp, 1, Q0(A, W1, W2, W3))
    data_obs <- as_tibble(
      list(W1 = W1, W2 = W2, W3 = W3, A = A, Y = Y)
    )
  }

  # output
  dgp_funs <- list(g0 = g0, Q0 = Q0)
  return(list(data_obs = data_obs, dgp_funs = dgp_funs))
}

# get truth for a given data-generating mechanism
# NOTE: g0 and Q0 are functions of {A, W1, W2, W3} only
get_truth <- function(n_samp = 1e7,
                      delta,
                      dgp_type = c("1a", "1b", "2a"),
                      param_type = c("pie", "tsm")) {
  # input check
  dgp_type <- match.arg(dgp_type)
  param_type <- match.arg(param_type)

  # compute large data set from data-generating mechanism
  dgp <- make_sim_data(n_samp = n_samp, dgp_type = dgp_type)
  wow_so_much_data <- dgp$data_obs

  # extract helper functions
  g0 <- dgp$dgp_funs$g0
  Q0 <- dgp$dgp_funs$Q0

  # compute treatment mechanism in large sample
  gn_mech <- with(wow_so_much_data, g0(W1, W2, W3))

  # compute likelihood factors
  if (dgp_type %in% c("1a", "1b")) {
    gn_Aminusdelta <- dnorm(wow_so_much_data$A - delta,
      mean = gn_mech$mu,
      sd = sqrt(gn_mech$sigma2)
    )
    gn_Anat <- dnorm(wow_so_much_data$A,
      mean = gn_mech$mu,
      sd = sqrt(gn_mech$sigma2)
    )
  } else if (dgp_type == "2a") {
    gn_Aminusdelta <- dpois(wow_so_much_data$A - delta,
      lambda = gn_mech$rate
    )
    gn_Anat <- dpois(wow_so_much_data$A,
      lambda = gn_mech$rate
    )
  }
  Qbar_noshift <- with(wow_so_much_data, Q0(A, W1, W2, W3))
  Qbar_shift <- with(wow_so_much_data, Q0(A + delta, W1, W2, W3))

  # EIF for shift parameter
  Hn <- gn_Aminusdelta / gn_Anat
  eif_shift <- Hn * (wow_so_much_data$Y - Qbar_noshift) +
    (Qbar_shift - mean(Qbar_shift))

  if (param_type == "pie") {
    # truth in terms of PIE = EY - EY_delta
    EY <- mean(wow_so_much_data$Y)
    eif_EY <- wow_so_much_data$Y - EY
    eif_pie <- eif_EY - eif_shift
    true_psi <- EY - mean(Qbar_shift)
    eff_bound <- var(eif_pie)
  } else if (param_type == "tsm") {
    # truth in terms of TSM EY_delta
    true_psi <- mean(Qbar_shift)
    eff_bound <- var(eif_shift)
  }

  # output
  return(list(
    psi = true_psi,
    eff_bound = eff_bound
  ))
}
