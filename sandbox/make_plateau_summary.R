# packages and settings
library(here)
library(tidyverse)
library(latex2exp)
library(patchwork)
conflicted::conflict_prefer("filter", "dplyr")
set.seed(72439)

# preset number of example experimental datasets
n_exp <- 4

# make bootstrap MSE plot
p_mse <- lapply(seq_len(n_exp), function(iter) {
  # load example analysis dataset
  file_name <- paste0("dgp1a_haldensify_bootstrap_testing123_",
                      "exp", iter, "_summary", ".rds")
  sim_metrics <- readRDS(here("data", file_name))

  # find the MSE minimizer and its index
  mse_min_idx <- which.min(sim_metrics$boot$mse_boot)
  mse_min <- min(sim_metrics$boot$mse_boot)
  psi_minmse <- sim_metrics$psi[1, mse_min_idx]
  lambda_minmse <- sim_metrics$boot$lambda[mse_min_idx]

  # extract MSE metrics for better plotting
  mse_max <- max(sim_metrics$boot$mse_boot)
  mse_sd <- sd(sim_metrics$boot$mse_boot)
  plot_max <- mse_max + 2 * mse_sd
  plot_min <- mse_min - 2 * mse_sd

  # create a summary plot
  p_mse <- sim_metrics$boot %>%
    ggplot(aes(x = -log(lambda), y = mse_boot)) +
    geom_point(alpha = 0.5) +
    geom_vline(xintercept = -log(lambda_minmse), linetype = "dashed") +
    geom_hline(yintercept = mse_min, linetype = "dashed") +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    geom_smooth(method = "gam", se = FALSE, color = "red") +
    labs(
      x = TeX("$-\\log(\\lambda)$"),
      y = "Mean squared error",
      title = TeX(
        paste0("Exp. ", iter, ": ", "Bootstrap MSE of IPW in $\\lambda$")
      ),
      subtitle = TeX(
        paste("min. MSE = ", round(mse_min, 5), "at",
              paste0(mse_min_idx, "th"), "$\\lambda$,", "with $\\psi_n =$",
              round(psi_minmse, 3)
        )
      )
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 32),
      axis.text.x = element_text(angle = 20, colour = "black",
                                 size = 30, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 30)
    ) +
    ylim(plot_min, plot_max)
})

# save summary plot
p_mse_all <- p_mse[[1]] + p_mse[[2]] + p_mse[[3]] + p_mse[[4]]
ggsave(filename = here("graphs", "plateau_exp", "mse_bootstrap.pdf"),
       plot = p_mse_all, height = 28, width = 33)

###############################################################################

# make bootstrap variance plot
p_var <- lapply(seq_len(n_exp), function(iter) {
  # load example analysis dataset
  file_name <- paste0("dgp1a_haldensify_bootstrap_testing123_",
                      "exp", iter, "_summary", ".rds")
  sim_metrics <- readRDS(here("data", file_name))

  # find the MSE minimizer and its index
  var_min_idx <- which.min(sim_metrics$boot$var_boot)
  var_min <- min(sim_metrics$boot$var_boot)

  # create a summary plot
  p_var <- sim_metrics$boot %>%
    ggplot(aes(x = -log(lambda), y = var_boot)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = var_min, linetype = "dashed") +
    #geom_line(linetype = "dotted") +
    #geom_smooth(method = "loess", se = FALSE, color = "blue") +
    #geom_smooth(method = "gam", se = FALSE, color = "red") +
    labs(
      x = TeX("$-\\log(\\lambda)$"),
      y = "Variance",
      title = TeX(
        paste0("Exp. ", iter, ": ", "Bootstrap Variance of IPW in $\\lambda$")
      ),
      subtitle = TeX(
        paste("min. variance = ", round(var_min, 6), "at",
              paste0(var_min_idx, "th"), "$\\lambda$")
      )
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 32),
      axis.text.x = element_text(angle = 20, colour = "black",
                                 size = 30, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 30)
    )
})

# save summary plot
p_var_all <- p_var[[1]] + p_var[[2]] + p_var[[3]] + p_var[[4]]
ggsave(filename = here("graphs", "plateau_exp", "var_bootstrap.pdf"),
       plot = p_var_all, height = 28, width = 33)

###############################################################################

# make Lepski-based selector plots
p_lepski <- lapply(seq_len(n_exp), function(iter) {
  # load example analysis dataset
  file_name <- paste0("dgp1a_haldensify_bootstrap_testing123_",
                      "exp", iter, "_summary", ".rds")
  sim_metrics <- readRDS(here("data", file_name))

  # get z_{1-alpha/2} and length of lambda sequence
  wald_mult <- abs(qnorm(p = (1 - 0.95) / 2))
  n_lambda <- nrow(sim_metrics$boot)

  # find the variance minimizer and its index
  var_min_idx <- which.min(sim_metrics$boot$var_boot)
  var_min <- min(sim_metrics$boot$var_boot)

  # reduce candidate lambdas by starting at variance minimizer
  boot_metrics_reduced <- sim_metrics$boot[var_min_idx:n_lambda, ]

  # compute SE based on bootstrap MSE and smooth via LOESS
  #se_boot <- sqrt(boot_metrics_reduced$mse_boot)
  #se_boot_smooth <- lowess(x = se_boot)$y
  se_boot <- sqrt(boot_metrics_reduced$var_boot)

  # where do the point estimate and standard error change sign?
  psi_orig <- sim_metrics$psi[1, var_min_idx:n_lambda]
  lepski_crit <- abs(diff(psi_orig)) <= wald_mult * diff(se_boot)
  lepski_idx <- which(lepski_crit)[1] + var_min_idx - 1
  psi_lepski <- sim_metrics$psi[1, lepski_idx]
  lambda_lepski <- sim_metrics$boot$lambda[lepski_idx]

  # create a summary plot
  p_lepski <- sim_metrics$boot %>%
    mutate(
      psi = sim_metrics$psi[1, ]
    ) %>%
    ggplot(aes(x = -log(lambda), y = psi)) +
    geom_point(alpha = 0.5) +
    geom_vline(xintercept = -log(lambda_lepski), linetype = "dashed") +
    geom_hline(yintercept = psi_lepski, linetype = "dashed") +
    geom_line(linetype = "dotted") +
    labs(
      x = TeX("$-\\log(\\lambda)$"),
      y = TeX("$\\psi_n$"),
      title = TeX(
        paste0("Exp. ", iter, ": ", "IPW estimate in $\\lambda$")
      ),
      subtitle = TeX(
        paste("Lepski-type criterion satisified at",
              paste0(lepski_idx, "th"), "$\\lambda$,", "with $\\psi_n =$",
              round(psi_lepski, 3)
        )
      )
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 32),
      axis.text.x = element_text(angle = 20, colour = "black",
                                 size = 30, hjust = 1),
      axis.text.y = element_text(colour = "black", size = 30)
    )
})

# save summary plot
p_lepski_all <- p_lepski[[1]] + p_lepski[[2]] + p_lepski[[3]] + p_lepski[[4]]
ggsave(filename = here("graphs", "plateau_exp", "lepski_bootstrap.pdf"),
       plot = p_lepski_all, height = 28, width = 33)
