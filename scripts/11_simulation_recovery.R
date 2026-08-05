# =============================================================================
# 11_simulation_recovery.R - Synthetic parameter-recovery validation (Fig S2)
# =============================================================================
# Simulates noisy normalised-O2 trajectories from the model across grids of true
# r, true per-cell R, sampling interval, duration and noise, mimics trimming
# (stop at 40% of start), re-fits with nlsLM, and quantifies recovery of r and R.
# Fully synthetic - needs no input data.
#
#   Output: results/tables/synthetic_parameter_recovery_results.csv
#           results/tables/Table_S2_synthetic_parameter_recovery_summary.csv
#   Figure: results/figures/Fig_S2_synthetic_param_recovery_combined.(pdf|png)
# =============================================================================

# ---- locate scripts/ dir and source shared config ---------------------------
.this_dir <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  d <- if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE)) else
         tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (is.null(d) || is.na(d) || !nzchar(d)) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
        nzchar(rstudioapi::getActiveDocumentContext()$path))
      d <- dirname(rstudioapi::getActiveDocumentContext()$path) else d <- getwd()
  }
  d
})
source(file.path(.this_dir, "config.R"))
suppressPackageStartupMessages({
  library(minpack.lm)
  library(patchwork)
})

set.seed(123)

# ---- true parameter ranges & simulation design ------------------------------
r_vals <- c(0.01, 0.02, 0.035)       # min^-1 (growth rates)
R_vals <- c(2e-12, 5e-12, 8e-12)     # mg O2 cell^-1 min^-1 (per-cell respiration)

N0_true     <- 1e9   # cells L^-1
O2_ref_true <- 1     # mg O2 L^-1
O2_0_true   <- 1     # initial normalised O2

noise_sd_vec <- c(0.002, 0.005, 0.01)
dt_vec       <- c(1, 3, 5)
dur_vec      <- c(60, 120, 180)
n_reps       <- 50

# ---- simulate one series & refit --------------------------------------------
simulate_and_fit <- function(r_true, R_true, N0_true, O2_0_true, O2_ref_true,
                             noise_sd, dt, duration) {
  K_true <- (R_true * N0_true) / O2_ref_true
  t_full <- seq(0, duration, by = dt)
  O2_true_full <- resp_model(r = r_true, K = K_true, t = t_full, O2_0 = O2_0_true)

  target_O2 <- 0.4 * O2_0_true
  idx_trim  <- which(O2_true_full <= target_O2)[1]
  if (!is.na(idx_trim) && idx_trim > 5) {
    t <- t_full[1:idx_trim]; O2_true <- O2_true_full[1:idx_trim]
  } else {
    t <- t_full; O2_true <- O2_true_full
  }

  O2_obs <- O2_true + rnorm(length(t), mean = 0, sd = noise_sd)
  O2_obs[O2_obs <= 0] <- 1e-6
  df_sim <- tibble::tibble(Time = t, Oxygen_norm = O2_obs)

  r_start <- r_true * runif(1, 0.5, 1.5)
  R_start <- R_true * runif(1, 0.5, 1.5)
  K_start <- (R_start * N0_true) / O2_ref_true

  fit <- tryCatch(
    nlsLM(Oxygen_norm ~ resp_model(r, K, Time, O2_0), data = df_sim,
          start = list(r = r_start, K = K_start, O2_0 = 1),
          lower = c(r = 1e-4, K = 1e-7, O2_0 = 0.8),
          upper = c(r = 0.1,  K = 0.5,  O2_0 = 1.2),
          control = nls.lm.control(maxiter = 300)),
    error = function(e) NULL)

  if (is.null(fit)) {
    return(tibble::tibble(
      r_true = r_true, K_true = K_true, R_true = R_true,
      noise_sd = noise_sd, dt = dt, duration = duration, converged = FALSE,
      r_est = NA_real_, K_est = NA_real_, R_est = NA_real_,
      O2_0_est = NA_real_, pseudo_R2 = NA_real_))
  }

  pred <- predict(fit, df_sim)
  pseudo_R2 <- 1 - sum((df_sim$Oxygen_norm - pred)^2) /
    sum((df_sim$Oxygen_norm - mean(df_sim$Oxygen_norm))^2)
  coefs <- coef(fit)
  r_est <- coefs[["r"]]; K_est <- coefs[["K"]]
  R_est <- K_est * O2_ref_true / N0_true

  tibble::tibble(
    r_true = r_true, K_true = K_true, R_true = R_true,
    noise_sd = noise_sd, dt = dt, duration = duration, converged = TRUE,
    r_est = r_est, K_est = K_est, R_est = R_est,
    O2_0_est = coefs[["O2_0"]], pseudo_R2 = pseudo_R2)
}

# ---- run all scenarios ------------------------------------------------------
sim_grid <- tidyr::expand_grid(
  r_true = r_vals, R_true = R_vals, noise_sd = noise_sd_vec,
  dt = dt_vec, duration = dur_vec, rep_id = seq_len(n_reps))

sim_results <- sim_grid %>%
  purrr::pmap_dfr(function(r_true, R_true, noise_sd, dt, duration, rep_id) {
    simulate_and_fit(r_true, R_true, N0_true, O2_0_true, O2_ref_true,
                     noise_sd, dt, duration)
  })

readr::write_csv(sim_results, tbl("synthetic_parameter_recovery_results.csv"))

# ---- summarise recovery -----------------------------------------------------
sim_summary <- sim_results %>%
  dplyr::filter(converged) %>%
  dplyr::mutate(
    r_error = r_est - r_true, K_error = K_est - K_true, R_error = R_est - R_true,
    r_rel_err = (r_est - r_true) / r_true, K_rel_err = (K_est - K_true) / K_true,
    R_rel_err = (R_est - R_true) / R_true
  ) %>%
  dplyr::group_by(noise_sd, dt, duration) %>%
  dplyr::summarise(
    n_fits = dplyr::n(),
    r_bias = mean(r_error, na.rm = TRUE), r_rmse = sqrt(mean(r_error^2, na.rm = TRUE)),
    K_bias = mean(K_error, na.rm = TRUE), K_rmse = sqrt(mean(K_error^2, na.rm = TRUE)),
    R_bias = mean(R_error, na.rm = TRUE), R_rmse = sqrt(mean(R_error^2, na.rm = TRUE)),
    r_rel_bias = mean(r_rel_err, na.rm = TRUE), K_rel_bias = mean(K_rel_err, na.rm = TRUE),
    R_rel_bias = mean(R_rel_err, na.rm = TRUE), mean_R2 = mean(pseudo_R2, na.rm = TRUE),
    .groups = "drop")

readr::write_csv(sim_summary, tbl("Table_S2_synthetic_parameter_recovery_summary.csv"))

# ---- Fig S2: true vs estimated r and R --------------------------------------
plot_subset <- sim_results %>% dplyr::filter(converged, dt == 1, duration == 180)

plot_r <- ggplot2::ggplot(plot_subset,
                          ggplot2::aes(x = r_true, y = r_est, colour = factor(noise_sd))) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::labs(x = expression("True " * r ~ (min^{-1})),
                y = expression("Estimated " * r ~ (min^{-1})), colour = "Noise SD") +
  ggplot2::theme_classic(base_size = 12)

plot_R <- ggplot2::ggplot(plot_subset,
                          ggplot2::aes(x = R_true, y = R_est, colour = factor(noise_sd))) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::geom_point(alpha = 0.6) +
  ggplot2::labs(x = expression("True " * R ~ (mg~O[2]~cell^{-1}~min^{-1})),
                y = expression("Estimated " * R ~ (mg~O[2]~cell^{-1}~min^{-1})), colour = "Noise SD") +
  ggplot2::theme_classic(base_size = 12)

synthetic_fig <- plot_r + plot_R + patchwork::plot_layout(nrow = 1)

ggplot2::ggsave(fig("Fig_S2_synthetic_param_recovery_combined.pdf"), synthetic_fig,
                width = 9, height = 4, dpi = 300)
ggplot2::ggsave(fig("Fig_S2_synthetic_param_recovery_combined.png"), synthetic_fig,
                width = 9, height = 4, dpi = 300)

message("09_simulation_recovery: ", sum(sim_results$converged), "/", nrow(sim_results),
        " fits converged. Fig S2 written.")
