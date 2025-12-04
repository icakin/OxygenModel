############################################################
# README – Synthetic-data validation for the O₂ growth–respiration model
#
# This script:
#   1) Uses the same normalised O₂ model as the main analysis:
#        O2_norm(t) = O2_0 + (resp_tot / r) * (1 - exp(r * t))
#      where resp_tot = (R * N0) / O2_ref.
#   2) Generates synthetic dissolved-oxygen time series across grids of:
#        - growth rate r
#        - per-cell respiration R
#        - noise level (Gaussian SD)
#        - sampling interval (dt)
#        - total duration (minutes)
#   3) Adds noise, mimics trimming to the main depletion phase, and refits the
#      same nlsLM model used for empirical data to recover r and R.
#   4) Quantifies parameter recovery (bias, RMSE, relative error, pseudo-R²)
#      as a function of noise, sampling frequency, and time-series length.
#   5) Produces parameter-recovery plots (true vs estimated r and R) for a
#      representative “realistic” scenario.
#
# Outputs:
#   - Tables (in `Tables/`):
#       * synthetic_parameter_recovery_results.csv
#       * synthetic_parameter_recovery_summary.csv
#   - Plots (in `plots/`):
#       * Fig_S2_synthetic_param_recovery_combined.pdf
#       * Fig_S2_synthetic_param_recovery_combined.png
############################################################

# ---- 0. Libraries & project paths -------------------------------------
library(tidyverse)
library(minpack.lm)
library(patchwork)

data_dir   <- "data"
plots_dir  <- "plots"
tables_dir <- "Tables"

dir.create(data_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Model function (same as in manuscript, normalised form) -------
# O2_norm(t) = O2_0 + (resp_tot / r) * (1 - exp(r * t))
# where resp_tot = (R * N0) / O2_ref  (dimensionless)
resp_model_norm <- function(r, resp_tot, t, O2_0) {
  O2_0 + (resp_tot / r) * (1 - exp(r * t))
}

# ---- 2. True parameter ranges & simulation design ----------------------
set.seed(123)  # reproducible

# Choose realistic ranges (tweak if needed)
r_vals <- c(0.01, 0.02, 0.035)          # min^-1  (growth rates)
R_vals <- c(2e-12, 5e-12, 8e-12)        # mg O2 cell^-1 min^-1  (per-cell respiration)

N0_true     <- 1e9   # cells L^-1 (scale only; used to define resp_tot_true)
O2_ref_true <- 1     # mg O2 L^-1, reference O2 for normalisation
O2_0_true   <- 1     # initial normalised O2

noise_sd_vec <- c(0.002, 0.005, 0.01)   # SD of Gaussian noise on O2_norm
dt_vec       <- c(1, 3, 5)              # sampling interval (min)
dur_vec      <- c(60, 120, 180)         # total duration (min)

n_reps <- 50   # replicate simulations per scenario

# ---- 3. Function: simulate one DO series & refit model -----------------
simulate_and_fit <- function(r_true, R_true, N0_true,
                             O2_0_true, O2_ref_true,
                             noise_sd, dt, duration) {
  
  # True dimensionless total respiration scaling:
  # resp_tot_true = (R * N0) / O2_ref
  resp_tot_true <- (R_true * N0_true) / O2_ref_true
  
  # full time grid
  t_full <- seq(0, duration, by = dt)
  
  # noise-free normalised O2
  O2_true_full <- resp_model_norm(
    r        = r_true,
    resp_tot = resp_tot_true,
    t        = t_full,
    O2_0     = O2_0_true
  )
  
  # mimic trimming: stop when O2 has dropped to 40% of start, if reached
  target_O2 <- 0.4 * O2_0_true
  idx_trim  <- which(O2_true_full <= target_O2)[1]
  if (!is.na(idx_trim) && idx_trim > 5) {
    t       <- t_full[1:idx_trim]
    O2_true <- O2_true_full[1:idx_trim]
  } else {
    t       <- t_full
    O2_true <- O2_true_full
  }
  
  # add Gaussian noise
  O2_obs <- O2_true + rnorm(length(t), mean = 0, sd = noise_sd)
  O2_obs[O2_obs <= 0] <- 1e-6  # avoid non-positive values
  
  df_sim <- tibble(
    Time        = t,
    Oxygen_norm = O2_obs
  )
  
  # starting values (not equal to truth)
  r_start        <- r_true * runif(1, 0.5, 1.5)
  R_start        <- R_true * runif(1, 0.5, 1.5)
  resp_tot_start <- (R_start * N0_true) / O2_ref_true
  
  fit <- tryCatch(
    nlsLM(
      Oxygen_norm ~ resp_model_norm(r, resp_tot, Time, O2_0),
      data    = df_sim,
      start   = list(r = r_start, resp_tot = resp_tot_start, O2_0 = 1),
      lower   = c(r = 1e-4,  resp_tot = 1e-5, O2_0 = 0.8),
      upper   = c(r = 0.1,   resp_tot = 0.1,  O2_0 = 1.2),
      control = nls.lm.control(maxiter = 300)
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(tibble(
      r_true    = r_true,
      R_true    = R_true,
      noise_sd  = noise_sd,
      dt        = dt,
      duration  = duration,
      converged = FALSE,
      r_est     = NA_real_,
      R_est     = NA_real_,
      O2_0_est  = NA_real_,
      pseudo_R2 = NA_real_
    ))
  }
  
  # pseudo-R²
  pred <- predict(fit, df_sim)
  pseudo_R2 <- 1 - sum((df_sim$Oxygen_norm - pred)^2) /
    sum((df_sim$Oxygen_norm - mean(df_sim$Oxygen_norm))^2)
  
  coefs         <- coef(fit)
  r_est         <- coefs[["r"]]
  resp_tot_est  <- coefs[["resp_tot"]]
  
  # Map back to per-cell respiration R:
  # resp_tot = (R * N0) / O2_ref  =>  R_est = resp_tot_est * O2_ref_true / N0_true
  R_est <- resp_tot_est * O2_ref_true / N0_true
  
  tibble(
    r_true    = r_true,
    R_true    = R_true,
    noise_sd  = noise_sd,
    dt        = dt,
    duration  = duration,
    converged = TRUE,
    r_est     = r_est,
    R_est     = R_est,
    O2_0_est  = coefs[["O2_0"]],
    pseudo_R2 = pseudo_R2
  )
}

# ---- 4. Run all simulation scenarios -----------------------------------
sim_grid <- expand_grid(
  r_true   = r_vals,
  R_true   = R_vals,
  noise_sd = noise_sd_vec,
  dt       = dt_vec,
  duration = dur_vec,
  rep_id   = seq_len(n_reps)
)

sim_results <- sim_grid %>%
  pmap_dfr(function(r_true, R_true, noise_sd, dt, duration, rep_id) {
    simulate_and_fit(
      r_true      = r_true,
      R_true      = R_true,
      N0_true     = N0_true,
      O2_0_true   = O2_0_true,
      O2_ref_true = O2_ref_true,
      noise_sd    = noise_sd,
      dt          = dt,
      duration    = duration
    )
  })

# Save raw simulation results (TABLES)
write_csv(
  sim_results,
  file.path(tables_dir, "synthetic_parameter_recovery_results.csv")
)

# ---- 5. Summarise recovery (bias, RMSE etc.) ---------------------------
sim_summary <- sim_results %>%
  filter(converged) %>%
  mutate(
    r_error   = r_est - r_true,
    R_error   = R_est - R_true,
    r_rel_err = (r_est - r_true) / r_true,
    R_rel_err = (R_est - R_true) / R_true
  ) %>%
  group_by(noise_sd, dt, duration) %>%
  summarise(
    n_fits     = n(),
    r_bias     = mean(r_error, na.rm = TRUE),
    r_rmse     = sqrt(mean(r_error^2, na.rm = TRUE)),
    R_bias     = mean(R_error, na.rm = TRUE),
    R_rmse     = sqrt(mean(R_error^2, na.rm = TRUE)),
    r_rel_bias = mean(r_rel_err, na.rm = TRUE),
    R_rel_bias = mean(R_rel_err, na.rm = TRUE),
    mean_R2    = mean(pseudo_R2, na.rm = TRUE),
    .groups    = "drop"
  )

write_csv(
  sim_summary,
  file.path(tables_dir, "Table_S2_synthetic_parameter_recovery_summary.csv")
)

# ---- 6. Plots: true vs estimated r and R (for main scenario) ----------
# Choose a “realistic” scenario to plot, e.g. dt = 1, duration = 180
plot_subset <- sim_results %>%
  filter(converged, dt == 1, duration == 180)

# (A) r
plot_r <- ggplot(plot_subset,
                 aes(x = r_true, y = r_est, colour = factor(noise_sd))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.6) +
  labs(
    x = expression("True " * r ~ (min^{-1})),
    y = expression("Estimated " * r ~ (min^{-1})),
    colour = "Noise SD"
  ) +
  theme_classic(base_size = 12)

# (B) R
plot_R <- ggplot(plot_subset,
                 aes(x = R_true, y = R_est, colour = factor(noise_sd))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.6) +
  labs(
    x = expression("True " * R ~ (mg~O[2]~cell^{-1}~min^{-1})),
    y = expression("Estimated " * R ~ (mg~O[2]~cell^{-1}~min^{-1})),
    colour = "Noise SD"
  ) +
  theme_classic(base_size = 12)

# Combine into a single figure (r left, R right)
synthetic_fig <- plot_r + plot_R + plot_layout(nrow = 1)

ggsave(
  file.path(plots_dir, "Fig_S2_synthetic_param_recovery_combined.pdf"),
  synthetic_fig, width = 9, height = 4, dpi = 300
)
ggsave(
  file.path(plots_dir, "Fig_S2_synthetic_param_recovery_combined.png"),
  synthetic_fig, width = 9, height = 4, dpi = 300
)

############################################################
# End of script
############################################################
