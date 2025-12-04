################################################################################
# README – O₂-based growth & respiration pipeline
#
# This script:
#   1) Trims raw dissolved-oxygen time series to the main exponential-decline phase.
#   2) Fits the normalised O₂ model with N₀ fixed to 1:
#        O2_norm(t) = O2_0 + (resp_tot / r) * (1 - exp(r * t))
#   3) Reconstructs N₀ at O₂ start from inoculation densities and the fitted r.
#   4) Converts resp_tot to per-cell respiration rates R (mg O₂ cell⁻¹ min⁻¹).
#   5) Produces diagnostics, comparison plots (O₂ vs OD vs FC), and summary tables.
#
# Main inputs (in /data): Oxygen_Data_Long.csv, Inoculation_Density.csv, OD_r_FC_r.csv
# Main outputs: trimmed O₂ CSVs, model fits (oxygen_model_results*.csv),
#               growth-rate comparison plots, residual diagnostics, and Bland–Altman plots.
################################################################################

# ── Libraries ─────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)
  library(zoo)
  library(minpack.lm)  # Levenberg–Marquardt nls
  library(patchwork)   # for wrap_plots()
  library(ggsignif)
  library(lme4)
  library(multcomp)
  library(gridExtra)
  library(lmerTest)
  library(glue)
})

################################################################################
# 0. PROJECT PATHS & INOCULATION CSV
################################################################################

data_dir   <- "data"
plots_dir  <- "plots"
tables_dir <- "Tables"

dir.create(data_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# New inoculation density CSV (you provide this)
# Expected columns: Taxon, Replicate, inoc_cells_per_uL
INOC_CSV <- file.path(data_dir, "Inoculation_Density.csv")

################################################################################
# 1. TRIMMING OXYGEN TIME SERIES
################################################################################

# ── Tunables (trimming aggressiveness) ────────────────────────────────
SPLINE_SPAR    <- 0.4     # higher = smoother (often trims earlier)
RUN_LEN        <- 10      # consecutive calm points required for plateau
REL_PROP       <- 0.008   # fraction of max post-peak down-slope allowed
ABS_THR        <- 0.0003  # absolute slope threshold floor
WINDOW_LEN     <- 4       # rolling range window length for flatness
LEVEL_DELTA    <- 0.0025  # allowed vertical wiggle (range) within window
FLAT_RANGE_OK  <- 0.05    # skip series if total range below this

INPUT_CSV      <- file.path(data_dir,   "Oxygen_Data_Long.csv")
TRIMMED_CSV    <- file.path(tables_dir, "Oxygen_Data_Trimmed.csv")
FILTERED_CSV   <- file.path(tables_dir, "Oxygen_Data_Filtered.csv")
SKIPPED_CSV    <- file.path(tables_dir, "Skipped_Series_Log.csv")
DIAG_PDF       <- file.path(plots_dir,  "oxygen_trimming_diagnostics.pdf")

# ── Helper functions ──────────────────────────────────────────────────
find_plateau <- function(o2_vec, peak_idx,
                         run_len = 2, rel_prop = 0.01, abs_thr = 0.0003,
                         window_len = 4, level_delta = 0.0025) {
  n <- length(o2_vec); if (peak_idx >= n) return(n)
  slopes <- diff(o2_vec)
  max_down <- abs(min(slopes[peak_idx:(n - 1)], na.rm = TRUE))
  slope_thr <- max(abs_thr, rel_prop * max_down)
  
  rng <- rep(NA_real_, n)
  if (n >= window_len) {
    rng[1:(n - window_len + 1)] <-
      zoo::rollapply(o2_vec, window_len,
                     FUN = function(x) max(x) - min(x), align = "left")
  }
  
  for (i in seq(peak_idx, n - run_len)) {
    if (all(abs(slopes[i:(i + run_len - 1)]) <= slope_thr, na.rm = TRUE) &&
        all(rng[i:(i + run_len - 1)] <= level_delta, na.rm = TRUE)) {
      return(i)
    }
  }
  n
}

find_second_inflection <- function(o2_fit, idx_peak, idx_plate) {
  slopes <- diff(o2_fit); curves <- diff(slopes)
  post_pk <- seq(idx_peak + 1, length(curves))
  if (length(post_pk) == 0) return(NA_integer_)
  steepest <- which.min(slopes[post_pk]) + idx_peak
  win <- seq(steepest + 1, idx_plate - 1)
  win <- win[!is.na(curves[win])]
  if (length(win) < 3) return(NA_integer_)
  win[which.max(curves[win])]
}

# ── Load raw oxygen data ──────────────────────────────────────────────
raw <- read_csv(INPUT_CSV, show_col_types = FALSE) %>%
  dplyr::mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate),
    series_id = paste(Taxon, "Rep=", Replicate, sep = " ")
  )

series_ids  <- unique(raw$series_id)
trimmed_lst <- vector("list", length(series_ids))
names(trimmed_lst) <- series_ids
skipped_log <- tibble(series_id = character(), reason = character())

# ── Main trimming loop ────────────────────────────────────────────────
for (sid in series_ids) {
  df <- raw %>%
    dplyr::filter(series_id == sid) %>%
    dplyr::arrange(Time)
  
  if (nrow(df) < 5) {
    skipped_log <- dplyr::add_row(skipped_log, series_id = sid, reason = "Too few rows")
    next
  }
  
  # Smoothing used only for detection; not plotted
  df <- df %>%
    dplyr::mutate(O2_fit = suppressWarnings(
      predict(smooth.spline(Time, Oxygen, spar = SPLINE_SPAR), x = Time)$y
    ))
  
  if (all(is.na(df$O2_fit))) {
    skipped_log <- dplyr::add_row(skipped_log, series_id = sid, reason = "Spline failed")
    next
  }
  
  if ((max(df$O2_fit, na.rm = TRUE) - min(df$O2_fit, na.rm = TRUE)) < FLAT_RANGE_OK) {
    skipped_log <- dplyr::add_row(skipped_log, series_id = sid, reason = "Too flat")
    next
  }
  
  # Peak on raw Oxygen (true peak)
  idx_peak <- which.max(df$Oxygen)
  
  # Plateau & second inflection on smoothed curve
  idx_plate  <- min(
    find_plateau(df$O2_fit, idx_peak,
                 run_len = RUN_LEN, rel_prop = REL_PROP, abs_thr = ABS_THR,
                 window_len = WINDOW_LEN, level_delta = LEVEL_DELTA),
    nrow(df)
  )
  idx_second <- find_second_inflection(df$O2_fit, idx_peak, idx_plate)
  
  if (is.na(idx_second) || idx_second <= idx_peak || idx_second > nrow(df)) {
    idx_end <- nrow(df); used_full <- TRUE
  } else {
    idx_end <- idx_second; used_full <- FALSE
  }
  
  trimmed_lst[[sid]] <- df[idx_peak:idx_end, ] %>%
    dplyr::mutate(
      peak_time           = Time[idx_peak],
      plateau_time        = if (idx_plate <= nrow(df)) Time[idx_plate] else NA_real_,
      second_inflect_time = if (!is.na(idx_second)) Time[idx_second] else NA_real_,
      used_full_series    = used_full,
      series_id           = sid
    )
}

# ── Save trimmed and filtered data ────────────────────────────────────
trimmed <- dplyr::bind_rows(trimmed_lst)
write_csv(trimmed,  TRIMMED_CSV)

filtered <- trimmed %>%
  dplyr::select(Taxon, Replicate, Time, Oxygen)

write_csv(filtered, FILTERED_CSV)
write_csv(skipped_log, SKIPPED_CSV)

# ── Diagnostics PDF (highlighted trimmed zone) ────────────────────────
pdf(DIAG_PDF, 7, 5)

for (sid in names(trimmed_lst)) {
  df_trim <- trimmed_lst[[sid]]
  if (is.null(df_trim) || nrow(df_trim) == 0) next
  
  raw_df <- raw %>%
    dplyr::filter(series_id == sid) %>%
    dplyr::arrange(Time)
  
  xmin_val <- suppressWarnings(min(df_trim$Time, na.rm = TRUE))
  xmax_val <- suppressWarnings(max(df_trim$Time, na.rm = TRUE))
  if (!is.finite(xmin_val) || !is.finite(xmax_val)) next
  
  rect_df <- tibble(
    xmin = xmin_val,
    xmax = xmax_val,
    ymin = -Inf, ymax = Inf
  )
  
  p <- ggplot(raw_df, aes(Time, Oxygen)) +
    geom_rect(
      data = rect_df, inherit.aes = FALSE,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "orange", alpha = 0.08
    ) +
    geom_line(colour = "grey60", linewidth = 0.8) +
    labs(
      title = sid,
      subtitle = "Trimmed zone (orange)",
      x = "Time (min)", y = "O₂ (mg L⁻¹)"
    ) +
    theme_classic(base_size = 11)
  
  print(p)
}

dev.off()


################################################################################
# 2. O2‐BASED GROWTH & TOTAL RESPIRATION (N0 FIXED TO 1 IN MODEL)
################################################################################

## ──────────────────── 1. QC / Fit thresholds ───────────────────────────── ##
REL_SE_THRESHOLD <- 0.15
PVAL_THRESHOLD   <- 0.001
R2_THRESHOLD     <- 0.90
MAX_RESID_RANGE  <- 0.20
AIC_IMPROVEMENT  <- 10
MAPE_MAX         <- 0.15

## ──────────────────── 2. Timing info ───────────────────────────────────── ##
INOC_DELAY_MIN <- 45   # minutes between inoculation and first O2 reading

## ───────────────────────── 3.  Load trimmed oxygen data ─────────────────── ##
oxygen_data <- read_csv(FILTERED_CSV, show_col_types = FALSE) %>%
  dplyr::mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate)
  )

## ───────────────────────── 4.  Model function (N0 fixed to 1) ───────────── ##
# Oxygen_norm(t) = O2_0 + (resp_tot / r) * (1 - exp(r * t))
# resp_tot is a TOTAL respiration scaling parameter in normalised units;
# dimensional per-cell rates (mg O2 cell^-1 min^-1) are recovered later using
# inoculation density, r, and O2_ref.
resp_model <- function(r, resp_tot, t, O2_0) {
  O2_0 + (resp_tot / r) * (1 - exp(r * t))
}

## ───────────────────────── 5.  Plot theme ───────────────────────────────── ##
isme_theme <- function() {
  theme_classic(base_size = 14) +
    theme(
      axis.title = element_text(size = 16),
      axis.text  = element_text(size = 14),
      axis.line  = element_line(linewidth = .8),
      axis.ticks = element_line(linewidth = .6),
      legend.position = "none"
    )
}

## ───────────────────────── 6.  Prepare outputs ──────────────────────────── ##
results <- tibble(
  Taxon                = character(),
  Replicate            = character(),
  r_per_minute         = numeric(),
  r_ci_lower           = numeric(),
  r_ci_upper           = numeric(),
  r_per_hour           = numeric(),
  resp_tot             = numeric(),   # total respiration scaling (model param)
  N0_cells_per_L       = numeric(),   # will be filled later
  resp_rate            = numeric(),   # per-cell respiration (mg O2 cell^-1 min^-1)
  resp_ci_lower        = numeric(),   # left NA
  resp_ci_upper        = numeric(),   # left NA
  O2_0                 = numeric(),
  O2_ref               = numeric(),   # mean initial O2 (mg L^-1) used for normalisation
  AICc                 = numeric(),
  lnO2_change_per_min  = numeric(),
  pseudo_R2            = numeric(),
  fit_ok               = logical()
)

plots_list <- list()

resid_all <- tibble(
  Taxon     = character(),
  Replicate = character(),
  Time0     = numeric(),
  Fitted    = numeric(),
  Residual  = numeric()
)

## ───────────────────────── 7.  Fitting loop ─────────────────────────────── ##
combos <- dplyr::group_keys(dplyr::group_by(oxygen_data, Taxon, Replicate))

for (i in seq_len(nrow(combos))) {
  cfg <- combos[i, ]
  Tax <- cfg$Taxon
  Rep <- cfg$Replicate
  
  df <- oxygen_data %>%
    dplyr::filter(Taxon == Tax, Replicate == Rep) %>%
    dplyr::arrange(Time)
  if (nrow(df) < 5) next
  
  # derivative-based trimming of early flat bits
  df <- df %>%
    dplyr::mutate(
      dO2   = c(NA, diff(Oxygen)),
      dt    = c(NA, diff(Time)),
      dO2dt = dO2 / dt,
      sm    = zoo::rollmean(dO2dt, 3, fill = NA, align = "right")
    )
  idx <- which(df$sm < -1e-7)[1]
  idx <- if (is.na(idx) || idx > nrow(df) - 2) 1 else idx + 5
  df  <- df[idx:nrow(df), ] %>%
    dplyr::mutate(Time0 = Time - min(Time, na.rm = TRUE))
  
  # Normalise O2 for fitting
  O0 <- mean(head(df$Oxygen, 3), na.rm = TRUE)   # O2_ref in mg L^-1
  df <- df %>% dplyr::mutate(Oxygen_norm = Oxygen / O0)
  
  # Starting values for r and resp_tot
  r_start <- {
    seg    <- head(df, max(3, floor(.3 * nrow(df))))
    slopes <- abs(diff(log(pmax(seg$Oxygen_norm, 1e-6))) / diff(seg$Time0))
    pmin(pmax(max(slopes, na.rm = TRUE), 1e-4), 5e-2)
  }
  resp_tot_start <- {
    slope <- suppressWarnings(min(diff(df$Oxygen_norm) / diff(df$Time0), na.rm = TRUE))
    if (!is.finite(slope)) slope <- -1e-3
    k_guess <- abs(slope)            # since dO2_norm/dt|0 ≈ -resp_tot
    k_guess <- pmin(pmax(k_guess, 1e-5), 0.1)
    k_guess
  }
  
  fit <- tryCatch(
    nlsLM(
      Oxygen_norm ~ resp_model(r, resp_tot, Time0, O2_0),
      data    = df,
      start   = list(r = r_start, resp_tot = resp_tot_start, O2_0 = 1),
      lower   = c(r = 1e-4,  resp_tot = 1e-5, O2_0 = .8),
      upper   = c(r = .1,    resp_tot = 0.5,  O2_0 = 1.2),
      control = nls.lm.control(maxiter = 300)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) next
  
  # Outlier removal based on residuals and refit
  pred    <- predict(fit, df)
  keep_ix <- abs(df$Oxygen_norm - pred) < 2 * sd(df$Oxygen_norm - pred, na.rm = TRUE)
  df_kept <- df[keep_ix, ]
  fit     <- tryCatch(update(fit, data = df_kept), error = function(e) fit)
  
  pars  <- coef(summary(fit))
  r_est <- pars["r", "Estimate"]
  K_est <- pars["resp_tot", "Estimate"]   # total respiration scaling
  
  pseudo_R2 <- 1 - sum(residuals(fit)^2) /
    sum((df_kept$Oxygen_norm - mean(df_kept$Oxygen_norm))^2)
  
  # ── QC: based on fitted parameters (r, resp_tot) and residuals ───────
  fit_ok <- all(
    abs(pars[, "Std. Error"] / pars[, "Estimate"]) < REL_SE_THRESHOLD,
    pars[, "Pr(>|t|)"] < PVAL_THRESHOLD,
    pseudo_R2 >= R2_THRESHOLD,
    diff(range(residuals(fit), na.rm = TRUE)) < MAX_RESID_RANGE,
    is.finite(K_est),
    K_est > 0, K_est < 0.5,
    dplyr::between(r_est, 1e-4, 0.1),
    AIC(lm(Oxygen_norm ~ 1, data = df_kept)) - AIC(fit) >= AIC_IMPROVEMENT,
    mean(abs(residuals(fit) / df_kept$Oxygen_norm)) < MAPE_MAX
  )
  
  # ── Wald CIs for r (respiration CIs left NA) ─────────────────────────
  se_r  <- pars["r", "Std. Error"]
  z_975 <- 1.96
  
  r_ci_lower    <- if (is.finite(se_r)) r_est - z_975 * se_r else NA_real_
  r_ci_upper    <- if (is.finite(se_r)) r_est + z_975 * se_r else NA_real_
  resp_ci_lower <- NA_real_
  resp_ci_upper <- NA_real_
  
  df_kept <- df_kept %>% dplyr::mutate(Pred = predict(fit, df_kept))
  
  resid_all <- dplyr::bind_rows(
    resid_all,
    tibble(
      Taxon     = Tax,
      Replicate = Rep,
      Time0     = df_kept$Time0,
      Fitted    = df_kept$Pred,
      Residual  = df_kept$Oxygen_norm - df_kept$Pred
    )
  )
  
  results <- results %>%
    dplyr::add_row(
      Taxon                = Tax,
      Replicate            = Rep,
      r_per_minute         = r_est,
      r_ci_lower           = r_ci_lower,
      r_ci_upper           = r_ci_upper,
      r_per_hour           = r_est * 60,
      resp_tot             = K_est,
      N0_cells_per_L       = NA_real_,  # to be filled later
      resp_rate            = NA_real_,  # to be filled later (mg O2 cell^-1 min^-1)
      resp_ci_lower        = resp_ci_lower,
      resp_ci_upper        = resp_ci_upper,
      O2_0                 = coef(fit)["O2_0"],
      O2_ref               = O0,        # store normalisation constant (mg L^-1)
      AICc                 = AIC(fit),
      lnO2_change_per_min  = (log(dplyr::last(df_kept$Oxygen_norm)) -
                                log(dplyr::first(df_kept$Oxygen_norm))) /
        (dplyr::last(df_kept$Time0) - dplyr::first(df_kept$Time0)),
      pseudo_R2            = pseudo_R2,
      fit_ok               = fit_ok
    )
  
  plot_key <- paste(Tax, Rep, sep = "_")
  plots_list[[plot_key]] <-
    ggplot(df_kept, aes(Time0, Oxygen_norm)) +
    geom_point(size = 2.5) +
    geom_line(aes(y = Pred), linewidth = .8, colour = "red") +
    labs(
      title = plot_key,
      x     = expression(Time~"(min)"),
      y     = expression(Normalised~O[2])
    ) +
    isme_theme()
}

# Save raw results (without N0 + resp_rate filled)
write_csv(results, file.path(tables_dir, "oxygen_model_results_raw.csv"))

if (length(plots_list) > 0) {
  pdf(file.path(plots_dir, "oxygen_dynamics_all_models.pdf"), width = 14, height = 10)
  print(wrap_plots(plots_list))
  dev.off()
  
  pdf(file.path(plots_dir, "oxygen_dynamics_fullsize_per_page.pdf"), width = 8, height = 6)
  purrr::walk(plots_list, print)
  dev.off()
}

## ─────────────────── 8. Residual diagnostics: Supp Fig S1 ─────────────── ##
if (nrow(resid_all) > 0) {
  resid_all <- resid_all %>%
    dplyr::mutate(
      Taxon     = factor(Taxon),
      Replicate = factor(Replicate)
    )
  
  p_resid_time <- ggplot(resid_all, aes(x = Time0, y = Residual, colour = Replicate)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_point(alpha = 0.7, size = 0.8) +
    facet_wrap(~ Taxon, scales = "free_x") +
    labs(
      x = "Time (min, re-zeroed)",
      y = "Residual (normalised O\u2082)"
    ) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.text      = element_text(face = "italic"),
      axis.title      = element_text(size = 11),
      axis.text       = element_text(size = 9),
      legend.title    = element_text(size = 10),
      legend.text     = element_text(size = 9)
    )
  
  p_resid_fit <- ggplot(resid_all, aes(x = Fitted, y = Residual, colour = Replicate)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_point(alpha = 0.7, size = 0.8) +
    facet_wrap(~ Taxon, scales = "free_x") +
    labs(
      x = "Fitted normalised O\u2082",
      y = "Residual (normalised O\u2082)"
    ) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "none",
      strip.text      = element_text(face = "italic"),
      axis.title      = element_text(size = 11),
      axis.text       = element_text(size = 9)
    )
  
  S1_residuals <- p_resid_time / p_resid_fit +
    patchwork::plot_layout(ncol = 1, heights = c(1, 1)) +
    patchwork::plot_annotation(
      title      = "Supplementary Fig. S1. Residual diagnostics for O\u2082-based growth–respiration fits",
      tag_levels = "A"
    )
  
  ggsave(
    filename = file.path(plots_dir, "Supp_Fig_S1_residuals.png"),
    plot     = S1_residuals,
    width    = 12,
    height   = 10,
    dpi      = 300
  )
}


################################################################################
# 3. N0 FROM INOCULATION DENSITY CSV + PER-CELL RESPIRATION
################################################################################

# 3.1 Extract O2-based growth rates
oxygen_r <- results %>%
  dplyr::mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate)
  ) %>%
  dplyr::select(Taxon, Replicate, r_per_minute, fit_ok)

# 3.2 Load inoculation density (cells/µL at t = inoculation)
# Expected columns in INOC_CSV: Taxon, Replicate, inoc_cells_per_uL
inoc_input <- read_csv(INOC_CSV, show_col_types = FALSE) %>%
  dplyr::mutate(
    Taxon             = as.character(Taxon),
    Replicate         = as.character(Replicate),
    inoc_cells_per_uL = as.numeric(inoc_cells_per_uL)
  )

# 3.3 Join and compute N0 at start of O2 measurements
# Exponential growth: cells(t) = cells_inoc * exp(r * t)
# So N0 at O2 start (t = INOC_DELAY_MIN) is:
#   N0_O2start = inoc_cells_per_uL * exp(r * INOC_DELAY_MIN)
startO2_density <- inoc_input %>%
  dplyr::left_join(oxygen_r, by = c("Taxon", "Replicate")) %>%
  dplyr::mutate(
    N0_O2start_cells_per_uL = inoc_cells_per_uL * exp(r_per_minute * INOC_DELAY_MIN),
    N0_O2start_cells_per_mL = N0_O2start_cells_per_uL * 1e3,
    N0_O2start_cells_per_L  = N0_O2start_cells_per_uL * 1e6
  )

# Optional QC
startO2_density_qc <- startO2_density %>% dplyr::filter(fit_ok)

startO2_summary <- startO2_density_qc %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(
    mean_inoc_cells_per_uL       = mean(inoc_cells_per_uL, na.rm = TRUE),
    sd_inoc_cells_per_uL         = sd(inoc_cells_per_uL,   na.rm = TRUE),
    mean_N0_O2start_cells_per_uL = mean(N0_O2start_cells_per_uL, na.rm = TRUE),
    sd_N0_O2start_cells_per_uL   = sd(N0_O2start_cells_per_uL,   na.rm = TRUE),
    n                            = dplyr::n(),
    .groups                      = "drop"
  )

write_csv(startO2_density,     file.path(tables_dir, "O2start_density_estimates_all.csv"))
write_csv(startO2_density_qc,  file.path(tables_dir, "O2start_density_estimates_QC.csv"))
write_csv(startO2_summary,     file.path(tables_dir, "O2start_density_summary_by_taxon.csv"))

# 3.4 Respiration rate per cell (using inoculation-based N0 at O2 start)
# Here we convert resp_tot (normalised units) to dimensional mg O2 cell^-1 min^-1:
#   R = resp_tot * O2_ref / N0
resp_per_cell <- results %>%
  dplyr::mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate)
  ) %>%
  dplyr::left_join(
    startO2_density %>%
      dplyr::mutate(
        Taxon     = as.character(Taxon),
        Replicate = as.character(Replicate)
      ) %>%
      dplyr::select(Taxon, Replicate, N0_O2start_cells_per_L),
    by = c("Taxon", "Replicate")
  ) %>%
  dplyr::mutate(
    respiration_per_cell = dplyr::if_else(
      is.finite(N0_O2start_cells_per_L) & N0_O2start_cells_per_L > 0 &
        is.finite(O2_ref) & O2_ref > 0,
      resp_tot * O2_ref / N0_O2start_cells_per_L,  # mg O2 cell^-1 min^-1
      NA_real_
    )
  )

resp_per_cell_qc <- resp_per_cell %>%
  dplyr::filter(
    fit_ok,
    is.finite(respiration_per_cell),
    respiration_per_cell > 0
  )

resp_per_cell_summary <- resp_per_cell_qc %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(
    mean_respiration_per_cell = mean(respiration_per_cell, na.rm = TRUE),
    sd_respiration_per_cell   = sd(respiration_per_cell,   na.rm = TRUE),
    n                         = dplyr::n(),
    .groups                   = "drop"
  )

write_csv(resp_per_cell,        file.path(tables_dir, "respiration_rate_per_cell_all.csv"))
write_csv(resp_per_cell_qc,     file.path(tables_dir, "respiration_rate_per_cell_QC.csv"))
write_csv(resp_per_cell_summary,file.path(tables_dir, "respiration_rate_per_cell_summary_by_taxon.csv"))

# 3.5 Update `results` with N0 and per-cell respiration, and save final table
results <- results %>%
  dplyr::mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate)
  ) %>%
  dplyr::left_join(
    resp_per_cell %>%
      dplyr::select(
        Taxon, Replicate,
        N0_cells_per_L = N0_O2start_cells_per_L,
        resp_rate      = respiration_per_cell
      ),
    by = c("Taxon", "Replicate")
  )

write_csv(results, file.path(tables_dir, "oxygen_model_results.csv"))


################################################################################
# 4. 4-PANEL FACET PLOT (LABELS UNDER CURVE, USING r & R)
################################################################################

TEXT_SIZE    <- 3
X_ANCHOR     <- 0.00001
X_MARGIN     <- 0.05
Y_GAP_FRAC   <- 0.80
Y_MARGIN     <- 0.05
ONE_LINE     <- FALSE

sci_pm <- function(x, digits = 2) {
  out <- rep(NA_character_, length(x))
  ok  <- is.finite(x)
  if (any(ok)) {
    s <- formatC(x[ok], format = "e", digits = digits - 1)
    m <- sub("e[+-].*$", "", s)
    e <- sub("^.*e([+-]?)([0-9]+)$", "\\1\\2", s)
    e <- sub("^\\+", "", e)
    out[ok] <- sprintf("%s%%*%%10^{%s}", m, e)
  }
  out
}

selected_combos <- tibble::tribble(
  ~Taxon,          ~Replicate,
  "Bacillus",      "R4",
  "Burkholderia",  "R2",
  "Arthrobacter",  "R2",
  "Yersinia",      "R1"
)
letters_vec <- letters[seq_len(nrow(selected_combos))]

# Rebuild panel data from plots_list
facet_data <- purrr::map2_dfr(
  selected_combos$Taxon, selected_combos$Replicate,
  ~{
    key <- paste(.x, .y, sep = "_")
    p   <- plots_list[[key]]
    if (is.null(p)) return(NULL)
    
    pts <- ggplot2::layer_data(p, 1) %>%
      dplyr::select(x, y) %>%
      dplyr::rename(Time = x, Oxygen = y)
    
    fit <- ggplot2::layer_data(p, 2) %>%
      dplyr::select(x, y) %>%
      dplyr::rename(Time = x, Predicted_O2 = y)
    
    dplyr::full_join(pts, fit, by = "Time") %>%
      dplyr::mutate(Taxon = .x, Replicate = .y, series_id = key)
  }
)

stopifnot(all(c("Time","Oxygen","Predicted_O2","Taxon","Replicate") %in% names(facet_data)))

## --- 4-panel: build resp_rate fresh, don't assume it's in `results` ---

# 1) Make sure we have per-cell respiration available for labels (mg O2 cell^-1 min^-1)
results_for_labels <- results %>%
  dplyr::mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate)
  ) %>%
  dplyr::left_join(
    startO2_density %>%
      dplyr::mutate(
        Taxon     = as.character(Taxon),
        Replicate = as.character(Replicate)
      ) %>%
      dplyr::select(Taxon, Replicate, N0_O2start_cells_per_L),
    by = c("Taxon", "Replicate")
  ) %>%
  dplyr::mutate(
    resp_rate = dplyr::if_else(
      is.finite(N0_O2start_cells_per_L) & N0_O2start_cells_per_L > 0 &
        is.finite(O2_ref) & O2_ref > 0,
      resp_tot * O2_ref / N0_O2start_cells_per_L,
      NA_real_
    )
  )

# 2) Build label info from results_for_labels (which now DEFINITELY has resp_rate)
annot_info <- results_for_labels %>%
  dplyr::inner_join(selected_combos, by = c("Taxon","Replicate")) %>%
  dplyr::mutate(
    FacetLabel = paste0(
      "(", letters_vec[
        match(paste(Taxon, Replicate),
              paste(selected_combos$Taxon, selected_combos$Replicate))
      ], ")~italic('", Taxon, "')"
    ),
    label_text = if (ONE_LINE) {
      paste0(
        "italic(r)==", sprintf("%.2f", r_per_hour), "~h^{-1}*','~~",
        "italic(R)==", sci_pm(resp_rate, digits = 2),
        "~mg~O[2]~cell^{-1}~min^{-1}"
      )
    } else {
      paste0(
        "atop(",
        "italic(r)==", sprintf("%.2f", r_per_hour), "~h^{-1},",
        "italic(R)==", sci_pm(resp_rate, digits = 2),
        "~mg~O[2]~cell^{-1}~min^{-1})"
      )
    }
  ) %>%
  dplyr::select(Taxon, Replicate, FacetLabel, label_text)

facet_data <- facet_data %>%
  dplyr::inner_join(annot_info, by = c("Taxon","Replicate"))

# Label positions
label_positions <- facet_data %>%
  dplyr::group_by(FacetLabel) %>%
  dplyr::group_modify(~{
    df <- .x %>%
      dplyr::arrange(Time) %>%
      dplyr::distinct(Time, .keep_all = TRUE)
    
    xmin <- min(df$Time,    na.rm = TRUE); xmax <- max(df$Time,    na.rm = TRUE)
    ymin <- min(df$Oxygen,  na.rm = TRUE); ymax <- max(df$Oxygen,  na.rm = TRUE)
    xr   <- xmax - xmin;    yr   <- ymax - ymin
    
    x_pos <- xmin + X_ANCHOR * xr
    x_pos <- max(xmin + X_MARGIN * xr, min(xmax - X_MARGIN * xr, x_pos))
    
    ok <- is.finite(df$Predicted_O2)
    curveY <- if (sum(ok) >= 2) approx(df$Time[ok], df$Predicted_O2[ok],
                                       xout = x_pos, rule = 2)$y
    else ymin + 0.90 * yr
    
    gap   <- Y_GAP_FRAC * yr
    y_raw <- curveY - gap
    y_pos <- max(ymin + Y_MARGIN * yr, min(ymax - Y_MARGIN * yr, y_raw))
    
    tibble(x_pos = x_pos, y_pos = y_pos)
  }) %>%
  dplyr::ungroup()

annotations <- facet_data %>%
  dplyr::distinct(FacetLabel, label_text) %>%
  dplyr::left_join(label_positions, by = "FacetLabel")

facet_plot <- ggplot(facet_data, aes(x = Time)) +
  geom_point(aes(y = Oxygen), size = 2.2, alpha = 0.85) +
  geom_line(aes(y = Predicted_O2), linewidth = 1.2) +
  facet_wrap(~ FacetLabel, nrow = 1, labeller = label_parsed) +
  geom_text(
    data = annotations,
    aes(x = x_pos, y = y_pos, label = label_text),
    inherit.aes = FALSE, parse = TRUE, hjust = 0, vjust = 1.05, size = TEXT_SIZE
  ) +
  coord_cartesian(clip = "on") +
  labs(
    x = "Time (minutes)",
    y = expression("Normalised Oxygen (" * O[2] / O[2*","*0] * ")")
  ) +
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 16, face = "italic"),
    axis.title       = element_text(size = 17),
    axis.text        = element_text(size = 14),
    panel.spacing    = unit(1.2, "lines"),
    plot.margin      = margin(12, 20, 12, 20),
    axis.line        = element_line(linewidth = 0.7),
    axis.ticks       = element_line(linewidth = 0.6)
  )

ggsave(
  filename = file.path(plots_dir, "Fig_2_oxygen_dynamics_facet4_no_overflow.pdf"),
  plot     = facet_plot,
  width    = 14,
  height   = 4.5
)


################################################################################
# 5. Supplementary (NORMALISED): colour by Replicate, facet by TaxonFull
################################################################################

all_fits_df_norm <- purrr::imap_dfr(
  plots_list,
  ~{
    key <- .y
    p   <- .x
    
    pts <- tryCatch(
      ggplot2::layer_data(p, 1) %>%
        dplyr::select(x, y) %>%
        dplyr::rename(Time = x, Oxygen = y),
      error = function(e) NULL
    )
    fit <- tryCatch(
      ggplot2::layer_data(p, 2) %>%
        dplyr::select(x, y) %>%
        dplyr::rename(Time = x, Predicted = y),
      error = function(e) NULL
    )
    if (is.null(pts) || is.null(fit)) return(NULL)
    
    df <- dplyr::full_join(pts, fit, by = "Time") %>%
      dplyr::arrange(Time)
    
    m1 <- stringr::str_match(key, "^(.*)_(R[0-9]+)$")
    TaxonFull <- m1[,2]
    Replicate <- m1[,3]
    
    O0 <- mean(head(df$Oxygen, 3), na.rm = TRUE)
    if (is.finite(O0) && O0 > 0) {
      df <- df %>%
        dplyr::mutate(
          Oxygen_n = Oxygen   / O0,
          Pred_n   = Predicted / O0
        )
    } else {
      df <- df %>%
        dplyr::mutate(
          Oxygen_n = Oxygen,
          Pred_n   = Predicted
        )
    }
    
    tibble(
      TaxonFull = TaxonFull,
      Replicate = Replicate,
      Time      = df$Time,
      Oxygen_n  = df$Oxygen_n,
      Pred_n    = df$Pred_n
    )
  }
) %>%
  dplyr::filter(!is.na(TaxonFull), !is.na(Replicate)) %>%
  dplyr::mutate(
    Replicate = factor(Replicate, levels = paste0("R", 1:5))
  )

supp_plot_norm_rep <- ggplot(all_fits_df_norm, aes(x = Time, y = Oxygen_n)) +
  geom_point(color = "grey60", size = 1, alpha = 0.55) +
  geom_line(
    data = all_fits_df_norm %>% dplyr::arrange(TaxonFull, Replicate, Time),
    aes(y = Pred_n, color = Replicate, group = Replicate),
    linewidth = 0.9,
    na.rm = TRUE
  ) +
  facet_wrap(~ TaxonFull, scales = "free_y") +
  labs(
    title = "Supplementary: Normalised O₂ Dynamics (colour = Replicate)",
    x = "Time (minutes)",
    y = expression("Normalised Oxygen ("*O[2]*"/"*O[2][0]*")"),
    color = "Replicate"
  ) +
  scale_color_brewer(palette = "Dark2", drop = FALSE) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold"),
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 10)
  )

ggsave(
  filename = file.path(plots_dir, "supp_Fig_3_oxygen_all_replicates_NORMALISED_by_replicate.pdf"),
  plot     = supp_plot_norm_rep,
  width    = 14,
  height   = 10
)


################################################################################
# 6. BOXPLOT + MIXED-EFFECTS ANOVA (O2 vs OD600 vs FC growth rates)
################################################################################

od_fc_data <- read_csv(file.path(data_dir, "OD_r_FC_r.csv"), show_col_types = FALSE) %>%
  dplyr::mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate),
    Duration  = Time
  )

OD_FC <- od_fc_data %>%
  dplyr::mutate(
    r_OD600 = (log(OD_Final) - log(OD_Initial)) / Duration,
    r_FC    = (log(FC_Final) - log(FC_Initial)) / Duration
  ) %>%
  dplyr::select(Taxon, Replicate, Duration, r_OD600, r_FC)

growth <- results %>%
  dplyr::filter(fit_ok) %>%
  dplyr::select(Taxon, Replicate, r_per_minute) %>%
  dplyr::rename(r_O2 = r_per_minute) %>%
  dplyr::left_join(OD_FC, by = c("Taxon", "Replicate")) %>%
  dplyr::select(Taxon, Replicate, r_O2, r_OD600, r_FC)

growth <- growth %>%
  dplyr::mutate(
    doubling_time_O2    = log(2) / r_O2,
    doubling_time_OD600 = log(2) / r_OD600,
    doubling_time_FC    = log(2) / r_FC
  )

write_csv(growth, file.path(tables_dir, "growth_rates_combined.csv"))

cb_colors <- c(
  "Oxygen"         = "#E69F00",
  "OD600"          = "#56B4E9",
  "Flow Cytometry" = "#009E73"
)

growth_long <- growth %>%
  tidyr::pivot_longer(
    cols = c(r_O2, r_OD600, r_FC),
    names_to = "Method",
    values_to = "Growth_Rate"
  ) %>%
  dplyr::mutate(
    Method = factor(
      Method,
      levels = c("r_O2", "r_OD600", "r_FC"),
      labels = c("Oxygen", "OD600", "Flow Cytometry")
    )
  )

growth_comparison <- ggplot(growth_long, aes(x = Taxon, y = Growth_Rate, fill = Method)) +
  geom_vline(
    xintercept = seq(1.5, length(unique(growth_long$Taxon)) - 0.5),
    color      = "gray80", linetype = "dashed"
  ) +
  geom_boxplot(outlier.size = 1.2, outlier.shape = 16) +
  scale_fill_manual(values = cb_colors) +
  theme_classic(base_size = 16) +
  labs(
    x    = "Taxon",
    y    = expression(paste("Growth rate (min"^{-1},")")),
    fill = "Method"
  ) +
  theme(
    axis.title      = element_text(size = 18),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 14, face = "italic"),
    axis.text.y     = element_text(size = 14),
    legend.title    = element_text(size = 16),
    legend.text     = element_text(size = 14),
    legend.position = "top",
    strip.text      = element_text(size = 16, face = "italic")
  ) +
  facet_grid(. ~ "By Taxon")

global_comparison <- ggplot(growth_long, aes(x = Method, y = Growth_Rate, fill = Method)) +
  geom_boxplot(outlier.size = 1.2, outlier.shape = 16) +
  geom_signif(
    comparisons = list(
      c("Oxygen", "OD600"),
      c("Oxygen", "Flow Cytometry"),
      c("OD600", "Flow Cytometry")
    ),
    annotations = c("ns", "ns", "ns"),
    step_increase = 0.1,
    tip_length = 0.01
  ) +
  scale_fill_manual(values = cb_colors) +
  theme_classic(base_size = 16) +
  labs(
    x = "Method",
    y = expression(paste("Growth rate (min"^{-1},")"))
  ) +
  theme(
    axis.title  = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "plain"),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  ) +
  facet_grid(. ~ "All Taxa")

combined_comparison <- growth_comparison + global_comparison +
  patchwork::plot_layout(widths = c(3, 1))

ggsave(
  filename = file.path(plots_dir, "Fig_3_growth_rate_comparison.pdf"),
  plot     = combined_comparison,
  width    = 14,
  height   = 6
)

mm_method   <- lmer(Growth_Rate ~ 1 + Method + (1 | Taxon), REML = FALSE, data = growth_long)
mm_method_1 <- lmer(Growth_Rate ~ 1 + (1 | Taxon), REML = FALSE, data = growth_long)

anova(mm_method, mm_method_1)

aic_comparison <- AIC(mm_method, mm_method_1)
print(aic_comparison)

print(summary(mm_method))
ci_method   <- confint(mm_method, method = "Wald")
coefficients <- fixef(mm_method)

mc      <- multcomp::glht(mm_method, linfct = multcomp::mcp(Method = "Tukey"),
                          test = multcomp::adjusted("bonferroni"))
summary_mc <- summary(mc)
mc_ci      <- confint(mc)

method_effects <- data.frame(
  Method = c("Oxygen", "OD600", "Flow Cytometry"),
  Estimate = c(
    coefficients[1],
    coefficients[1] + coefficients["MethodOD600"],
    coefficients[1] + coefficients["MethodFlow Cytometry"]
  ),
  CI_lower = c(
    ci_method["(Intercept)", 1],
    ci_method["(Intercept)", 1] + ci_method["MethodOD600", 1],
    ci_method["(Intercept)", 1] + ci_method["MethodFlow Cytometry", 1]
  ),
  CI_upper = c(
    ci_method["(Intercept)", 2],
    ci_method["(Intercept)", 2] + ci_method["MethodOD600", 2],
    ci_method["(Intercept)", 2] + ci_method["MethodFlow Cytometry", 2]
  )
)

pairs <- data.frame(
  y.position = max(method_effects$CI_upper, na.rm = TRUE) + c(0.0005, 0.001, 0.0015),
  xmin = c(1, 1, 2),
  xmax = c(2, 3, 3),
  p.value = summary_mc$test$pvalues
) %>%
  dplyr::mutate(stars = dplyr::case_when(
    p.value > 0.05                     ~ "ns",
    p.value <= 0.05 & p.value > 0.01  ~ "*",
    p.value <= 0.01 & p.value > 0.001 ~ "**",
    p.value <= 0.001                  ~ "***"
  ))

method_comparison <- ggplot(method_effects, aes(x = Method, y = Estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  geom_segment(data = pairs, aes(x = xmin, xend = xmax, y = y.position, yend = y.position)) +
  geom_text(data = pairs, aes(x = (xmin + xmax) / 2, y = y.position, label = stars),
            vjust = -0.5, size = 5) +
  labs(
    x = "Method", 
    y = expression("Growth rate ("*min^{-1}*")"),
    title = "Comparison of Method Effects (Mixed Model)"
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(plots_dir, "method_effects_CI.pdf"),
  plot     = method_comparison,
  width    = 8,
  height   = 6
)

write_csv(method_effects, file.path(tables_dir, "method_effects_estimates.csv"))
write_csv(pairs,          file.path(tables_dir, "method_effects_significance.csv"))

capture.output(summary(mm_OD  <- lmer(r_O2 ~ r_OD600 + (1 | Taxon), data = growth)),
               file = file.path(tables_dir, "mixed_model_OD600_summary.txt"))
capture.output(summary(mm_FC  <- lmer(r_O2 ~ r_FC    + (1 | Taxon), data = growth)),
               file = file.path(tables_dir, "mixed_model_FC_summary.txt"))

print(method_effects)
print(summary_mc)


################################################################################
# 7. LME REGRESSIONS (O2 vs OD600 / FC)
################################################################################

create_norm_data <- function(model, x_var) {
  ranef_taxon <- ranef(model)$Taxon
  fixed_coef  <- fixef(model)
  
  rand_effects <- data.frame(
    Taxon    = rownames(ranef_taxon),
    rand_int = ranef_taxon[, 1]
  )
  
  growth_norm <- growth %>%
    dplyr::left_join(rand_effects, by = "Taxon") %>%
    dplyr::mutate(r_O2_norm = r_O2 - rand_int)
  
  r2_mixed <- cor(fitted(model), growth$r_O2, use = "complete.obs")^2
  
  eq_text <- sprintf("y = %.3f + %.3fx\nR² = %.3f", 
                     fixed_coef[1], 
                     fixed_coef[2],
                     r2_mixed)
  
  list(
    data      = growth_norm,
    fixed_coef = fixed_coef,
    eq_text   = eq_text,
    x_var     = x_var
  )
}

od_norm <- create_norm_data(mm_OD, "r_OD600")
fc_norm <- create_norm_data(mm_FC, "r_FC")

rate_range_od <- range(c(od_norm$data$r_O2_norm, od_norm$data$r_OD600), na.rm = TRUE)
rate_range_fc <- range(c(fc_norm$data$r_O2_norm, fc_norm$data$r_FC),    na.rm = TRUE)
rate_range    <- range(c(rate_range_od, rate_range_fc))

okabe_ito_extended <- c(
  "#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
  "#D55E00","#CC79A7","#000000","#E69F00B0","#56B4E9B0",
  "#009E73B0","#0072B2B0","#D55E00B0","#999999","#AA4499"
)

p_od <- ggplot(od_norm$data, aes(x = r_OD600, y = r_O2_norm)) +
  geom_point(aes(color = Taxon), size = 3) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_abline(intercept = od_norm$fixed_coef[1], slope = od_norm$fixed_coef[2],
              color = "black", linetype = "solid") +
  scale_color_manual(values = okabe_ito_extended) +
  theme_classic() +
  labs(
    x = expression(paste("Growth Rate - OD"[600], " (min"^-1, ")")),
    y = expression(paste("Growth Rate - O"[2], " (min"^-1, ")"))
  ) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  ) +
  coord_cartesian(xlim = rate_range, ylim = rate_range) +
  annotate("text",
           x = min(rate_range) + 0.1 * diff(rate_range),
           y = max(rate_range) - 0.1 * diff(rate_range),
           label = od_norm$eq_text, hjust = 0, vjust = 1, size = 5) +
  ggtitle("(A)")

p_fc <- ggplot(fc_norm$data, aes(x = r_FC, y = r_O2_norm)) +
  geom_point(aes(color = Taxon), size = 3) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_abline(intercept = fc_norm$fixed_coef[1], slope = fc_norm$fixed_coef[2],
              color = "black", linetype = "solid") +
  scale_color_manual(values = okabe_ito_extended) +
  theme_classic() +
  labs(
    x = expression(paste("Growth Rate - FC (min"^-1, ")")),
    y = expression(paste("Growth Rate - O"[2], " (min"^-1, ")"))
  ) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = rate_range, ylim = rate_range) +
  annotate("text",
           x = min(rate_range) + 0.1 * diff(rate_range),
           y = max(rate_range) - 0.1 * diff(rate_range),
           label = fc_norm$eq_text, hjust = 0, vjust = 1, size = 5) +
  ggtitle("(B)")

combined_norm <- p_od / p_fc

ggsave(
  filename = file.path(plots_dir, "Fig_4_growth_rate_regression_normalized_combined.pdf"),
  plot     = combined_norm,
  width    = 10,
  height   = 16
)


################################################################################
# 8. BLAND–ALTMAN PLOTS (All replicates) WITH LME REGRESSION
################################################################################

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && is.finite(a)) a else b

df_raw <- read_csv(file.path(tables_dir, "growth_rates_combined.csv"),
                   show_col_types = FALSE) %>%
  dplyr::rename(
    Oxygen_r = r_O2,
    OD_r     = r_OD600,
    FC_r     = r_FC
  ) %>%
  dplyr::mutate(Taxon = as.factor(Taxon))

get_lmer_slope_p <- function(fit) {
  tab <- coef(summary(fit))
  if (is.null(tab)) return(NA_real_)
  rn <- rownames(tab)
  if (is.null(rn)) return(NA_real_)
  idx <- which(rn == "avg")
  if (!length(idx)) return(NA_real_)
  val <- tab[idx, "Pr(>|t|)"]
  if (length(val) == 0) return(NA_real_) else as.numeric(val)
}

bland_altman_plot_lme <- function(df, method1, method2, label1, label2, panel_letter) {
  df_ba <- df %>%
    dplyr::mutate(
      method1 = {{ method1 }},
      method2 = {{ method2 }},
      avg  = (method1 + method2) / 2,
      diff = method1 - method2
    ) %>%
    dplyr::filter(is.finite(avg), is.finite(diff))
  
  bias    <- mean(df_ba$diff, na.rm = TRUE)
  sd_diff <- sd(df_ba$diff,   na.rm = TRUE)
  upper   <- bias + 1.96 * sd_diff
  lower   <- bias - 1.96 * sd_diff
  within_limits <- mean(df_ba$diff >= lower & df_ba$diff <= upper, na.rm = TRUE) * 100
  
  df_reg <- df_ba %>% dplyr::filter(diff >= lower, diff <= upper)
  
  use_lmer <- dplyr::n_distinct(df_reg$Taxon) >= 2 && nrow(df_reg) >= 5
  
  if (use_lmer) {
    mm_fit <- lmer(diff ~ avg + (1 | Taxon), data = df_reg, REML = TRUE)
    fe <- fixef(mm_fit)
    intercept <- unname(fe[1]) %||% NA_real_
    slope     <- unname(fe[2]) %||% NA_real_
    r2_mixed  <- cor(fitted(mm_fit), df_reg$diff, use = "complete.obs")^2
    p_value   <- get_lmer_slope_p(mm_fit)
    model_used <- "lmer(diff ~ avg + (1 | Taxon))"
  } else {
    mm_fit <- lm(diff ~ avg, data = df_reg)
    coefs <- coef(mm_fit)
    intercept <- unname(coefs[1]) %||% NA_real_
    slope     <- unname(coefs[2]) %||% NA_real_
    r2_mixed  <- cor(fitted(mm_fit), df_reg$diff, use = "complete.obs")^2
    slope_row <- summary(mm_fit)$coefficients
    p_value   <- if (!is.null(slope_row) && nrow(slope_row) >= 2) slope_row["avg", "Pr(>|t|)"] else NA_real_
    model_used <- "lm(diff ~ avg)"
  }
  
  cat(glue::glue(
    "\n[{panel_letter}] Bland–Altman {label1} vs {label2}\n",
    "  Model: {model_used}\n",
    "  N (all points): {nrow(df_ba)},  N (within LoA for regression): {nrow(df_reg)}\n",
    "  Bias: {round(bias, 6)}\n",
    "  SD(diff): {round(sd_diff, 6)}\n",
    "  Upper LoA: {round(upper, 6)}\n",
    "  Lower LoA: {round(lower, 6)}\n",
    "  % within LoA: {round(within_limits, 1)}%\n",
    "  Intercept: {round(intercept, 6)}\n",
    "  Slope: {round(slope, 6)}\n",
    "  R^2 (fitted vs. observed diffs, within LoA): {round(r2_mixed, 6)}\n",
    "  p-value for slope: {ifelse(is.na(p_value), 'NA', format.pval(p_value, digits = 4))}\n"
  ))
  
  ggplot(df_ba, aes(x = avg, y = diff)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper,
             fill = "#D55E00", alpha = 0.15) +
    geom_point(size = 2.5, alpha = 0.9, color = "black") +
    geom_hline(yintercept = bias,  color = "black",    linewidth = 0.7) +
    geom_hline(yintercept = upper, color = "firebrick", linewidth = 0.8, linetype = "dashed") +
    geom_hline(yintercept = lower, color = "firebrick", linewidth = 0.8, linetype = "dashed") +
    annotate("text", x = min(df_ba$avg, na.rm = TRUE), y = upper,
             label = sprintf("Upper LoA = %.3f", upper),
             hjust = 0, vjust = -0.8, size = 4, color = "black", fontface = "bold") +
    annotate("text", x = min(df_ba$avg, na.rm = TRUE), y = lower,
             label = sprintf("Lower LoA = %.3f", lower),
             hjust = 0, vjust = 1.8, size = 4, color = "black", fontface = "bold") +
    annotate("text", x = min(df_ba$avg, na.rm = TRUE), y = bias,
             label = sprintf("Bias = %.3f", bias),
             hjust = 0, vjust = -2.6, size = 4, color = "black", fontface = "bold") +
    annotate("text", x = min(df_ba$avg, na.rm = TRUE), y = lower - 0.005,
             label = sprintf("%% within limits: %.1f%%", within_limits),
             hjust = 0, size = 4, color = "black", fontface = "italic") +
    labs(
      title = glue::glue("({panel_letter}) {label1} vs {label2} — Bland–Altman (all replicates)"),
      x = "Mean growth rate (per replicate)",
      y = "Difference in growth rate (per replicate)"
    ) +
    theme_classic(base_size = 14)
}

plot1 <- bland_altman_plot_lme(df_raw, Oxygen_r, OD_r, "Oxygen", "OD600", "A")
plot2 <- bland_altman_plot_lme(df_raw, Oxygen_r, FC_r, "Oxygen", "Flow Cytometry", "B")

combined_plot <- gridExtra::grid.arrange(plot1, plot2, ncol = 2)

ggsave(
  filename = file.path(plots_dir, "Fig_5_BlandAltman_AllReplicates_LME.pdf"),
  plot     = combined_plot,
  width    = 14,
  height   = 6
)

################################################################################
# END OF SCRIPT
################################################################################
