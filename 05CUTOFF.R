################################################################################
# README – O₂ trimming + O₂-based growth–respiration + O₂ ≥ 0.5 sensitivity
#
# This script does three things:
#   1) Trims raw dissolved-oxygen (DO) time series:
#        - Automatically detects the post-peak exponential-decline window
#        - Removes very flat / noisy series
#        - Saves the trimmed trajectories and a log of skipped series
#
#   2) Fits the normalised O₂ growth–respiration model with N₀ fixed to 1:
#        O2_norm(t) = O2_0 + (resp_tot / r) * (1 - exp(r * t))
#      for each Taxon × Replicate, on the full trimmed window
#      → returns r (min⁻¹), resp_tot, O₂₀, O₂_ref and QC metrics.
#
#   3) Reviewer sensitivity test (O₂_norm ≥ 0.5 only):
#        - Re-fits the same model using only points with O₂_norm ≥ 0.5
#        - Compares r and resp_tot between the two fits
#        - Summarises ratios and differences (r₀.₅ / r, R₀.₅ / R)
#
# Outputs:
#   Tables:
#     - Tables/Oxygen_Data_Trimmed.csv
#     - Tables/Oxygen_Data_Filtered.csv
#     - Tables/Skipped_Series_Log.csv
#     - Tables/oxygen_model_results_main.csv
#     - Tables/oxygen_model_results_O2_ge_0.5.csv
#     - Tables/Table_S3_oxygen_model_comparison_O2_ge_0.5.csv
#     - Tables/oxygen_model_comparison_O2_ge_0.5_summary.csv
#
#   Plots:
#     - plots/oxygen_trimming_diagnostics.pdf
#     - plots/oxygen_dynamics_all_models.pdf
#     - plots/oxygen_dynamics_fullsize_per_page.pdf
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(zoo)
  library(minpack.lm)
  library(patchwork)
})

################################################################################
# 0. PROJECT PATHS
################################################################################

data_dir   <- "data"
plots_dir  <- "plots"
tables_dir <- "Tables"

dir.create(data_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

INPUT_CSV    <- file.path(data_dir,   "Oxygen_Data_Long.csv")
TRIMMED_CSV  <- file.path(tables_dir, "Oxygen_Data_Trimmed.csv")
FILTERED_CSV <- file.path(tables_dir, "Oxygen_Data_Filtered.csv")
SKIPPED_CSV  <- file.path(tables_dir, "Skipped_Series_Log.csv")
DIAG_PDF     <- file.path(plots_dir,  "oxygen_trimming_diagnostics.pdf")

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
readr::write_csv(trimmed, TRIMMED_CSV)

filtered <- trimmed %>%
  dplyr::select(Taxon, Replicate, Time, Oxygen)
readr::write_csv(filtered, FILTERED_CSV)

readr::write_csv(skipped_log, SKIPPED_CSV)

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
  
  p <- ggplot2::ggplot(raw_df, ggplot2::aes(Time, Oxygen)) +
    ggplot2::geom_rect(
      data = rect_df, inherit.aes = FALSE,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "orange", alpha = 0.08
    ) +
    ggplot2::geom_line(colour = "grey60", linewidth = 0.8) +
    ggplot2::labs(
      title = sid,
      subtitle = "Trimmed zone (orange)",
      x = "Time (min)", y = "O₂ (mg L⁻¹)"
    ) +
    ggplot2::theme_classic(base_size = 11)
  
  print(p)
}

dev.off()

################################################################################
# 2. O2-BASED GROWTH & TOTAL RESPIRATION (NEW MODEL, N0 FIXED TO 1)
#    + PARALLEL FIT USING ONLY OXYGEN_NORM ≥ 0.5
################################################################################

## ────────────── QC / Fit thresholds ───────────────────────────────── ##
REL_SE_THRESHOLD <- 0.15
PVAL_THRESHOLD   <- 0.001
R2_THRESHOLD     <- 0.90
MAX_RESID_RANGE  <- 0.20
AIC_IMPROVEMENT  <- 10
MAPE_MAX         <- 0.15

O2_THRESHOLD     <- 0.5   # reviewer: use only Oxygen_norm ≥ 0.5

## ────────────── Load trimmed oxygen data ───────────────────────────── ##
oxygen_data <- read_csv(FILTERED_CSV, show_col_types = FALSE) %>%
  dplyr::mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate)
  )

## ────────────── Model function (N0 fixed to 1) ─────────────────────── ##
# Oxygen_norm(t) = O2_0 + (resp_tot / r) * (1 - exp(r * t))
resp_model <- function(r, resp_tot, t, O2_0) {
  O2_0 + (resp_tot / r) * (1 - exp(r * t))
}

## ────────────── Plot theme ─────────────────────────────────────────── ##
isme_theme <- function() {
  ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 16),
      axis.text  = ggplot2::element_text(size = 14),
      axis.line  = ggplot2::element_line(linewidth = .8),
      axis.ticks = ggplot2::element_line(linewidth = .6),
      legend.position = "none"
    )
}

## ────────────── Output objects ─────────────────────────────────────── ##
results_main <- tibble(
  Taxon               = character(),
  Replicate           = character(),
  r_per_minute        = numeric(),
  r_per_hour          = numeric(),
  resp_tot            = numeric(),
  O2_0                = numeric(),
  O2_ref              = numeric(),
  AICc                = numeric(),
  lnO2_change_per_min = numeric(),
  pseudo_R2           = numeric(),
  fit_ok              = logical()
)

results_O2_ge_0.5 <- tibble(
  Taxon               = character(),
  Replicate           = character(),
  r_per_minute        = numeric(),
  r_per_hour          = numeric(),
  resp_tot            = numeric(),
  O2_0                = numeric(),
  O2_ref              = numeric(),
  AICc                = numeric(),
  lnO2_change_per_min = numeric(),
  pseudo_R2           = numeric(),
  fit_ok              = logical()
)

plots_list <- list()

## ────────────── Fitting loop ───────────────────────────────────────── ##
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
  
  if (nrow(df) < 5) next
  
  # Normalise O2 for fitting
  O0 <- mean(head(df$Oxygen, 3), na.rm = TRUE)   # O2_ref (mg L^-1)
  df <- df %>% dplyr::mutate(Oxygen_norm = Oxygen / O0)
  
  # starting values for main fit (full trimmed series)
  r_start <- {
    seg    <- head(df, max(3, floor(.3 * nrow(df))))
    slopes <- abs(diff(log(pmax(seg$Oxygen_norm, 1e-6))) / diff(seg$Time0))
    pmin(pmax(max(slopes, na.rm = TRUE), 1e-4), 5e-2)
  }
  resp_tot_start <- {
    slope <- suppressWarnings(min(diff(df$Oxygen_norm) / diff(df$Time0), na.rm = TRUE))
    if (!is.finite(slope)) slope <- -1e-3
    k_guess <- abs(slope)            # because dO2_norm/dt|0 ≈ -resp_tot
    k_guess <- pmin(pmax(k_guess, 1e-5), 0.1)
    k_guess
  }
  
  # ── MAIN FIT: full trimmed series ────────────────────────────────────
  fit <- tryCatch(
    minpack.lm::nlsLM(
      Oxygen_norm ~ resp_model(r, resp_tot, Time0, O2_0),
      data    = df,
      start   = list(r = r_start, resp_tot = resp_tot_start, O2_0 = 1),
      lower   = c(r = 1e-4,  resp_tot = 1e-5, O2_0 = .8),
      upper   = c(r = .1,    resp_tot = 0.5,  O2_0 = 1.2),
      control = minpack.lm::nls.lm.control(maxiter = 300)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) next
  
  # Outlier removal and refit
  pred    <- stats::predict(fit, df)
  keep_ix <- abs(df$Oxygen_norm - pred) <
    2 * stats::sd(df$Oxygen_norm - pred, na.rm = TRUE)
  df_kept <- df[keep_ix, ]
  if (nrow(df_kept) < 5) next
  
  fit <- tryCatch(stats::update(fit, data = df_kept), error = function(e) fit)
  
  pars  <- summary(fit)$coefficients
  r_est <- pars["r", "Estimate"]
  K_est <- pars["resp_tot", "Estimate"]
  
  pseudo_R2 <- 1 - sum(stats::residuals(fit)^2) /
    sum((df_kept$Oxygen_norm - mean(df_kept$Oxygen_norm))^2)
  
  fit_ok_main <- all(
    abs(pars[, "Std. Error"] / pars[, "Estimate"]) < REL_SE_THRESHOLD,
    pars[, "Pr(>|t|)"] < PVAL_THRESHOLD,
    pseudo_R2 >= R2_THRESHOLD,
    diff(range(stats::residuals(fit), na.rm = TRUE)) < MAX_RESID_RANGE,
    is.finite(K_est),
    K_est > 0, K_est < 0.5,
    dplyr::between(r_est, 1e-4, 0.1),
    AIC(stats::lm(Oxygen_norm ~ 1, data = df_kept)) - AIC(fit) >= AIC_IMPROVEMENT,
    mean(abs(stats::residuals(fit) / df_kept$Oxygen_norm)) < MAPE_MAX
  )
  
  # store main fit
  results_main <- results_main %>%
    dplyr::add_row(
      Taxon               = Tax,
      Replicate           = Rep,
      r_per_minute        = r_est,
      r_per_hour          = r_est * 60,
      resp_tot            = K_est,
      O2_0                = coef(fit)["O2_0"],
      O2_ref              = O0,
      AICc                = AIC(fit),
      lnO2_change_per_min = (log(dplyr::last(df_kept$Oxygen_norm)) -
                               log(dplyr::first(df_kept$Oxygen_norm))) /
        (dplyr::last(df_kept$Time0) - dplyr::first(df_kept$Time0)),
      pseudo_R2           = pseudo_R2,
      fit_ok              = fit_ok_main
    )
  
  # store plot
  df_kept <- df_kept %>%
    dplyr::mutate(Pred = stats::predict(fit, df_kept))
  
  plot_key <- paste(Tax, Rep, sep = "_")
  plots_list[[plot_key]] <-
    ggplot2::ggplot(df_kept, ggplot2::aes(Time0, Oxygen_norm)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_line(ggplot2::aes(y = Pred), linewidth = .8, colour = "red") +
    ggplot2::labs(
      title = plot_key,
      x     = expression(Time~"(min)"),
      y     = expression(Normalised~O[2])
    ) +
    isme_theme()
  
  # ── REVIEWER TEST: FIT USING ONLY OXYGEN_NORM ≥ 0.5 ──────────────────
  df_05 <- df %>% dplyr::filter(Oxygen_norm >= O2_THRESHOLD)
  if (nrow(df_05) >= 5) {
    
    r_start_05 <- {
      seg    <- head(df_05, max(3, floor(.3 * nrow(df_05))))
      slopes <- abs(diff(log(pmax(seg$Oxygen_norm, 1e-6))) / diff(seg$Time0))
      pmin(pmax(max(slopes, na.rm = TRUE), 1e-4), 5e-2)
    }
    resp_tot_start_05 <- {
      slope <- suppressWarnings(min(diff(df_05$Oxygen_norm) / diff(df_05$Time0),
                                    na.rm = TRUE))
      if (!is.finite(slope)) slope <- -1e-3
      k_guess <- abs(slope)
      k_guess <- pmin(pmax(k_guess, 1e-5), 0.1)
      k_guess
    }
    
    fit_05 <- tryCatch(
      minpack.lm::nlsLM(
        Oxygen_norm ~ resp_model(r, resp_tot, Time0, O2_0),
        data    = df_05,
        start   = list(r = r_start_05, resp_tot = resp_tot_start_05, O2_0 = 1),
        lower   = c(r = 1e-4,  resp_tot = 1e-5, O2_0 = .8),
        upper   = c(r = .1,    resp_tot = 0.5,  O2_0 = 1.2),
        control = minpack.lm::nls.lm.control(maxiter = 300)
      ),
      error = function(e) NULL
    )
    
    if (!is.null(fit_05)) {
      pred_05    <- stats::predict(fit_05, df_05)
      keep_ix_05 <- abs(df_05$Oxygen_norm - pred_05) <
        2 * stats::sd(df_05$Oxygen_norm - pred_05, na.rm = TRUE)
      df_05_kept <- df_05[keep_ix_05, ]
      
      if (nrow(df_05_kept) >= 5) {
        fit_05 <- tryCatch(stats::update(fit_05, data = df_05_kept),
                           error = function(e) fit_05)
        
        pars_05  <- summary(fit_05)$coefficients
        r_est_05 <- pars_05["r", "Estimate"]
        K_est_05 <- pars_05["resp_tot", "Estimate"]
        
        pseudo_R2_05 <- 1 - sum(stats::residuals(fit_05)^2) /
          sum((df_05_kept$Oxygen_norm -
                 mean(df_05_kept$Oxygen_norm))^2)
        
        fit_ok_05 <- all(
          abs(pars_05[, "Std. Error"] / pars_05[, "Estimate"]) < REL_SE_THRESHOLD,
          pars_05[, "Pr(>|t|)"] < PVAL_THRESHOLD,
          pseudo_R2_05 >= R2_THRESHOLD,
          diff(range(stats::residuals(fit_05), na.rm = TRUE)) < MAX_RESID_RANGE,
          is.finite(K_est_05),
          K_est_05 > 0, K_est_05 < 0.5,
          dplyr::between(r_est_05, 1e-4, 0.1),
          AIC(stats::lm(Oxygen_norm ~ 1, data = df_05_kept)) - AIC(fit_05) >= AIC_IMPROVEMENT,
          mean(abs(stats::residuals(fit_05) / df_05_kept$Oxygen_norm)) < MAPE_MAX
        )
        
        results_O2_ge_0.5 <- results_O2_ge_0.5 %>%
          dplyr::add_row(
            Taxon               = Tax,
            Replicate           = Rep,
            r_per_minute        = r_est_05,
            r_per_hour          = r_est_05 * 60,
            resp_tot            = K_est_05,
            O2_0                = coef(fit_05)["O2_0"],
            O2_ref              = O0,
            AICc                = AIC(fit_05),
            lnO2_change_per_min = (log(dplyr::last(df_05_kept$Oxygen_norm)) -
                                     log(dplyr::first(df_05_kept$Oxygen_norm))) /
              (dplyr::last(df_05_kept$Time0) - dplyr::first(df_05_kept$Time0)),
            pseudo_R2           = pseudo_R2_05,
            fit_ok              = fit_ok_05
          )
      }
    }
  }
}

################################################################################
# 3. SAVE MODEL RESULTS
################################################################################

readr::write_csv(results_main,
                 file.path(tables_dir, "oxygen_model_results_main.csv"))
readr::write_csv(results_O2_ge_0.5,
                 file.path(tables_dir, "oxygen_model_results_O2_ge_0.5.csv"))

# save diagnostic curves
if (length(plots_list) > 0) {
  pdf(file.path(plots_dir, "oxygen_dynamics_all_models.pdf"),
      width = 14, height = 10)
  print(patchwork::wrap_plots(plots_list))
  dev.off()
  
  pdf(file.path(plots_dir, "oxygen_dynamics_fullsize_per_page.pdf"),
      width = 8, height = 6)
  purrr::walk(plots_list, print)
  dev.off()
}

################################################################################
# 4. COMPARISON: GROWTH RATE & RESPIRATION
#    (FULL TRIMMED) vs (OXYGEN_NORM ≥ 0.5)
################################################################################

comparison_05 <- results_main %>%
  dplyr::filter(fit_ok) %>%
  dplyr::select(
    Taxon,
    Replicate,
    r_main = r_per_minute,
    R_main = resp_tot
  ) %>%
  dplyr::inner_join(
    results_O2_ge_0.5 %>%
      dplyr::filter(fit_ok) %>%
      dplyr::select(
        Taxon,
        Replicate,
        r_05 = r_per_minute,
        R_05 = resp_tot
      ),
    by = c("Taxon", "Replicate")
  ) %>%
  dplyr::mutate(
    r_ratio = r_05 / r_main,
    r_diff  = r_05 - r_main,
    R_ratio = R_05 / R_main,
    R_diff  = R_05 - R_main
  )

readr::write_csv(
  comparison_05,
  file.path(tables_dir, "Table_S3_oxygen_model_comparison_O2_ge_0.5.csv")
)

cat("\n==== Comparison: full trimmed fit vs O2_norm ≥ 0.5 fit ====\n")
if (nrow(comparison_05) > 0) {
  cat("N paired fits: ", nrow(comparison_05), "\n")
  
  ## --- Growth rate r ---
  cat("Correlation r_main vs r_05: ",
      cor(comparison_05$r_main, comparison_05$r_05, use = "complete.obs"), "\n")
  cat("Median r_ratio (r_05 / r_main): ",
      median(comparison_05$r_ratio, na.rm = TRUE), "\n")
  cat("IQR r_ratio: ",
      paste(signif(quantile(comparison_05$r_ratio,
                            probs = c(0.25, 0.75),
                            na.rm = TRUE), 3),
            collapse = " – "),
      "\n")
  cat("Max |r_ratio - 1|: ",
      max(abs(comparison_05$r_ratio - 1), na.rm = TRUE), "\n")
  
  ## --- Respiration R (resp_tot) ---
  cat("Correlation R_main vs R_05: ",
      cor(comparison_05$R_main, comparison_05$R_05, use = "complete.obs"), "\n")
  cat("Median R_ratio (R_05 / R_main): ",
      median(comparison_05$R_ratio, na.rm = TRUE), "\n")
  cat("IQR R_ratio: ",
      paste(signif(quantile(comparison_05$R_ratio,
                            probs = c(0.25, 0.75),
                            na.rm = TRUE), 3),
            collapse = " – "),
      "\n")
  cat("Max |R_ratio - 1|: ",
      max(abs(comparison_05$R_ratio - 1), na.rm = TRUE), "\n")
  
} else {
  cat("No overlapping good fits between main and O2 ≥ 0.5 analyses.\n")
}

# ───────────────── 5. Save global summary of r_main vs r_05 and R_main vs R_05 ───────────── #

if (nrow(comparison_05) > 0) {
  comparison_05_summary <- comparison_05 %>%
    dplyr::summarise(
      n_paired_fits             = dplyr::n(),
      
      # Growth rate r
      cor_r_main_r_05           = cor(r_main, r_05, use = "complete.obs"),
      median_r_ratio            = median(r_ratio, na.rm = TRUE),
      q25_r_ratio               = quantile(r_ratio, 0.25, na.rm = TRUE),
      q75_r_ratio               = quantile(r_ratio, 0.75, na.rm = TRUE),
      max_abs_r_ratio_minus1    = max(abs(r_ratio - 1), na.rm = TRUE),
      
      # Respiration R (resp_tot)
      cor_R_main_R_05           = cor(R_main, R_05, use = "complete.obs"),
      median_R_ratio            = median(R_ratio, na.rm = TRUE),
      q25_R_ratio               = quantile(R_ratio, 0.25, na.rm = TRUE),
      q75_R_ratio               = quantile(R_ratio, 0.75, na.rm = TRUE),
      max_abs_R_ratio_minus1    = max(abs(R_ratio - 1), na.rm = TRUE)
    )
  
  readr::write_csv(
    comparison_05_summary,
    file.path(tables_dir, "oxygen_model_comparison_O2_ge_0.5_summary.csv")
  )
}

################################################################################
# END OF SCRIPT
################################################################################
