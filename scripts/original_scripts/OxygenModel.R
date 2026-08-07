################################################################################
# FULL PIPELINE (ALL TAXA) — CORE KEPT, FIGURE ORDER UPDATED, TIFF OUTPUTS
#
# Requested changes:
# 1) RIS plot becomes Fig 6 and is moved to the end (after Fig 5).
# 2) Figures are saved as TIFF (instead of PDF). (Multi-page diagnostics remain PDF.)
# 3) Fig 6 legend uses 3 rows.
# 4) Fig 5 Bland–Altman long title is wrapped so it fits.
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(zoo)
  library(minpack.lm)
  library(lme4)
  library(lmerTest)
  library(grid)
  library(viridis)
  
  
  # appended analyses
  library(patchwork)
  library(ggsignif)
  library(multcomp)
  library(gridExtra)
  library(glue)
  library(stringr)
})



################################################################################
# CONSTANTS for carbon-unit fluxes
################################################################################
O2_to_C_mass  <- 12/32   # mg C per mg O2 (1 mol O2 : 1 mol C respired)
C_per_cell_fg <- 100    # fg C per cell (constant; your chosen value)



################################################################################
# 0) PATHS
################################################################################

data_dir   <- "data"
plots_dir  <- "plots"
tables_dir <- "Tables"

dir.create(data_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

OXYGEN_CSV <- file.path(data_dir, "Oxygen_Data_Long.csv")

# <-- YOUR FILE
NINOC_TABLE_CSV <- file.path(data_dir, "Ninoc_and_deltaTime_to_N0.csv")

TRIMMED_CSV  <- file.path(tables_dir, "Oxygen_Data_Trimmed.csv")
FILTERED_CSV <- file.path(tables_dir, "Oxygen_Data_Filtered.csv")
SKIPPED_CSV  <- file.path(tables_dir, "Skipped_Series_Log.csv")
DIAG_PDF     <- file.path(plots_dir,  "oxygen_trimming_diagnostics.pdf")

RESULTS_FIT_CSV    <- file.path(tables_dir, "oxygen_fit_results.csv")
RESULTS_FINAL_CSV  <- file.path(tables_dir, "oxygen_results_with_R.csv")

# fitted curves + per-series fits (multi-page PDFs kept as PDF)
FITCURVES_CSV <- file.path(tables_dir, "oxygen_fit_curves.csv")
FITS_PDF      <- file.path(plots_dir,  "oxygen_model_fit_curves.pdf")

# ---- FIGURE OUTPUTS (TIFF) ----
FIG2_FACET4_TIF     <- file.path(plots_dir, "Fig_2_oxygen_dynamics_facet4_no_overflow.tiff")
FIG3_GROWTH_TIF     <- file.path(plots_dir, "Fig_3_growth_rate_comparison.tiff")
FIG4_REGRESS_TIF    <- file.path(plots_dir, "Fig_4_growth_rate_regression_normalized_combined.tiff")
FIG5_BA_TIF         <- file.path(plots_dir, "Fig_5_BlandAltman_AllReplicates_LME.tiff")
FIG6_RIS_MAIN_TIF   <- file.path(plots_dir, "Fig_6_r_vs_R_RIS_derand_full_MAIN.tiff")

SUPP_S1_RESID_TIF   <- file.path(plots_dir, "Supp_Fig_S1_residuals.tiff")
SUPP_FIG3_NORM_TIF  <- file.path(plots_dir, "supp_Fig_3_oxygen_all_replicates_NORMALISED_by_replicate.tiff")
METHOD_EFFECTS_TIF  <- file.path(plots_dir, "method_effects_CI.tiff")

################################################################################
# ============================ CORE (ADAPTED) ==================================
# (Core computations are the same; only RIS figure moved to Fig 6 and TIFFs used)
################################################################################

################################################################################
# 1) LOAD RAW OXYGEN
################################################################################

raw <- readr::read_csv(OXYGEN_CSV, show_col_types = FALSE) %>%
  dplyr::mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate)
  ) %>%
  dplyr::arrange(Taxon, Replicate, Time)

stopifnot(all(c("Taxon","Replicate","Time","Oxygen") %in% names(raw)))

raw <- raw %>%
  dplyr::mutate(series_id = paste(Taxon, Replicate, sep = " | "))

################################################################################
# 2) TRIM OXYGEN TIME SERIES (peak → 2nd inflection)
################################################################################

SPLINE_SPAR   <- 0.40
RUN_LEN       <- 10
REL_PROP      <- 0.008
ABS_THR       <- 0.0003
WINDOW_LEN    <- 4
LEVEL_DELTA   <- 0.0025
FLAT_RANGE_OK <- 0.05

find_plateau <- function(o2_vec, peak_idx,
                         run_len = 2, rel_prop = 0.01, abs_thr = 0.0003,
                         window_len = 4, level_delta = 0.0025) {
  n <- length(o2_vec)
  if (peak_idx >= n) return(n)
  
  slopes   <- diff(o2_vec)
  max_down <- abs(min(slopes[peak_idx:(n - 1)], na.rm = TRUE))
  slope_thr <- max(abs_thr, rel_prop * max_down)
  
  rng <- rep(NA_real_, n)
  if (n >= window_len) {
    rng[1:(n - window_len + 1)] <-
      zoo::rollapply(o2_vec, window_len,
                     FUN = function(x) max(x) - min(x),
                     align = "left", fill = NA_real_)
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
  slopes <- diff(o2_fit)
  curves <- diff(slopes)
  post_pk <- seq(idx_peak + 1, length(curves))
  if (length(post_pk) == 0) return(NA_integer_)
  
  steepest <- which.min(slopes[post_pk]) + idx_peak
  win <- seq(steepest + 1, idx_plate - 1)
  win <- win[win >= 1 & win <= length(curves)]
  win <- win[!is.na(curves[win])]
  if (length(win) < 3) return(NA_integer_)
  
  win[which.max(curves[win])]
}

series_ids  <- unique(raw$series_id)
trimmed_lst <- vector("list", length(series_ids))
names(trimmed_lst) <- series_ids
skipped_log <- tibble::tibble(series_id = character(), reason = character())

for (sid in series_ids) {
  
  df <- raw %>%
    dplyr::filter(series_id == sid) %>%
    dplyr::arrange(Time)
  
  if (nrow(df) < 5) {
    skipped_log <- dplyr::add_row(skipped_log, series_id = sid, reason = "Too few rows")
    next
  }
  
  ss <- tryCatch(stats::smooth.spline(df$Time, df$Oxygen, spar = SPLINE_SPAR),
                 error = function(e) NULL)
  if (is.null(ss)) {
    skipped_log <- dplyr::add_row(skipped_log, series_id = sid, reason = "Spline failed")
    next
  }
  
  df$O2_fit <- stats::predict(ss, x = df$Time)$y
  
  if ((max(df$O2_fit, na.rm = TRUE) - min(df$O2_fit, na.rm = TRUE)) < FLAT_RANGE_OK) {
    skipped_log <- dplyr::add_row(skipped_log, series_id = sid, reason = "Too flat")
    next
  }
  
  idx_peak  <- which.max(df$Oxygen)
  idx_plate <- min(
    find_plateau(df$O2_fit, idx_peak,
                 run_len = RUN_LEN, rel_prop = REL_PROP, abs_thr = ABS_THR,
                 window_len = WINDOW_LEN, level_delta = LEVEL_DELTA),
    nrow(df)
  )
  idx_second <- find_second_inflection(df$O2_fit, idx_peak, idx_plate)
  
  if (is.na(idx_second) || idx_second <= idx_peak || idx_second > nrow(df)) {
    idx_end <- nrow(df)
    used_full <- TRUE
  } else {
    idx_end <- idx_second
    used_full <- FALSE
  }
  
  trimmed_lst[[sid]] <- df[idx_peak:idx_end, ] %>%
    dplyr::mutate(used_full_series = used_full)
}

trimmed <- dplyr::bind_rows(trimmed_lst)
readr::write_csv(trimmed, TRIMMED_CSV)
readr::write_csv(skipped_log, SKIPPED_CSV)

filtered <- trimmed %>% dplyr::select(Taxon, Replicate, Time, Oxygen)
readr::write_csv(filtered, FILTERED_CSV)

# Diagnostics PDF (multi-page; kept as PDF)
grDevices::pdf(DIAG_PDF, width = 7, height = 5, useDingbats = FALSE)
for (sid in names(trimmed_lst)) {
  df_trim <- trimmed_lst[[sid]]
  if (is.null(df_trim) || nrow(df_trim) == 0) next
  
  raw_df <- raw %>% dplyr::filter(series_id == sid) %>% dplyr::arrange(Time)
  xmin_val <- min(df_trim$Time, na.rm = TRUE)
  xmax_val <- max(df_trim$Time, na.rm = TRUE)
  
  p <- ggplot2::ggplot(raw_df, ggplot2::aes(Time, Oxygen)) +
    ggplot2::geom_rect(
      data = tibble::tibble(xmin = xmin_val, xmax = xmax_val, ymin = -Inf, ymax = Inf),
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "orange", alpha = 0.10
    ) +
    ggplot2::geom_line(colour = "grey40", linewidth = 0.8) +
    ggplot2::labs(title = sid, x = "Time (min)", y = expression(O[2]~"(mg L"^{-1}*")")) +
    ggplot2::theme_classic(base_size = 11)
  
  print(p)
}
grDevices::dev.off()

################################################################################
# 3) FIT O2 MODEL ON TRIMMED DATA (per Taxon × Replicate)
#    + EXPORT FITTED CURVES (CSV) + PER-SERIES FIT PDF
################################################################################

oxygen_data <- readr::read_csv(FILTERED_CSV, show_col_types = FALSE) %>%
  dplyr::mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate)
  ) %>%
  dplyr::arrange(Taxon, Replicate, Time)

resp_model <- function(r, K, t, O2_0) {
  O2_0 + (K / r) * (1 - exp(r * t))
}

REL_SE_THRESHOLD <- 0.15
PVAL_THRESHOLD   <- 0.001
R2_THRESHOLD     <- 0.90
MAX_RESID_RANGE  <- 0.20
AIC_IMPROVEMENT  <- 10
MAPE_MAX         <- 0.15

results_fit <- tibble::tibble(
  Taxon           = character(),
  Replicate       = character(),
  r_per_minute    = numeric(),
  K               = numeric(),
  O2_0            = numeric(),
  O2_ref          = numeric(),
  T_end_min       = numeric(),
  C_tot_mg_per_L  = numeric(),
  pseudo_R2       = numeric(),
  AIC             = numeric(),
  fit_ok          = logical()
)

fit_curves_lst <- list()

combos <- oxygen_data %>%
  dplyr::group_by(Taxon, Replicate) %>%
  dplyr::group_keys()

# per-series fit plots PDF (multi-page; kept as PDF)
grDevices::pdf(FITS_PDF, width = 7.2, height = 5.2, useDingbats = FALSE)
on.exit(grDevices::dev.off(), add = TRUE)

for (i in seq_len(nrow(combos))) {
  
  Tax <- combos$Taxon[i]
  Rep <- combos$Replicate[i]
  
  df0 <- oxygen_data %>%
    dplyr::filter(Taxon == Tax, Replicate == Rep) %>%
    dplyr::arrange(Time)
  
  if (nrow(df0) < 5) next
  
  # onset trim via derivative
  df0 <- df0 %>%
    dplyr::mutate(
      dO2   = c(NA, diff(Oxygen)),
      dt    = c(NA, diff(Time)),
      dO2dt = dO2 / dt,
      sm    = zoo::rollmean(dO2dt, 3, fill = NA, align = "right")
    )
  
  idx <- which(df0$sm < -1e-7)[1]
  idx <- if (is.na(idx) || idx > nrow(df0) - 2) 1 else min(idx + 5, nrow(df0) - 2)
  
  df <- df0[idx:nrow(df0), ] %>%
    dplyr::mutate(Time0 = Time - min(Time, na.rm = TRUE))
  
  T_end <- max(df$Time0, na.rm = TRUE)
  if (!is.finite(T_end) || T_end <= 0) next
  
  O0 <- mean(head(df$Oxygen, 3), na.rm = TRUE)
  if (!is.finite(O0) || O0 <= 0) next
  
  df <- df %>% dplyr::mutate(Oxygen_norm = Oxygen / O0)
  
  r_start <- {
    seg <- head(df, max(3, floor(0.3 * nrow(df))))
    slopes <- abs(diff(log(pmax(seg$Oxygen_norm, 1e-6))) / diff(seg$Time0))
    pmin(pmax(max(slopes, na.rm = TRUE), 1e-4), 5e-2)
  }
  
  K_start <- {
    slope <- suppressWarnings(min(diff(df$Oxygen_norm) / diff(df$Time0), na.rm = TRUE))
    if (!is.finite(slope)) slope <- -1e-3
    k_guess <- abs(slope)
    pmin(pmax(k_guess, 1e-5), 0.1)
  }
  
  fit <- tryCatch(
    minpack.lm::nlsLM(
      Oxygen_norm ~ resp_model(r, K, Time0, O2_0),
      data    = df,
      start   = list(r = r_start, K = K_start, O2_0 = 1),
      lower   = c(r = 1e-4,  K = 1e-5, O2_0 = 0.8),
      upper   = c(r = 0.1,   K = 0.5,  O2_0 = 1.2),
      control = minpack.lm::nls.lm.control(maxiter = 300)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) next
  
  pars <- summary(fit)$coefficients
  r_est    <- pars["r", "Estimate"]
  K_est    <- pars["K", "Estimate"]
  O2_0_est <- pars["O2_0", "Estimate"]
  
  pred  <- stats::predict(fit, df)
  resid <- df$Oxygen_norm - pred
  
  pseudo_R2 <- 1 - sum(resid^2, na.rm = TRUE) /
    sum((df$Oxygen_norm - mean(df$Oxygen_norm, na.rm = TRUE))^2, na.rm = TRUE)
  
  fit_ok <- all(
    abs(pars[, "Std. Error"] / pars[, "Estimate"]) < REL_SE_THRESHOLD,
    pars[, "Pr(>|t|)"] < PVAL_THRESHOLD,
    pseudo_R2 >= R2_THRESHOLD,
    diff(range(resid, na.rm = TRUE)) < MAX_RESID_RANGE,
    is.finite(K_est), K_est > 0, K_est < 0.5,
    is.finite(r_est), dplyr::between(r_est, 1e-4, 0.1),
    stats::AIC(stats::lm(Oxygen_norm ~ 1, data = df)) - stats::AIC(fit) >= AIC_IMPROVEMENT,
    mean(abs(resid / df$Oxygen_norm), na.rm = TRUE) < MAPE_MAX
  )
  
  C_tot <- (K_est / r_est) * (exp(r_est * T_end) - 1) * O0
  
  results_fit <- results_fit %>%
    dplyr::add_row(
      Taxon          = Tax,
      Replicate      = Rep,
      r_per_minute   = r_est,
      K              = K_est,
      O2_0           = O2_0_est,
      O2_ref         = O0,
      T_end_min      = T_end,
      C_tot_mg_per_L = C_tot,
      pseudo_R2      = pseudo_R2,
      AIC            = stats::AIC(fit),
      fit_ok         = fit_ok
    )
  
  curve_df <- tibble::tibble(
    Taxon        = Tax,
    Replicate    = Rep,
    Time_min     = df$Time,
    Time0_min    = df$Time0,
    Oxygen_raw   = df$Oxygen,
    Oxygen_norm  = df$Oxygen_norm,
    Fit_norm     = pred,
    Fit_raw      = pred * O0,
    r_per_minute = r_est,
    K            = K_est,
    O2_0         = O2_0_est,
    O2_ref       = O0,
    pseudo_R2    = pseudo_R2,
    AIC          = stats::AIC(fit),
    fit_ok       = fit_ok
  )
  fit_curves_lst[[paste(Tax, Rep, sep = " | ")]] <- curve_df
  
  df_plot_all <- df0 %>% dplyr::mutate(Time0 = Time - min(df$Time, na.rm = TRUE))
  
  p_fit <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = df_plot_all,
      ggplot2::aes(x = Time0, y = Oxygen),
      colour = "grey70", linewidth = 0.9
    ) +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = Time0, y = Oxygen),
      size = 1.9, alpha = 0.9
    ) +
    ggplot2::geom_line(
      data = curve_df,
      ggplot2::aes(x = Time0_min, y = Fit_raw),
      colour = "black", linewidth = 1.1
    ) +
    ggplot2::labs(
      title = paste0(Tax, " | ", Rep),
      subtitle = sprintf("fit_ok=%s   r=%.4g min^-1   K=%.4g   R2=%.3f",
                         fit_ok, r_est, K_est, pseudo_R2),
      x = "Time since fit start (min)",
      y = expression(O[2]~"(mg L"^{-1}*")")
    ) +
    ggplot2::theme_classic(base_size = 11)
  
  print(p_fit)
}

try(grDevices::dev.off(), silent = TRUE)

readr::write_csv(results_fit, RESULTS_FIT_CSV)

fit_curves <- dplyr::bind_rows(fit_curves_lst)
readr::write_csv(fit_curves, FITCURVES_CSV)

################################################################################
# 4) READ Ninoc_and_deltaTime_to_N0.csv  ->  compute N0  ->  compute R
################################################################################

ninoc_tbl <- readr::read_csv(NINOC_TABLE_CSV, show_col_types = FALSE)
names(ninoc_tbl) <- trimws(names(ninoc_tbl))

expected <- c("Taxon","Replicate","N_inoculation_cells_per_L","delta_Ninoc_to_N0_min")
if (!all(expected %in% names(ninoc_tbl)) && ncol(ninoc_tbl) >= 4) {
  names(ninoc_tbl)[1:4] <- expected
}
stopifnot(all(expected %in% names(ninoc_tbl)))

ninoc_tbl <- ninoc_tbl %>%
  dplyr::transmute(
    Taxon  = as.character(Taxon),
    Replicate = as.character(Replicate),
    N_inoculation_cells_per_L = as.numeric(N_inoculation_cells_per_L),
    delta_Ninoc_to_N0_min     = as.numeric(delta_Ninoc_to_N0_min)
  ) %>%
  dplyr::distinct()

resp <- results_fit %>%
  dplyr::left_join(ninoc_tbl, by = c("Taxon","Replicate")) %>%
  dplyr::mutate(
    N0_cells_per_L = dplyr::if_else(
      is.finite(N_inoculation_cells_per_L) & N_inoculation_cells_per_L > 0 &
        is.finite(delta_Ninoc_to_N0_min) & delta_Ninoc_to_N0_min >= 0 &
        is.finite(r_per_minute) & r_per_minute > 0,
      N_inoculation_cells_per_L * exp(r_per_minute * delta_Ninoc_to_N0_min),
      NA_real_
    ),
    biomass_integral_cells_min_per_L = dplyr::if_else(
      is.finite(N0_cells_per_L) & N0_cells_per_L > 0 &
        is.finite(r_per_minute) & r_per_minute > 0 &
        is.finite(T_end_min) & T_end_min > 0,
      N0_cells_per_L * (exp(r_per_minute * T_end_min) - 1) / r_per_minute,
      NA_real_
    ),
    R = dplyr::if_else(
      is.finite(C_tot_mg_per_L) & C_tot_mg_per_L > 0 &
        is.finite(biomass_integral_cells_min_per_L) & biomass_integral_cells_min_per_L > 0,
      C_tot_mg_per_L / biomass_integral_cells_min_per_L,
      NA_real_
    ),
    
    # ---- carbon fluxes per cell per hour (fg C cell^-1 h^-1) ----
    R_C_fg_cell_h = R * O2_to_C_mass * 1e12 * 60,          # mg O2 -> fg C, per hour
    G_C_fg_cell_h = r_per_minute * 60 * C_per_cell_fg,     # (h^-1) * (fg C cell^-1)
    
    r_per_hour = r_per_minute * 60
  )

readr::write_csv(resp, RESULTS_FINAL_CSV)

################################################################################
# ========================== FIGURES / ANALYSES ================================
################################################################################

# Reload fit_curves for safety in case you run sections separately
fit_curves <- readr::read_csv(FITCURVES_CSV, show_col_types = FALSE)
results    <- resp

################################################################################
# Supp Fig S1 — Residual diagnostics (TIFF)
################################################################################

resid_all <- fit_curves %>%
  dplyr::mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate),
    Time0     = Time0_min,
    Fitted    = Fit_norm,
    Residual  = Oxygen_norm - Fit_norm
  ) %>%
  dplyr::filter(is.finite(Time0), is.finite(Fitted), is.finite(Residual))

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
    labs(x = "Time (min, re-zeroed)", y = "Residual (normalised O\u2082)") +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.text      = element_text(face = "italic"),
      axis.title      = element_text(size = 11),
      axis.text       = element_text(size = 9)
    )
  
  p_resid_fit <- ggplot(resid_all, aes(x = Fitted, y = Residual, colour = Replicate)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_point(alpha = 0.7, size = 0.8) +
    facet_wrap(~ Taxon, scales = "free_x") +
    labs(x = "Fitted normalised O\u2082", y = "Residual (normalised O\u2082)") +
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
    filename = SUPP_S1_RESID_TIF,
    plot     = S1_residuals,
    width    = 12,
    height   = 10,
    dpi      = 600,
    device   = "tiff",
    compression = "lzw"
  )
}

################################################################################
# Supp Fig 3 — Normalised O2 dynamics, colour by Replicate (TIFF)
################################################################################

all_fits_df_norm <- fit_curves %>%
  dplyr::mutate(
    TaxonFull = as.character(Taxon),
    Replicate = factor(Replicate)
  ) %>%
  dplyr::select(TaxonFull, Replicate, Time0_min, Oxygen_norm, Fit_norm) %>%
  dplyr::rename(Time = Time0_min, Oxygen_n = Oxygen_norm, Pred_n = Fit_norm) %>%
  dplyr::arrange(TaxonFull, Replicate, Time)

supp_plot_norm_rep <- ggplot(all_fits_df_norm, aes(x = Time, y = Oxygen_n)) +
  geom_point(color = "grey60", size = 1, alpha = 0.55) +
  geom_line(aes(y = Pred_n, color = Replicate, group = Replicate),
            linewidth = 0.9, na.rm = TRUE) +
  facet_wrap(~ TaxonFull, scales = "free_y") +
  labs(
    title = "Supplementary: Normalised O₂ Dynamics (colour = Replicate)",
    x = "Time since fit start (minutes)",
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
  filename = SUPP_FIG3_NORM_TIF,
  plot     = supp_plot_norm_rep,
  width    = 14,
  height   = 10,
  dpi      = 600,
  device   = "tiff",
  compression = "lzw"
)

################################################################################
# Fig 2 — 4-panel facet plot with r and R labels (TIFF)
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

facet_data <- fit_curves %>%
  dplyr::mutate(
    Time = Time0_min,
    Predicted_O2 = Fit_norm,
    Oxygen = Oxygen_norm
  ) %>%
  dplyr::inner_join(selected_combos, by = c("Taxon","Replicate")) %>%
  dplyr::select(Taxon, Replicate, Time, Oxygen, Predicted_O2)

if (nrow(facet_data) > 0) {
  
  annot_info <- results %>%
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
          "italic(R)==", sci_pm(R, digits = 2),
          "~mg~O[2]~cell^{-1}~min^{-1}"
        )
      } else {
        paste0(
          "atop(",
          "italic(r)==", sprintf("%.2f", r_per_hour), "~h^{-1},",
          "italic(R)==", sci_pm(R, digits = 2),
          "~mg~O[2]~cell^{-1}~min^{-1})"
        )
      }
    ) %>%
    dplyr::select(Taxon, Replicate, FacetLabel, label_text)
  
  facet_data <- facet_data %>%
    dplyr::inner_join(annot_info, by = c("Taxon","Replicate"))
  
  label_positions <- facet_data %>%
    dplyr::group_by(FacetLabel) %>%
    dplyr::group_modify(~{
      df <- .x %>% dplyr::arrange(Time) %>% dplyr::distinct(Time, .keep_all = TRUE)
      xmin <- min(df$Time, na.rm = TRUE); xmax <- max(df$Time, na.rm = TRUE)
      ymin <- min(df$Oxygen, na.rm = TRUE); ymax <- max(df$Oxygen, na.rm = TRUE)
      xr <- xmax - xmin; yr <- ymax - ymin
      
      x_pos <- xmin + X_ANCHOR * xr
      x_pos <- max(xmin + X_MARGIN * xr, min(xmax - X_MARGIN * xr, x_pos))
      
      ok <- is.finite(df$Predicted_O2)
      curveY <- if (sum(ok) >= 2) approx(df$Time[ok], df$Predicted_O2[ok], xout = x_pos, rule = 2)$y
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
      x = "Time since fit start (minutes)",
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
    filename = FIG2_FACET4_TIF,
    plot     = facet_plot,
    width    = 14,
    height   = 4.5,
    dpi      = 600,
    device   = "tiff",
    compression = "lzw"
  )
}

################################################################################
# Fig 3 — Boxplot + mixed model (TIFF)
################################################################################

OD_FC_CSV <- file.path(data_dir, "OD_r_FC_r.csv")

od_fc_data <- readr::read_csv(OD_FC_CSV, show_col_types = FALSE) %>%
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

readr::write_csv(growth, file.path(tables_dir, "growth_rates_combined.csv"))

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
  filename = FIG3_GROWTH_TIF,
  plot     = combined_comparison,
  width    = 14,
  height   = 6,
  dpi      = 600,
  device   = "tiff",
  compression = "lzw"
)

mm_method   <- lmer(Growth_Rate ~ 1 + Method + (1 | Taxon), REML = FALSE, data = growth_long)
mm_method_1 <- lmer(Growth_Rate ~ 1 + (1 | Taxon), REML = FALSE, data = growth_long)

anova(mm_method, mm_method_1)

aic_comparison <- AIC(mm_method, mm_method_1)
print(aic_comparison)

print(summary(mm_method))
ci_method    <- confint(mm_method, method = "Wald")
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = METHOD_EFFECTS_TIF,
  plot     = method_comparison,
  width    = 8,
  height   = 6,
  dpi      = 600,
  device   = "tiff",
  compression = "lzw"
)

readr::write_csv(method_effects, file.path(tables_dir, "method_effects_estimates.csv"))
readr::write_csv(pairs,          file.path(tables_dir, "method_effects_significance.csv"))

capture.output(summary(mm_OD  <- lmer(r_O2 ~ r_OD600 + (1 | Taxon), data = growth)),
               file = file.path(tables_dir, "mixed_model_OD600_summary.txt"))
capture.output(summary(mm_FC  <- lmer(r_O2 ~ r_FC    + (1 | Taxon), data = growth)),
               file = file.path(tables_dir, "mixed_model_FC_summary.txt"))

################################################################################
# Fig 4 — LME regressions (TIFF)
################################################################################

create_norm_data <- function(model) {
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
  eq_text  <- sprintf("y = %.3f + %.3fx\nR² = %.3f", fixed_coef[1], fixed_coef[2], r2_mixed)
  
  list(data = growth_norm, fixed_coef = fixed_coef, eq_text = eq_text)
}

od_norm <- create_norm_data(mm_OD)
fc_norm <- create_norm_data(mm_FC)

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
  theme(legend.position = "top", legend.title = element_blank()) +
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
  filename = FIG4_REGRESS_TIF,
  plot     = combined_norm,
  width    = 10,
  height   = 16,
  dpi      = 600,
  device   = "tiff",
  compression = "lzw"
)

################################################################################
# Fig 5 — Bland–Altman (title WRAPPED so it fits) (TIFF)
################################################################################

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && is.finite(a)) a else b

df_raw <- readr::read_csv(file.path(tables_dir, "growth_rates_combined.csv"),
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
  
  # Wrapped title so it fits
  title_raw <- glue::glue("({panel_letter}) {label1} vs {label2} — Bland–Altman (all replicates)")
  title_wrapped <- stringr::str_wrap(title_raw, width = 55)
  
  ggplot(df_ba, aes(x = avg, y = diff)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper,
             fill = "#D55E00", alpha = 0.15) +
    geom_point(size = 2.5, alpha = 0.9, color = "black") +
    geom_hline(yintercept = bias,  color = "black",     linewidth = 0.7) +
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
      title = title_wrapped,
      x = "Mean growth rate (per replicate)",
      y = "Difference in growth rate (per replicate)"
    ) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(lineheight = 1.05))
}

plot1 <- bland_altman_plot_lme(df_raw, Oxygen_r, OD_r, "Oxygen", "OD600", "A")
plot2 <- bland_altman_plot_lme(df_raw, Oxygen_r, FC_r, "Oxygen", "Flow cytometry", "B")

combined_plot <- gridExtra::grid.arrange(plot1, plot2, ncol = 2)

ggsave(
  filename = FIG5_BA_TIF,
  plot     = combined_plot,
  width    = 14,
  height   = 6,
  dpi      = 600,
  device   = "tiff",
  compression = "lzw"
)

################################################################################
# Fig 6 — RIS + slope intensity (publication combined A/B using patchwork)
#  - Fig 6a (A): ONE-PANEL RIS (collapsed) + taxon-specific fitted lines
#  - Fig 6b (B): Taxon-specific slope “intensity” (forest/caterpillar)
#  - Combined export: Fig_6_COMBINED.tiff with tags (A) and (B)
#  - Same Taxon colours used in BOTH panels (locked via named palette)
#
# USER REQUESTED CHANGES:
# 1) Panel tags as (A) and (B) (parentheses)
# 2) Panel (A) annotation uses y = m x + n (NO (x - x̄) text)
# 3) Panel (B) x-axis text: "Taxon-specific slope" ONLY (no beta1+b1_taxon text)
################################################################################

# ---- OUTPUTS (add these paths) ----
FIG6B_SLOPES_TIF   <- file.path(plots_dir, "Fig_6b_taxon_specific_slopes.tiff")
FIG6_COMBINED_TIF  <- file.path(plots_dir, "Fig_6_COMBINED.tiff")

rm_df <- results %>%
  dplyr::filter(
    fit_ok,
    is.finite(G_C_fg_cell_h), G_C_fg_cell_h > 0,
    is.finite(R_C_fg_cell_h), R_C_fg_cell_h > 0
  ) %>%
  dplyr::mutate(
    Taxon = factor(Taxon),
    log_G = log10(G_C_fg_cell_h),
    log_R = log10(R_C_fg_cell_h)
  )

if (nrow(rm_df) < 8 || dplyr::n_distinct(rm_df$Taxon) < 2) {
  warning("Not enough data / taxa for Fig 6 random-slope RIS model.")
} else {
  
  # -----------------------------
  # 1) Center x for the mixed model (critical)
  # -----------------------------
  x0 <- mean(rm_df$log_R, na.rm = TRUE)
  
  rm_df2 <- rm_df %>%
    dplyr::mutate(log_Rc = log_R - x0)
  
  mm <- suppressWarnings(
    lmerTest::lmer(
      log_G ~ log_Rc + (1 + log_Rc | Taxon),
      data = rm_df2,
      REML = FALSE
    )
  )
  
  utils::capture.output(
    summary(mm),
    file = file.path(tables_dir, "mixed_model_Fig6_RIS_centeredX.txt")
  )
  
  fe <- lme4::fixef(mm)  # (Intercept), log_Rc
  
  re_tbl <- lme4::ranef(mm)$Taxon %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Taxon") %>%
    dplyr::transmute(
      Taxon = factor(Taxon, levels = levels(rm_df2$Taxon)),
      b0    = .data[["(Intercept)"]],
      b1    = .data[["log_Rc"]]
    )
  
  # -----------------------------
  # 2) Collapse: remove ONLY taxon intercept b0
  # -----------------------------
  rm_collapsed <- rm_df2 %>%
    dplyr::left_join(re_tbl, by = "Taxon") %>%
    dplyr::mutate(log_G_collapsed = log_G - b0)
  
  # -----------------------------
  # 3) Global fixed-effect line + 95% CI (fixed effects only)
  # -----------------------------
  x_seq <- seq(min(rm_df2$log_R, na.rm = TRUE),
               max(rm_df2$log_R, na.rm = TRUE),
               length.out = 250)
  x_seq_c <- x_seq - x0
  
  V  <- as.matrix(vcov(mm))
  X  <- cbind(1, x_seq_c)
  mu <- as.numeric(X %*% fe)
  se <- sqrt(pmax(0, rowSums((X %*% V) * X)))
  
  pred_global <- tibble::tibble(
    log_R = x_seq,
    mu    = mu,
    lo    = mu - 1.96 * se,
    hi    = mu + 1.96 * se
  )
  
  # -----------------------------
  # 4) Taxon-specific fitted lines: extend ranges + minimum span
  # -----------------------------
  EXPAND_FRAC <- 0.25
  MIN_SPAN    <- 0.35
  
  x_global_min <- min(rm_df2$log_R, na.rm = TRUE)
  x_global_max <- max(rm_df2$log_R, na.rm = TRUE)
  
  taxon_ranges <- rm_df2 %>%
    dplyr::group_by(Taxon) %>%
    dplyr::summarise(
      xmin = min(log_R, na.rm = TRUE),
      xmax = max(log_R, na.rm = TRUE),
      xmid = mean(log_R, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      span  = pmax(xmax - xmin, MIN_SPAN),
      xmin2 = pmax(x_global_min, xmid - 0.5 * span),
      xmax2 = pmin(x_global_max, xmid + 0.5 * span),
      xmin3 = pmax(x_global_min, xmin2 - EXPAND_FRAC * span),
      xmax3 = pmin(x_global_max, xmax2 + EXPAND_FRAC * span)
    )
  
  taxon_lines <- re_tbl %>%
    dplyr::left_join(taxon_ranges, by = "Taxon") %>%
    dplyr::mutate(
      slope_tax           = unname(fe[2]) + b1,
      intercept_collapsed = unname(fe[1])
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(log_R_seq = list(seq(xmin3, xmax3, length.out = 120))) %>%
    dplyr::ungroup() %>%
    tidyr::unnest(log_R_seq) %>%
    dplyr::mutate(log_Rc_seq = log_R_seq - x0) %>%
    dplyr::transmute(
      Taxon,
      log_R = log_R_seq,
      mu_tax_collapsed = intercept_collapsed + slope_tax * log_Rc_seq
    )
  
  # -----------------------------
  # 5) Colours (same mapping used for Fig 6a + Fig 6b)
  # -----------------------------
  leg_nrow <- 3
  n_tax    <- nlevels(rm_df2$Taxon)
  leg_ncol <- ceiling(n_tax / leg_nrow)
  
  okabe_ito_extended_distinct <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000",
    "#332288", "#88CCEE", "#117733", "#999933",
    "#CC6677", "#AA4499", "#661100", "#BBBBBB"
  )
  
  pal <- if (n_tax <= length(okabe_ito_extended_distinct)) {
    okabe_ito_extended_distinct[seq_len(n_tax)]
  } else {
    rep(okabe_ito_extended_distinct, length.out = n_tax)
  }
  
  pal_named <- stats::setNames(pal, levels(rm_df2$Taxon))
  
  isme_theme_pub <- function() {
    ggplot2::theme_classic(base_size = 14) +
      ggplot2::theme(
        axis.title    = ggplot2::element_text(size = 15),
        axis.text     = ggplot2::element_text(size = 13),
        axis.line     = ggplot2::element_line(linewidth = 0.9),
        axis.ticks    = ggplot2::element_line(linewidth = 0.7),
        legend.position   = "top",
        legend.direction  = "horizontal",
        legend.justification = "center",
        legend.title      = ggplot2::element_blank(),
        legend.text       = ggplot2::element_text(face = "italic", size = 11),
        legend.key.width  = grid::unit(1.00, "lines"),
        legend.key.height = grid::unit(0.75, "lines"),
        legend.spacing.x  = grid::unit(0.30, "lines"),
        legend.box        = "vertical",
        legend.margin     = ggplot2::margin(b = 6, unit = "pt"),
        plot.margin = grid::unit(c(10, 42, 10, 22), "pt")
      )
  }
  
  r2_mixed <- cor(stats::fitted(mm), rm_df2$log_G, use = "complete.obs")^2
  
  # For annotation: express the fitted fixed-effects line back on UNCENTERED x:
  # log_G = (beta0 - beta1*x0) + beta1*log_R  => y = m x + n
  m_fixed <- unname(fe[2])
  n_fixed <- unname(fe[1] - fe[2] * x0)
  
  ann_text <- sprintf(
    "Fixed effect: y = %.2f x + %.2f\nFixed slope = %.2f   R² (fitted vs observed) = %.2f",
    m_fixed, n_fixed, m_fixed, r2_mixed
  )
  
  # -----------------------------
  # 6a) Fig 6a — ONE-PANEL RIS (collapsed)
  # -----------------------------
  p6 <- ggplot2::ggplot(rm_collapsed, ggplot2::aes(x = log_R, y = log_G_collapsed, colour = Taxon)) +
    ggplot2::geom_ribbon(
      data = pred_global,
      ggplot2::aes(x = log_R, ymin = lo, ymax = hi),
      inherit.aes = FALSE,
      alpha = 0.06,
      colour = NA
    ) +
    ggplot2::geom_line(
      data = pred_global,
      ggplot2::aes(x = log_R, y = mu),
      inherit.aes = FALSE,
      linewidth = 1.15,
      colour = "black"
    ) +
    ggplot2::geom_point(size = 2.2, alpha = 0.55) +
    ggplot2::geom_line(
      data = taxon_lines,
      ggplot2::aes(x = log_R, y = mu_tax_collapsed, colour = Taxon, group = Taxon),
      inherit.aes = FALSE,
      linewidth = 1.55,
      alpha = 0.98
    ) +
    ggplot2::scale_color_manual(values = pal_named, labels = function(x) gsub("_", " ", x)) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        nrow = leg_nrow, ncol = leg_ncol, byrow = TRUE,
        override.aes = list(size = 3.1, alpha = 1, linewidth = 1.2)
      )
    ) +
    ggplot2::labs(
      x = expression(log[10](R~"(fg C cell"^{-1}~" h"^{-1}*")")),
      y = expression("Taxon-collapsed " * log[10](G~"(fg C cell"^{-1}~" h"^{-1}*")"))
    ) +
    ggplot2::annotate(
      "text",
      x = min(rm_collapsed$log_R, na.rm = TRUE) + 0.03 * diff(range(rm_collapsed$log_R, na.rm = TRUE)),
      y = max(rm_collapsed$log_G_collapsed, na.rm = TRUE) - 0.05 * diff(range(rm_collapsed$log_G_collapsed, na.rm = TRUE)),
      label = ann_text,
      hjust = 0, vjust = 1, size = 4.1, lineheight = 1.05
    ) +
    isme_theme_pub()
  
  # Save Fig 6a (optional)
  ggsave(
    filename = FIG6_RIS_MAIN_TIF,
    plot     = p6,
    width    = 9.0,
    height   = 6.9 + 0.55 * (leg_nrow - 1),
    dpi      = 600,
    device   = "tiff",
    compression = "lzw"
  )
  
  readr::write_csv(rm_collapsed, file.path(tables_dir, "data_Fig6_RIS_COLLAPSED_centeredX.csv"))
  readr::write_csv(taxon_lines,  file.path(tables_dir, "Fig6_taxon_specific_lines_COLLAPSED_centeredX_EXTENDED.csv"))
  
  # -----------------------------
  # 6b) Fig 6b — Taxon-specific slope “intensity”
  # -----------------------------
  beta1    <- lme4::fixef(mm)[["log_Rc"]]
  se_beta1 <- sqrt(as.matrix(vcov(mm))[["log_Rc", "log_Rc"]])
  
  re_tax <- lme4::ranef(mm, condVar = TRUE)$Taxon
  pv     <- attr(re_tax, "postVar")
  
  taxa  <- rownames(re_tax)
  b1    <- re_tax[, "log_Rc"]
  se_b1 <- sqrt(pmax(0, pv[2, 2, ]))
  
  slopes_df <- tibble::tibble(
    Taxon     = factor(taxa, levels = levels(rm_df2$Taxon)),
    slope_tax = as.numeric(beta1 + b1),
    se_tax    = sqrt(se_beta1^2 + se_b1^2),
    lo        = slope_tax - 1.96 * se_tax,
    hi        = slope_tax + 1.96 * se_tax
  ) %>%
    dplyr::filter(!is.na(Taxon)) %>%
    dplyr::arrange(slope_tax) %>%
    dplyr::mutate(Taxon_ord = factor(as.character(Taxon), levels = as.character(Taxon)))
  
  readr::write_csv(slopes_df, file.path(tables_dir, "Fig6b_taxon_specific_slopes.csv"))
  
  p6b <- ggplot2::ggplot(slopes_df, ggplot2::aes(x = slope_tax, y = Taxon_ord, colour = Taxon)) +
    ggplot2::geom_vline(xintercept = beta1, linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0.0, linewidth = 0.9) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::scale_color_manual(values = pal_named, guide = "none") +
    ggplot2::labs(
      x = "Taxon-specific slope",
      y = NULL
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_text(face = "italic", size = 11),
      axis.title.x = ggplot2::element_text(size = 15)
    )
  
  # Save Fig 6b (optional)
  ggplot2::ggsave(
    filename = FIG6B_SLOPES_TIF,
    plot     = p6b,
    width    = 7.8,
    height   = max(4.8, 0.28 * nrow(slopes_df) + 1.2),
    dpi      = 600,
    device   = "tiff",
    compression = "lzw"
  )
  
  # -----------------------------
  # 6 combined (A/B) — publication combined figure with tags "(A)" "(B)"
  # -----------------------------
  fig6_combined <- (p6 + p6b) +
    patchwork::plot_layout(widths = c(2.2, 1)) +
    patchwork::plot_annotation(tag_levels = list(c("(A)", "(B)")))
  
  ggplot2::ggsave(
    filename = FIG6_COMBINED_TIF,
    plot     = fig6_combined,
    width    = 16,
    height   = 7.2,
    dpi      = 600,
    device   = "tiff",
    compression = "lzw"
  )
}

