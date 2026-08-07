# =============================================================================
# 05_oxygen_fits.R - Fit the normalised O2 model + reconstruct N0 and R
# =============================================================================
# Per (Taxon, Replicate) series: apply any manual fit window / exclusion set in
# the 03 trim-selector app (loaded by config), normalise O2, fit
#   O2_norm(t) = O2_0 + (K / r) * (1 - exp(r * t))   via minpack.lm::nlsLM,
# run the QC gate, export fitted curves, then reconstruct N0 by the
# depletion-anchored method (FC_Final projected back to the fit start over the
# fit-start -> 90% O2 depletion interval; config N0_METHOD), and per-cell
# respiration R plus carbon fluxes (using per-taxon cell carbon).
#
#   Input : results/tables/Oxygen_Data_Filtered.csv            (from 02)
#           results/tables/Oxygen_All_Long.csv                 (full trace -> depletion time)
#           data/OD_r_FC_r.csv (FC_Final); config: N0_METHOD, cell carbon,
#                   MANUAL_FIT_WINDOWS, PLOT_EXCLUDE_POINTS
#   Output: results/tables/oxygen_fit_results.csv
#           results/tables/oxygen_fit_curves.csv
#           results/tables/oxygen_results_with_R.csv
#   Figure: results/figures/oxygen_model_fit_curves.pdf        (multi-page)
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
suppressPackageStartupMessages({ library(zoo); library(minpack.lm) })

FITS_PDF <- fig("oxygen_model_fit_curves.pdf")

# ---- load filtered oxygen (from 02); apply taxon + curve exclusions ----------
oxygen_data <- readr::read_csv(FILTERED_CSV, show_col_types = FALSE) %>%
  dplyr::mutate(Taxon = as.character(Taxon), Replicate = as.character(Replicate)) %>%
  dplyr::arrange(Taxon, Replicate, Time)
oxygen_data <- drop_excluded_taxa(oxygen_data, "filtered")
oxygen_data <- apply_excluded_curves(oxygen_data, "filtered curves")

# ---- per-curve delay (inoculation -> O2 start) from 02 metadata --------------
delta_tbl <- NULL
if (file.exists(TRIM_META_CSV)) {
  delta_tbl <- tryCatch(readr::read_csv(TRIM_META_CSV, show_col_types = FALSE) %>%
    dplyr::transmute(Taxon = as.character(Taxon), Replicate = as.character(Replicate),
                     delta_Ninoc_to_N0_min = as.numeric(delta_Ninoc_to_N0_min)),
    error = function(e) NULL)
}
delta_of <- function(tax, rep) {
  if (is.null(delta_tbl)) return(NA_real_)
  hit <- delta_tbl$delta_Ninoc_to_N0_min[delta_tbl$Taxon == tax & delta_tbl$Replicate == rep]
  if (length(hit)) hit[1] else NA_real_
}

# ---- optional per-curve inoculation table (overrides the config scalar) ------
ninoc_tbl <- load_ninoc_table()
ninoc_of <- function(tax, rep) {
  # returns c(N_inoc, delta_override) for this curve; NA delta = use 02's delay.
  if (is.null(ninoc_tbl)) return(c(N_inoculation_cells_per_L, NA_real_))
  hit <- ninoc_tbl[ninoc_tbl$Taxon == tax & ninoc_tbl$Replicate == rep, , drop = FALSE]
  if (nrow(hit) == 0) return(c(N_inoculation_cells_per_L, NA_real_))
  n <- hit$N_inoculation_cells_per_L[1]
  c(if (is.finite(n) && n > 0) n else N_inoculation_cells_per_L, hit$delta_Ninoc_to_N0_min[1])
}

# ---- manual fit-window lookup (from the 03 app, loaded by config) -----------
mw_of <- function(tax, rep) {
  if (is.null(MANUAL_FIT_WINDOWS) || nrow(MANUAL_FIT_WINDOWS) == 0) return(c(NA_real_, NA_real_))
  hit <- MANUAL_FIT_WINDOWS[as.character(MANUAL_FIT_WINDOWS$Taxon) == tax &
                              toupper(as.character(MANUAL_FIT_WINDOWS$Replicate)) == toupper(rep), , drop = FALSE]
  if (nrow(hit) == 0) return(c(NA_real_, NA_real_))
  c(suppressWarnings(as.numeric(hit$fit_start[1])), suppressWarnings(as.numeric(hit$fit_end[1])))
}

results_fit <- tibble::tibble(
  Taxon = character(), Replicate = character(),
  r_per_minute = numeric(), K = numeric(), O2_0 = numeric(), O2_ref = numeric(),
  T_end_min = numeric(), fit_start_min = numeric(), C_tot_mg_per_L = numeric(),
  pseudo_R2 = numeric(), AIC = numeric(), fit_ok = logical()
)
fit_curves_lst <- list()

combos <- oxygen_data %>% dplyr::group_by(Taxon, Replicate) %>% dplyr::group_keys()

grDevices::pdf(FITS_PDF, width = 7.2, height = 5.2, useDingbats = FALSE)
on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)

for (i in seq_len(nrow(combos))) {
  Tax <- combos$Taxon[i]; Rep <- combos$Replicate[i]

  df0 <- oxygen_data %>% dplyr::filter(Taxon == Tax, Replicate == Rep) %>% dplyr::arrange(Time)
  if (nrow(df0) < 5) next

  mw <- mw_of(Tax, Rep)
  if (is.finite(mw[1]) || is.finite(mw[2])) {
    # Manual fit window from the app: use it directly (no derivative onset-trim).
    dfw <- df0
    if (is.finite(mw[1])) dfw <- dfw %>% dplyr::filter(Time >= mw[1])
    if (is.finite(mw[2])) dfw <- dfw %>% dplyr::filter(Time <= mw[2])
    if (nrow(dfw) < 5) next
    df <- dfw %>% dplyr::mutate(Time0 = Time - min(Time, na.rm = TRUE))
  } else {
    # Automatic onset trim via derivative.
    df0 <- df0 %>% dplyr::mutate(dO2 = c(NA, diff(Oxygen)), dt = c(NA, diff(Time)),
                                 dO2dt = dO2 / dt, sm = zoo::rollmean(dO2dt, 3, fill = NA, align = "right"))
    idx <- which(df0$sm < -1e-7)[1]
    idx <- if (is.na(idx) || idx > nrow(df0) - 2) 1 else min(idx + 5, nrow(df0) - 2)
    df <- df0[idx:nrow(df0), ] %>% dplyr::mutate(Time0 = Time - min(Time, na.rm = TRUE))
  }

  O0 <- mean(head(df$Oxygen, 3), na.rm = TRUE)
  if (!is.finite(O0) || O0 <= 0) next
  df <- df %>% dplyr::mutate(Oxygen_norm = Oxygen / O0)

  # Optional fractional-drawdown window standardisation.
  if (isTRUE(USE_DRAWDOWN_WINDOW) && nrow(df) >= 5) {
    y1 <- df$Oxygen_norm[1]; ymin <- min(df$Oxygen_norm, na.rm = TRUE)
    if (is.finite(y1) && is.finite(ymin) && y1 > ymin) {
      target <- y1 - FIT_DRAWDOWN_FRAC * (y1 - ymin)
      cut <- which(df$Oxygen_norm <= target)[1]
      if (is.finite(cut) && cut >= 5) df <- df[1:cut, ]
    }
  }

  T_end <- max(df$Time0, na.rm = TRUE)
  if (!is.finite(T_end) || T_end <= 0) next
  fit_start_min <- min(df$Time, na.rm = TRUE)   # time (from t0) of the FIRST point in the model

  r_start <- {
    seg <- head(df, max(3, floor(0.3 * nrow(df))))
    slopes <- abs(diff(log(pmax(seg$Oxygen_norm, 1e-6))) / diff(seg$Time0))
    pmin(pmax(max(slopes, na.rm = TRUE), 1e-4), 5e-2)
  }
  K_start <- {
    slope <- suppressWarnings(min(diff(df$Oxygen_norm) / diff(df$Time0), na.rm = TRUE))
    if (!is.finite(slope)) slope <- -1e-3
    pmin(pmax(abs(slope), 1e-5), 0.1)
  }

  fit <- tryCatch(
    minpack.lm::nlsLM(Oxygen_norm ~ resp_model(r, K, Time0, O2_0), data = df,
                      start = list(r = r_start, K = K_start, O2_0 = 1),
                      lower = FIT_LOWER, upper = FIT_UPPER,
                      control = minpack.lm::nls.lm.control(maxiter = 300)),
    error = function(e) NULL)
  if (is.null(fit)) next

  pars <- summary(fit)$coefficients
  r_est <- pars["r", "Estimate"]; K_est <- pars["K", "Estimate"]; O2_0_est <- pars["O2_0", "Estimate"]
  pred <- stats::predict(fit, df); resid <- df$Oxygen_norm - pred
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
    dplyr::add_row(Taxon = Tax, Replicate = Rep, r_per_minute = r_est, K = K_est,
                   O2_0 = O2_0_est, O2_ref = O0, T_end_min = T_end, fit_start_min = fit_start_min,
                   C_tot_mg_per_L = C_tot,
                   pseudo_R2 = pseudo_R2, AIC = stats::AIC(fit), fit_ok = fit_ok)

  curve_df <- tibble::tibble(
    Taxon = Tax, Replicate = Rep, Time_min = df$Time, Time0_min = df$Time0,
    Oxygen_raw = df$Oxygen, Oxygen_norm = df$Oxygen_norm, Fit_norm = pred, Fit_raw = pred * O0,
    r_per_minute = r_est, K = K_est, O2_0 = O2_0_est, O2_ref = O0,
    pseudo_R2 = pseudo_R2, AIC = stats::AIC(fit), fit_ok = fit_ok)
  fit_curves_lst[[paste(Tax, Rep, sep = " | ")]] <- curve_df

  df_plot_all <- df0 %>% dplyr::mutate(Time0 = Time - min(df$Time, na.rm = TRUE))
  p_fit <- ggplot2::ggplot() +
    ggplot2::geom_line(data = df_plot_all, ggplot2::aes(Time0, Oxygen), colour = "grey70", linewidth = 0.9) +
    ggplot2::geom_point(data = df, ggplot2::aes(Time0, Oxygen), size = 1.9, alpha = 0.9) +
    ggplot2::geom_line(data = curve_df, ggplot2::aes(Time0_min, Fit_raw), colour = "black", linewidth = 1.1) +
    ggplot2::labs(title = paste0(Tax, " | ", Rep),
                  subtitle = sprintf("fit_ok=%s   r=%.4g min^-1   K=%.4g   R2=%.3f", fit_ok, r_est, K_est, pseudo_R2),
                  x = "Time since fit start (min)", y = expression(O[2] ~ "(mg L"^{-1} * ")")) +
    ggplot2::theme_classic(base_size = 11)
  print(p_fit)
}

try(grDevices::dev.off(), silent = TRUE)
readr::write_csv(results_fit, RESULTS_FIT_CSV)
fit_curves <- dplyr::bind_rows(fit_curves_lst)
readr::write_csv(fit_curves, FITCURVES_CSV)

# =============================================================================
# Reconstruct N0 (depletion-anchored FC_Final; config N0_METHOD) -> per-cell R
# =============================================================================
# DEFAULT (N0_METHOD = "depletion"): anchor N0 to the FINAL flow-cytometry count
# (FC_Final), projected back to the fit start over the interval during which growth
# actually happened (fit start -> 90% O2 depletion). FALLBACK ("initial"): the old
# N0 = N_inoc * exp(r * delta) route. r and the fit windows are untouched either way.
.dep <- load_depletion_table(); .fcf <- load_fc_final()
if (identical(N0_METHOD, "depletion") && !is.null(.dep) && !is.null(.fcf)) {
  resp <- results_fit %>%
    dplyr::left_join(.dep, by = c("Taxon", "Replicate")) %>%
    dplyr::left_join(.fcf, by = c("Taxon", "Replicate")) %>%
    dplyr::mutate(
      cell_carbon_fg = cell_carbon_of(Taxon),
      t_growth_min   = t_depletion_min - fit_start_min,   # fit start -> O2 depletion
      N0_cells_per_L = dplyr::if_else(
        is.finite(FC_Final) & FC_Final > 0 & is.finite(r_per_minute) & r_per_minute > 0 &
          is.finite(t_growth_min) & t_growth_min > 0,
        n0_depletion(FC_Final, r_per_minute, t_depletion_min, fit_start_min), NA_real_))
  .nmiss <- sum(!resp$depleted, na.rm = TRUE)
  if (.nmiss > 0) message("05_oxygen_fits: ", .nmiss, " curve(s) did not reach ",
    round(100 * (1 - DEPLETION_FRAC)), "% depletion; used recording end as growth-stop.")
} else {
  .taxa  <- results_fit$Taxon; .reps <- results_fit$Replicate
  .nin   <- if (length(.taxa)) t(vapply(seq_along(.taxa),
               function(i) ninoc_of(.taxa[i], .reps[i]), numeric(2))) else matrix(numeric(0), 0, 2)
  .dtrim <- results_fit$fit_start_min
  resp <- results_fit %>%
    dplyr::mutate(
      N_inoc_cells_per_L    = .nin[, 1],
      delta_Ninoc_to_N0_min = dplyr::if_else(is.finite(.nin[, 2]), .nin[, 2], .dtrim + INOC_DELAY_MIN),
      cell_carbon_fg        = cell_carbon_of(Taxon),
      N0_cells_per_L = dplyr::case_when(
        !isTRUE(N0_BACKPROJECT) ~ N_inoc_cells_per_L,
        is.finite(delta_Ninoc_to_N0_min) & is.finite(r_per_minute) & r_per_minute > 0 ~
          N_inoc_cells_per_L * exp(r_per_minute * delta_Ninoc_to_N0_min),
        TRUE ~ N_inoc_cells_per_L))
}

resp <- resp %>%
  dplyr::mutate(
    biomass_integral_cells_min_per_L = dplyr::if_else(
      is.finite(N0_cells_per_L) & N0_cells_per_L > 0 & is.finite(r_per_minute) &
        r_per_minute > 0 & is.finite(T_end_min) & T_end_min > 0,
      N0_cells_per_L * (exp(r_per_minute * T_end_min) - 1) / r_per_minute, NA_real_),
    R = dplyr::if_else(
      is.finite(C_tot_mg_per_L) & C_tot_mg_per_L > 0 &
        is.finite(biomass_integral_cells_min_per_L) & biomass_integral_cells_min_per_L > 0,
      C_tot_mg_per_L / biomass_integral_cells_min_per_L, NA_real_),
    R_C_fg_cell_h = R * O2_to_C_mass * MG_TO_FG * MIN_TO_H,
    G_C_fg_cell_h = r_per_minute * MIN_TO_H * cell_carbon_fg,
    r_per_hour    = r_per_minute * MIN_TO_H
  )

readr::write_csv(resp, RESULTS_FINAL_CSV)
message("05_oxygen_fits: ", sum(results_fit$fit_ok), "/", nrow(results_fit),
        " good fits. N0 method = ", N0_METHOD, ". Curves -> ", FITS_PDF)
