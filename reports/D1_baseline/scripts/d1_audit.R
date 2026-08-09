# =============================================================================
# d1_audit.R - PART E: data integrity, the FC/Ninoc cross-check, the QC gate
# =============================================================================
# Read-only over data/ and results/, plus this run's runs/D1_baseline/tables.
# Writes to runs/D1_baseline/comparisons/. Changes no analysis decision.
#
# Run:  Rscript reports/D1_baseline/scripts/d1_audit.R
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(minpack.lm); library(zoo) })

BASE <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", "..", ".."),
  mustWork = FALSE)
if (!dir.exists(file.path(BASE, "scripts"))) BASE <- getwd()

DATA <- file.path(BASE, "data")
RUN  <- file.path(BASE, "runs", "D1_baseline", "tables")
OUT  <- file.path(BASE, "runs", "D1_baseline", "comparisons")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

rd <- function(p) suppressWarnings(suppressMessages(
  readr::read_csv(p, show_col_types = FALSE, progress = FALSE)))

# ---------------------------------------------------------------------------
# E1. Audit every file in data/
# ---------------------------------------------------------------------------
CONSUMERS <- tribble(
  ~file,                            ~consumed_by,                                          ~units,
  "Oxygen_Data_Long.csv",           "01_longdata.R -> 02 -> 05/07 (config OXYGEN_LONG_CSV)", "Time: min; Oxygen: mg O2 / L",
  "Oxygen_Data_Filtered_CUE.csv",   "10_temperature_cue.R only (config CUE_CSV)",           "Temperature: deg C; Time: min; Oxygen: mg O2 / L",
  "Ninoc.csv",                      "05_oxygen_fits.R via config::load_ninoc_table()",      "N_inoculation_cells_per_L: cells / L; delta: min",
  "taxon_cell_sizes.csv",           "config::cell_carbon_of() -> 05_oxygen_fits.R",         "um, um^3, fg C / um^3, fg C / cell",
  "OD_r_FC_r.csv",                  "06_main_figures.R (Figs 3/4/5, method effects)",       "OD600: absorbance; FC: cells / uL (nominal); Time: min",
  "Cell_Counts.csv",                "04_experiment_inputs.R app only, as a taxa fallback",  "OD600: absorbance; FC: cells / uL (nominal)"
)

data_files <- sort(basename(Sys.glob(file.path(DATA, "*.csv"))))
audit <- purrr::map_dfr(data_files, function(f) {
  x <- rd(file.path(DATA, f))
  miss <- vapply(x, function(cc) sum(is.na(cc)), integer(1))
  tibble(
    file = f,
    n_rows = nrow(x), n_cols = ncol(x),
    columns = paste(names(x), collapse = ", "),
    n_missing_cells = sum(miss),
    cols_with_missing = if (any(miss > 0))
      paste(sprintf("%s(%d)", names(miss)[miss > 0], miss[miss > 0]), collapse = ", ") else "none",
    n_taxa = if ("Taxon" %in% names(x)) dplyr::n_distinct(x$Taxon) else NA_integer_,
    n_series = if (all(c("Taxon", "Replicate") %in% names(x)))
      dplyr::n_distinct(paste(x$Taxon, x$Replicate)) else NA_integer_,
    bytes = file.info(file.path(DATA, f))$size
  )
}) %>% dplyr::left_join(CONSUMERS, by = "file") %>%
  dplyr::mutate(consumed_by = dplyr::coalesce(consumed_by, "NOT CONSUMED BY ANY SCRIPT"),
                units = dplyr::coalesce(units, "unknown"))
readr::write_csv(audit, file.path(OUT, "E1_data_audit.csv"))

# ---------------------------------------------------------------------------
# E2. FC_Initial (data/Cell_Counts.csv and data/OD_r_FC_r.csv)
#     vs N_inoculation_cells_per_L (data/Ninoc.csv)
#     Same quantity, up to three sources. Report only; D2 acts on it.
# ---------------------------------------------------------------------------
cc  <- rd(file.path(DATA, "Cell_Counts.csv")) %>%
  dplyr::transmute(Taxon = as.character(Taxon), Replicate = as.character(Replicate),
                   FC_Initial_CellCounts = as.numeric(FC_Initial))
odfc <- rd(file.path(DATA, "OD_r_FC_r.csv")) %>%
  dplyr::transmute(Taxon = as.character(Taxon), Replicate = as.character(Replicate),
                   FC_Initial_ODFC = as.numeric(FC_Initial))
nin <- rd(file.path(DATA, "Ninoc.csv")) %>%
  dplyr::transmute(Taxon = as.character(Taxon), Replicate = as.character(Replicate),
                   N_inoc_cells_per_L = as.numeric(N_inoculation_cells_per_L),
                   delta_min = as.numeric(delta_Ninoc_to_N0_min))

# The 04 app's stated convention is counts as cells/uL, converted with 1e6 to cells/L.
FC_TO_CELLS_PER_L <- 1e6

fc_ratio <- nin %>%
  dplyr::full_join(cc,   by = c("Taxon", "Replicate")) %>%
  dplyr::full_join(odfc, by = c("Taxon", "Replicate")) %>%
  dplyr::mutate(
    FC_CellCounts_as_cells_per_L = FC_Initial_CellCounts * FC_TO_CELLS_PER_L,
    FC_ODFC_as_cells_per_L       = FC_Initial_ODFC * FC_TO_CELLS_PER_L,
    ratio_Ninoc_over_FC_CellCounts = N_inoc_cells_per_L / FC_CellCounts_as_cells_per_L,
    ratio_Ninoc_over_FC_ODFC       = N_inoc_cells_per_L / FC_ODFC_as_cells_per_L,
    FC_sources_agree = dplyr::case_when(
      is.na(FC_Initial_CellCounts) | is.na(FC_Initial_ODFC) ~ NA,
      TRUE ~ FC_Initial_CellCounts == FC_Initial_ODFC)
  ) %>%
  dplyr::arrange(Taxon, Replicate)
readr::write_csv(fc_ratio, file.path(OUT, "E2_FC_vs_Ninoc_per_replicate.csv"))

fc_ratio_taxon <- fc_ratio %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(
    n = dplyr::n(),
    n_FC_sources_disagree = sum(FC_sources_agree == FALSE, na.rm = TRUE),
    median_ratio_vs_CellCounts = stats::median(ratio_Ninoc_over_FC_CellCounts, na.rm = TRUE),
    min_ratio_vs_CellCounts    = suppressWarnings(min(ratio_Ninoc_over_FC_CellCounts, na.rm = TRUE)),
    max_ratio_vs_CellCounts    = suppressWarnings(max(ratio_Ninoc_over_FC_CellCounts, na.rm = TRUE)),
    median_ratio_vs_ODFC       = stats::median(ratio_Ninoc_over_FC_ODFC, na.rm = TRUE),
    .groups = "drop")
readr::write_csv(fc_ratio_taxon, file.path(OUT, "E2_FC_vs_Ninoc_by_taxon.csv"))

# ---------------------------------------------------------------------------
# E3. The QC gate: which of the criteria does each series pass, and by how much?
#     Re-runs the 05 fit per curve purely to recover the per-criterion values;
#     the fit itself is identical to the pipeline's (same config, same windows).
# ---------------------------------------------------------------------------
source(file.path(BASE, "scripts", "config.R"))

filtered <- rd(file.path(RUN, "Oxygen_Data_Filtered.csv")) %>%
  dplyr::mutate(Taxon = as.character(Taxon), Replicate = as.character(Replicate)) %>%
  dplyr::arrange(Taxon, Replicate, Time)

mw_of <- function(tax, rep) {
  if (is.null(MANUAL_FIT_WINDOWS) || nrow(MANUAL_FIT_WINDOWS) == 0) return(c(NA_real_, NA_real_))
  hit <- MANUAL_FIT_WINDOWS[as.character(MANUAL_FIT_WINDOWS$Taxon) == tax &
                              toupper(as.character(MANUAL_FIT_WINDOWS$Replicate)) == toupper(rep), , drop = FALSE]
  if (nrow(hit) == 0) return(c(NA_real_, NA_real_))
  c(suppressWarnings(as.numeric(hit$fit_start[1])), suppressWarnings(as.numeric(hit$fit_end[1])))
}

combos <- filtered %>% dplyr::group_by(Taxon, Replicate) %>% dplyr::group_keys()
qc <- purrr::map_dfr(seq_len(nrow(combos)), function(i) {
  Tax <- combos$Taxon[i]; Rep <- combos$Replicate[i]
  df0 <- filtered %>% dplyr::filter(Taxon == Tax, Replicate == Rep) %>% dplyr::arrange(Time)
  if (nrow(df0) < 5) return(tibble(Taxon = Tax, Replicate = Rep, fitted = FALSE))
  mw <- mw_of(Tax, Rep)
  if (is.finite(mw[1]) || is.finite(mw[2])) {
    dfw <- df0
    if (is.finite(mw[1])) dfw <- dfw %>% dplyr::filter(Time >= mw[1])
    if (is.finite(mw[2])) dfw <- dfw %>% dplyr::filter(Time <= mw[2])
    if (nrow(dfw) < 5) return(tibble(Taxon = Tax, Replicate = Rep, fitted = FALSE))
    df <- dfw %>% dplyr::mutate(Time0 = Time - min(Time, na.rm = TRUE))
  } else {
    df0 <- df0 %>% dplyr::mutate(dO2 = c(NA, diff(Oxygen)), dt = c(NA, diff(Time)),
                                 dO2dt = dO2 / dt, sm = zoo::rollmean(dO2dt, 3, fill = NA, align = "right"))
    idx <- which(df0$sm < -1e-7)[1]
    idx <- if (is.na(idx) || idx > nrow(df0) - 2) 1 else min(idx + 5, nrow(df0) - 2)
    df <- df0[idx:nrow(df0), ] %>% dplyr::mutate(Time0 = Time - min(Time, na.rm = TRUE))
  }
  O0 <- mean(head(df$Oxygen, 3), na.rm = TRUE)
  df <- df %>% dplyr::mutate(Oxygen_norm = Oxygen / O0)
  r_start <- { seg <- head(df, max(3, floor(0.3 * nrow(df))))
    s <- abs(diff(log(pmax(seg$Oxygen_norm, 1e-6))) / diff(seg$Time0))
    pmin(pmax(max(s, na.rm = TRUE), 1e-4), 5e-2) }
  K_start <- { s <- suppressWarnings(min(diff(df$Oxygen_norm) / diff(df$Time0), na.rm = TRUE))
    if (!is.finite(s)) s <- -1e-3; pmin(pmax(abs(s), 1e-5), 0.1) }
  fit <- tryCatch(minpack.lm::nlsLM(Oxygen_norm ~ resp_model(r, K, Time0, O2_0), data = df,
      start = list(r = r_start, K = K_start, O2_0 = 1), lower = FIT_LOWER, upper = FIT_UPPER,
      control = minpack.lm::nls.lm.control(maxiter = 300)), error = function(e) NULL)
  if (is.null(fit)) return(tibble(Taxon = Tax, Replicate = Rep, fitted = FALSE))

  pars <- summary(fit)$coefficients
  r_est <- pars["r", "Estimate"]; K_est <- pars["K", "Estimate"]
  pred <- stats::predict(fit, df); resid <- df$Oxygen_norm - pred
  pr2 <- 1 - sum(resid^2, na.rm = TRUE) /
    sum((df$Oxygen_norm - mean(df$Oxygen_norm, na.rm = TRUE))^2, na.rm = TRUE)

  tibble(
    Taxon = Tax, Replicate = Rep, fitted = TRUE, n_points = nrow(df),
    r_per_minute = r_est, K = K_est, pseudo_R2 = pr2,
    max_rel_SE   = max(abs(pars[, "Std. Error"] / pars[, "Estimate"])),
    max_pval     = max(pars[, "Pr(>|t|)"]),
    resid_range  = diff(range(resid, na.rm = TRUE)),
    aic_gain     = stats::AIC(stats::lm(Oxygen_norm ~ 1, data = df)) - stats::AIC(fit),
    mape         = mean(abs(resid / df$Oxygen_norm), na.rm = TRUE)
  ) %>% dplyr::mutate(
    pass_REL_SE   = max_rel_SE  < REL_SE_THRESHOLD,
    pass_PVAL     = max_pval    < PVAL_THRESHOLD,
    pass_R2       = pseudo_R2  >= R2_THRESHOLD,
    pass_RESID    = resid_range < MAX_RESID_RANGE,
    pass_K_BOUNDS = is.finite(K) & K > 0 & K < 0.5,
    pass_R_BOUNDS = is.finite(r_per_minute) & dplyr::between(r_per_minute, 1e-4, 0.1),
    pass_AIC      = aic_gain   >= AIC_IMPROVEMENT,
    pass_MAPE     = mape        < MAPE_MAX,
    fit_ok = pass_REL_SE & pass_PVAL & pass_R2 & pass_RESID &
             pass_K_BOUNDS & pass_R_BOUNDS & pass_AIC & pass_MAPE,
    # headroom: how far the binding criterion is from its threshold, as a fraction
    headroom_REL_SE = (REL_SE_THRESHOLD - max_rel_SE) / REL_SE_THRESHOLD,
    headroom_R2     = (pseudo_R2 - R2_THRESHOLD) / (1 - R2_THRESHOLD),
    headroom_MAPE   = (MAPE_MAX - mape) / MAPE_MAX,
    headroom_RESID  = (MAX_RESID_RANGE - resid_range) / MAX_RESID_RANGE
  )
})
readr::write_csv(qc, file.path(OUT, "E3_qc_gate_per_series.csv"))

crit <- c("pass_REL_SE", "pass_PVAL", "pass_R2", "pass_RESID",
          "pass_K_BOUNDS", "pass_R_BOUNDS", "pass_AIC", "pass_MAPE")
n_fail_by_crit <- vapply(crit, function(k) as.integer(sum(!qc[[k]], na.rm = TRUE)), integer(1))
failing_by_crit <- vapply(crit, function(k) {
  f <- qc[!is.na(qc[[k]]) & !qc[[k]], , drop = FALSE]
  if (!nrow(f)) "" else paste(paste(f$Taxon, f$Replicate), collapse = "; ")
}, character(1))

qc_summary <- tibble(
  criterion = crit,
  threshold = c(sprintf("max rel SE < %.2f", REL_SE_THRESHOLD),
                sprintf("max p < %g", PVAL_THRESHOLD),
                sprintf("pseudo R2 >= %.2f", R2_THRESHOLD),
                sprintf("residual range < %.2f", MAX_RESID_RANGE),
                "0 < K < 0.5",
                "1e-4 <= r <= 0.1",
                sprintf("AIC gain >= %g", AIC_IMPROVEMENT),
                sprintf("MAPE < %.2f", MAPE_MAX)),
  n_series = nrow(qc),
  n_fail   = n_fail_by_crit,
  failing_series = failing_by_crit
)
readr::write_csv(qc_summary, file.path(OUT, "E3_qc_gate_by_criterion.csv"))

# ---- automatic gate vs manual curation --------------------------------------
excl <- if (file.exists(PLOT_EXCLUDE_CSV)) rd(PLOT_EXCLUDE_CSV) else tibble()
mfw  <- if (file.exists(MANUAL_FIT_WINDOWS_CSV)) rd(MANUAL_FIT_WINDOWS_CSV) else tibble()
auto_fail <- qc %>% dplyr::filter(!fit_ok) %>% dplyr::transmute(key = paste(Taxon, Replicate))
man_excl  <- if (nrow(excl)) tibble(key = paste(excl$Taxon, toupper(excl$Replicate))) else tibble(key = character(0))

gate <- tibble(
  metric = c("series in the run",
             "series that produced a fit",
             "series passing fit_ok (automatic gate)",
             "series failing fit_ok (automatic gate)",
             "series removed by plot_exclude_points.csv (manual)",
             "overlap: manually excluded AND automatically failed",
             "manually excluded but automatically PASSING",
             "automatically failed but NOT manually excluded",
             "curves with a manual fit window in manual_fit_windows.csv",
             "  of which both fit_start and fit_end are set",
             "  of which only one side is set"),
  value = c(nrow(qc),
            sum(qc$fitted),
            sum(qc$fit_ok, na.rm = TRUE),
            sum(!qc$fit_ok, na.rm = TRUE),
            nrow(man_excl),
            length(intersect(auto_fail$key, man_excl$key)),
            length(setdiff(man_excl$key, auto_fail$key)),
            length(setdiff(auto_fail$key, man_excl$key)),
            nrow(mfw),
            if (nrow(mfw)) sum(is.finite(mfw$fit_start) & is.finite(mfw$fit_end)) else 0L,
            if (nrow(mfw)) sum(xor(is.finite(mfw$fit_start), is.finite(mfw$fit_end))) else 0L)
)
readr::write_csv(gate, file.path(OUT, "E3_gate_vs_manual_curation.csv"))

# How much does the manual window actually move the fit window?
if (nrow(mfw)) {
  auto_meta <- rd(file.path(RUN, "Oxygen_Trimmed_Series_Metadata.csv")) %>%
    dplyr::transmute(Taxon, Replicate, auto_start = peak_time, auto_end = chosen_end_time)
  win_shift <- mfw %>%
    dplyr::mutate(Taxon = as.character(Taxon), Replicate = as.character(Replicate)) %>%
    dplyr::left_join(auto_meta, by = c("Taxon", "Replicate")) %>%
    dplyr::mutate(shift_start_min = fit_start - auto_start,
                  shift_end_min   = fit_end   - auto_end,
                  manual_window_min = fit_end - fit_start,
                  auto_window_min   = auto_end - auto_start,
                  window_ratio      = manual_window_min / auto_window_min)
  readr::write_csv(win_shift, file.path(OUT, "E3_manual_vs_auto_windows.csv"))
}

# ---------------------------------------------------------------------------
# E4. N0_BACKPROJECT: the distribution of delta and exp(r * delta)
# ---------------------------------------------------------------------------
res <- rd(file.path(RUN, "oxygen_results_with_R.csv")) %>%
  dplyr::mutate(amplification = exp(r_per_minute * delta_Ninoc_to_N0_min))
readr::write_csv(
  res %>% dplyr::select(Taxon, Replicate, r_per_minute, delta_Ninoc_to_N0_min,
                        N_inoc_cells_per_L, N0_cells_per_L, amplification,
                        R, R_C_fg_cell_h, G_C_fg_cell_h, fit_ok),
  file.path(OUT, "E4_backprojection_per_series.csv"))

bp <- res %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(
    n = dplyr::n(),
    delta_min_min = min(delta_Ninoc_to_N0_min), delta_med = stats::median(delta_Ninoc_to_N0_min),
    delta_max = max(delta_Ninoc_to_N0_min),
    amp_min = min(amplification), amp_med = stats::median(amplification), amp_max = max(amplification),
    r_med_per_min = stats::median(r_per_minute), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(amp_med))
readr::write_csv(bp, file.path(OUT, "E4_backprojection_by_taxon.csv"))

bp_overall <- tibble(
  statistic = c("n series", "delta min (min)", "delta median (min)", "delta max (min)",
                "exp(r*delta) min", "exp(r*delta) median", "exp(r*delta) max",
                "N0/N_inoc median", "fraction of series with exp(r*delta) > 2",
                "fraction with exp(r*delta) > 5"),
  value = c(nrow(res), min(res$delta_Ninoc_to_N0_min), stats::median(res$delta_Ninoc_to_N0_min),
            max(res$delta_Ninoc_to_N0_min), min(res$amplification),
            stats::median(res$amplification), max(res$amplification),
            stats::median(res$N0_cells_per_L / res$N_inoc_cells_per_L),
            mean(res$amplification > 2), mean(res$amplification > 5)))
readr::write_csv(bp_overall, file.path(OUT, "E4_backprojection_overall.csv"))

message("d1_audit: wrote E1-E4 artefacts to ", OUT)
