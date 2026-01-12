###############################################################################
# MONTE CARLO SENSITIVITY (OUR MODEL) — MINIMAL OUTPUT VERSION
#
# Model:
#   O2_norm(t) = O2_0 + (K / r) * (1 - exp(r * t))
#
# Per-cell respiration:
#   R = K * O2_ref / N0
#
# Choice "B": each Taxon × Replicate uses its OWN fitted N0 as the mean (N0_mu).
#
# INPUTS:
#   - Tables/oxygen_results_with_R.csv
#       must contain: Taxon, Replicate, N0_cells_per_L, and K (or resp_tot/resp_rate)
#       may contain:  O2_ref (optional)
#
#   - Tables/Oxygen_Data_Filtered.csv  (ONLY if O2_ref missing/NA)
#
# OUTPUTS (ONLY 2 CSVs by default):
#   - Tables/N0_MC_ourmodel_taxon_summary_ALL.csv
#   - Tables/N0_MC_ourmodel_overall_summary_ALL.csv
#
# OPTIONAL (OFF by default):
#   - Tables/N0_MC_ourmodel_summary_per_series_ALL.csv
#   - Plots/N0_MC_ourmodel_R_rel_sd_density_ALL.(png/pdf)
#   - Plots/N0_MC_ourmodel_R_rel_sd_by_taxon_ALL.(png/pdf)
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(zoo)
  library(scales)
})

# ───────────────────────── 0. Project paths ──────────────────────────────── #
data_dir   <- "data"
plots_dir  <- "plots"
tables_dir <- "Tables"

dir.create(data_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# ───────────────────────── 0b. MC settings ───────────────────────────────── #
N_MC    <- 300
CV_GRID <- c(0.10, 0.20, 0.30, 0.40)

# Draw model:
#   "lognormal" (recommended) or "truncnorm"
DRAW_MODEL <- "lognormal"

# If using truncnorm, enforce a floor to prevent near-zero blowups:
TRUNC_FLOOR_FRAC <- 0.20   # N0_draw >= 0.20 * N0_mu

set.seed(123)

# ───────────────────────── 0c. OUTPUT CONTROLS ───────────────────────────── #
# (You said you generate too many CSVs -> defaults are MINIMAL)
SAVE_RAW_DRAWS         <- FALSE  # HUGE -> keep FALSE
SAVE_PER_SCENARIO_CSVS <- FALSE  # per-CV CSVs -> keep FALSE
SAVE_SERIES_ALL_CSV    <- FALSE  # 1 combined per-series file -> optional
SAVE_PLOTS             <- FALSE  # set TRUE if you want 2 combined plots total

# ───────────────────────── 1. Load results (K + N0) ─────────────────────── #
RESULTS_CSV <- file.path(tables_dir, "oxygen_results_with_R.csv")

model_results <- readr::read_csv(RESULTS_CSV, show_col_types = FALSE) %>%
  dplyr::mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate)
  )

# Map to K if needed (avoid masking errors by always using dplyr::select etc.)
if (!"K" %in% names(model_results)) {
  if ("resp_tot" %in% names(model_results)) {
    model_results <- model_results %>% dplyr::rename(K = resp_tot)
  } else if ("resp_rate" %in% names(model_results)) {
    model_results <- model_results %>% dplyr::rename(K = resp_rate)
  } else {
    stop("No column named 'K', 'resp_tot', or 'resp_rate' found in oxygen_results_with_R.csv")
  }
}

if (!"N0_cells_per_L" %in% names(model_results)) {
  stop("Column 'N0_cells_per_L' not found in oxygen_results_with_R.csv")
}

model_results <- model_results %>%
  dplyr::mutate(
    K              = as.numeric(K),
    N0_cells_per_L = as.numeric(N0_cells_per_L),
    O2_ref         = if ("O2_ref" %in% names(.)) as.numeric(O2_ref) else NA_real_
  )

# Optional: only good fits
if ("fit_ok" %in% names(model_results)) {
  model_results <- model_results %>% dplyr::filter(fit_ok %in% c(TRUE, "TRUE", 1))
}

model_results <- model_results %>%
  dplyr::filter(is.finite(K), is.finite(N0_cells_per_L), N0_cells_per_L > 0)

message("Rows after basic QC: ", nrow(model_results))

# ───────────────────────── 2. Ensure O2_ref exists (compute if needed) ───── #
need_o2ref <- (!("O2_ref" %in% names(model_results))) || any(!is.finite(model_results$O2_ref))

if (isTRUE(need_o2ref)) {
  
  message("O2_ref missing/NA → computing from Tables/Oxygen_Data_Filtered.csv")
  
  oxygen_filtered <- readr::read_csv(
    file.path(tables_dir, "Oxygen_Data_Filtered.csv"),
    show_col_types = FALSE
  ) %>%
    dplyr::mutate(
      Taxon     = as.character(Taxon),
      Replicate = as.character(Replicate),
      Time      = as.numeric(Time),
      Oxygen    = as.numeric(Oxygen)
    )
  
  # Compute O2_ref as mean of first 3 points AFTER a mild derivative onset detect
  O2_ref_series <- oxygen_filtered %>%
    dplyr::group_by(Taxon, Replicate) %>%
    dplyr::group_modify(~{
      df <- .x %>% dplyr::arrange(Time)
      if (nrow(df) < 5) return(tibble::tibble(O2_ref_calc = NA_real_))
      
      df <- df %>%
        dplyr::mutate(
          dO2   = c(NA_real_, diff(Oxygen)),
          dt    = c(NA_real_, diff(Time)),
          dO2dt = dO2 / dt,
          sm    = zoo::rollmean(dO2dt, 3, fill = NA_real_, align = "right")
        )
      
      idx <- which(df$sm < -1e-7)[1]
      idx <- if (is.na(idx) || idx > nrow(df) - 2) 1 else min(idx + 5, nrow(df) - 2)
      
      df_trim <- df[idx:nrow(df), , drop = FALSE]
      if (nrow(df_trim) < 3) return(tibble::tibble(O2_ref_calc = NA_real_))
      
      O0 <- mean(head(df_trim$Oxygen, 3), na.rm = TRUE)
      tibble::tibble(O2_ref_calc = as.numeric(O0))
    }) %>%
    dplyr::ungroup()
  
  model_results <- model_results %>%
    dplyr::left_join(O2_ref_series, by = c("Taxon", "Replicate")) %>%
    dplyr::mutate(
      O2_ref = dplyr::coalesce(O2_ref, O2_ref_calc),
      O2_ref = as.numeric(O2_ref)
    ) %>%
    dplyr::select(-dplyr::any_of("O2_ref_calc"))
}

model_results <- model_results %>%
  dplyr::filter(is.finite(O2_ref), O2_ref > 0)

message("Rows with usable O2_ref: ", nrow(model_results))

# ───────────────────────── 3. Vectorised helper: draw N0 ─────────────────── #
draw_N0 <- function(mu, cv, model = c("lognormal", "truncnorm"), floor_frac = 0.2) {
  model <- match.arg(model)
  
  mu <- as.numeric(mu)
  cv <- as.numeric(cv)
  
  n <- length(mu)
  if (length(cv) == 1) cv <- rep(cv, n)
  
  out <- rep(NA_real_, n)
  ok  <- is.finite(mu) & mu > 0 & is.finite(cv) & cv >= 0
  if (!any(ok)) return(out)
  
  if (model == "lognormal") {
    # mean = mu, sd = cv*mu
    sigma2  <- log1p(cv[ok]^2)
    sigma   <- sqrt(sigma2)
    meanlog <- log(mu[ok]) - 0.5 * sigma2
    out[ok] <- stats::rlnorm(sum(ok), meanlog = meanlog, sdlog = sigma)
    return(out)
  }
  
  # truncnorm (implemented via floor, no extra packages)
  sd <- cv[ok] * mu[ok]
  x  <- stats::rnorm(sum(ok), mean = mu[ok], sd = sd)
  floor_val <- floor_frac * mu[ok]
  out[ok] <- pmax(x, floor_val)
  out
}

# ───────────────────────── 4. Run scenarios ─────────────────────────────── #
all_taxon_summaries   <- list()
all_overall_summaries <- list()
all_series_summaries  <- list()  # optional 1-file output + plots

for (cv in CV_GRID) {
  
  tag <- paste0("cv", sprintf("%02d", round(100 * cv)), "_", DRAW_MODEL)
  message("Running scenario: ", tag)
  
  series_info <- model_results %>%
    dplyr::select(Taxon, Replicate, K, O2_ref, N0_cells_per_L) %>%
    dplyr::mutate(
      N0_mu = N0_cells_per_L,
      N0_cv = cv
    )
  
  mc_results <- series_info %>%
    tidyr::uncount(weights = N_MC, .id = "mc_id") %>%
    dplyr::mutate(
      N0_draw    = draw_N0(N0_mu, N0_cv, model = DRAW_MODEL, floor_frac = TRUNC_FLOOR_FRAC),
      R_draw     = K * O2_ref / N0_draw,
      N0_CV      = cv,
      draw_model = DRAW_MODEL
    ) %>%
    dplyr::filter(is.finite(N0_draw), N0_draw > 0, is.finite(R_draw), R_draw > 0)
  
  if (isTRUE(SAVE_RAW_DRAWS)) {
    readr::write_csv(
      mc_results,
      file.path(tables_dir, paste0("N0_MC_ourmodel_raw_R_draws_", tag, ".csv"))
    )
  }
  
  # per-series summary (small)
  mc_summary <- mc_results %>%
    dplyr::group_by(Taxon, Replicate, N0_CV, draw_model) %>%
    dplyr::summarise(
      n_MC       = dplyr::n(),
      N0_mu_used = dplyr::first(N0_mu),
      R_mean     = mean(R_draw, na.rm = TRUE),
      R_sd       = stats::sd(R_draw, na.rm = TRUE),
      R_q025     = stats::quantile(R_draw, 0.025, na.rm = TRUE),
      R_q975     = stats::quantile(R_draw, 0.975, na.rm = TRUE),
      R_rel_sd   = R_sd / R_mean,
      .groups    = "drop"
    ) %>%
    dplyr::mutate(N0_CV_pct = round(100 * N0_CV),
                  scenario  = tag)
  
  # taxon summary (small)
  mc_summary_taxon <- mc_summary %>%
    dplyr::group_by(Taxon, N0_CV, draw_model) %>%
    dplyr::summarise(
      n_reps       = dplyr::n(),
      R_rel_sd_med = stats::median(R_rel_sd, na.rm = TRUE),
      R_rel_sd_min = min(R_rel_sd, na.rm = TRUE),
      R_rel_sd_max = max(R_rel_sd, na.rm = TRUE),
      .groups      = "drop"
    ) %>%
    dplyr::mutate(N0_CV_pct = round(100 * N0_CV),
                  scenario  = tag)
  
  overall_stats <- mc_summary %>%
    dplyr::summarise(
      N0_CV        = cv,
      draw_model   = DRAW_MODEL,
      n_series     = dplyr::n(),
      R_rel_sd_med = stats::median(R_rel_sd, na.rm = TRUE),
      R_rel_sd_min = min(R_rel_sd, na.rm = TRUE),
      R_rel_sd_max = max(R_rel_sd, na.rm = TRUE)
    ) %>%
    dplyr::mutate(N0_CV_pct = round(100 * N0_CV),
                  scenario  = tag)
  
  # store for combined outputs
  all_taxon_summaries[[tag]]   <- mc_summary_taxon
  all_overall_summaries[[tag]] <- overall_stats
  all_series_summaries[[tag]]  <- mc_summary
  
  if (isTRUE(SAVE_PER_SCENARIO_CSVS)) {
    readr::write_csv(
      mc_summary,
      file.path(tables_dir, paste0("N0_MC_ourmodel_summary_per_series_", tag, ".csv"))
    )
    readr::write_csv(
      mc_summary_taxon,
      file.path(tables_dir, paste0("N0_MC_ourmodel_taxon_summary_", tag, ".csv"))
    )
  }
}

# ───────────────────────── 5. Write ONLY the combined outputs ───────────── #
taxon_summary_all <- dplyr::bind_rows(all_taxon_summaries) %>%
  dplyr::arrange(draw_model, N0_CV, Taxon)

readr::write_csv(
  taxon_summary_all,
  file.path(tables_dir, "N0_MC_ourmodel_taxon_summary_ALL.csv")
)

overall_summary_all <- dplyr::bind_rows(all_overall_summaries) %>%
  dplyr::arrange(draw_model, N0_CV)

readr::write_csv(
  overall_summary_all,
  file.path(tables_dir, "N0_MC_ourmodel_overall_summary_ALL.csv")
)

# Optional: ONE combined per-series file
series_summary_all <- dplyr::bind_rows(all_series_summaries) %>%
  dplyr::arrange(draw_model, N0_CV, Taxon, Replicate)

if (isTRUE(SAVE_SERIES_ALL_CSV)) {
  readr::write_csv(
    series_summary_all,
    file.path(tables_dir, "N0_MC_ourmodel_summary_per_series_ALL.csv")
  )
}

# ───────────────────────── 6. OPTIONAL: ONLY 2 PLOTS TOTAL ───────────────── #
if (isTRUE(SAVE_PLOTS)) {
  
  # Density by scenario (facet)
  p_rel_sd <- ggplot(series_summary_all, aes(x = R_rel_sd)) +
    geom_density(alpha = 0.35) +
    scale_x_continuous(labels = scientific) +
    facet_wrap(~ scenario, scales = "free_y") +
    labs(
      x = "Relative SD of R (SD/mean)",
      y = "Density",
      title = "MC sensitivity of R to N0 uncertainty (all scenarios)"
    ) +
    theme_classic(base_size = 12)
  
  ggsave(file.path(plots_dir, "N0_MC_ourmodel_R_rel_sd_density_ALL.png"),
         p_rel_sd, width = 10, height = 6, dpi = 300)
  ggsave(file.path(plots_dir, "N0_MC_ourmodel_R_rel_sd_density_ALL.pdf"),
         p_rel_sd, width = 10, height = 6)
  
  # Boxplot by taxon (facet by scenario)
  p_box <- ggplot(series_summary_all, aes(x = Taxon, y = R_rel_sd)) +
    geom_boxplot(outlier.size = 0.7, fill = "grey80") +
    scale_y_continuous(labels = scientific) +
    coord_flip() +
    facet_wrap(~ scenario, scales = "free_x") +
    labs(
      x = "Taxon",
      y = "Relative SD of R (SD/mean)",
      title = "Relative uncertainty in R by taxon (all scenarios)"
    ) +
    theme_classic(base_size = 11)
  
  ggsave(file.path(plots_dir, "N0_MC_ourmodel_R_rel_sd_by_taxon_ALL.png"),
         p_box, width = 12, height = 8, dpi = 300)
  ggsave(file.path(plots_dir, "N0_MC_ourmodel_R_rel_sd_by_taxon_ALL.pdf"),
         p_box, width = 12, height = 8)
}

message("DONE.")
message("Wrote:")
message(" - ", file.path(tables_dir, "N0_MC_ourmodel_taxon_summary_ALL.csv"))
message(" - ", file.path(tables_dir, "N0_MC_ourmodel_overall_summary_ALL.csv"))
if (isTRUE(SAVE_SERIES_ALL_CSV)) message(" - ", file.path(tables_dir, "N0_MC_ourmodel_summary_per_series_ALL.csv"))
if (isTRUE(SAVE_PLOTS)) message(" - 2 plots in: ", plots_dir)
###############################################################################
# End script
###############################################################################
