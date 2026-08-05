# =============================================================================
# 09_montecarlo_N0.R - Monte-Carlo sensitivity of respiration R to N0
# =============================================================================
# For each (Taxon, Replicate), treats the fitted N0 as the mean and samples N0
# across a grid of CV values, propagating the uncertainty into R = K*O2_ref/N0
# and summarising the relative spread. Reads the reconstructed results from 03.
#
#   Input : results/tables/oxygen_results_with_R.csv         (from 05)
#           results/tables/Oxygen_Data_Filtered.csv          (only if O2_ref missing)
#   Output: results/tables/N0_MC_ourmodel_taxon_summary_ALL.csv
#           results/tables/N0_MC_ourmodel_overall_summary_ALL.csv
#   Figures (SAVE_PLOTS = TRUE): results/figures/N0_MC_ourmodel_R_rel_sd_density_ALL.(png|pdf)
#                                results/figures/N0_MC_ourmodel_R_rel_sd_by_taxon_ALL.(png|pdf)
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
  library(zoo)
  library(scales)
})

# ---- MC settings ------------------------------------------------------------
N_MC             <- 300
CV_GRID          <- c(0.10, 0.20, 0.30, 0.40)
DRAW_MODEL       <- "lognormal"   # "lognormal" or "truncnorm"
TRUNC_FLOOR_FRAC <- 0.20
set.seed(123)

# Output controls (minimal by default; plots ON so the sensitivity figures exist)
SAVE_RAW_DRAWS         <- FALSE
SAVE_PER_SCENARIO_CSVS <- FALSE
SAVE_SERIES_ALL_CSV    <- FALSE
SAVE_PLOTS             <- TRUE

# ---- load results (K + N0) --------------------------------------------------
model_results <- readr::read_csv(RESULTS_FINAL_CSV, show_col_types = FALSE) %>%
  dplyr::mutate(Taxon = as.character(Taxon), Replicate = as.character(Replicate))

if (!"K" %in% names(model_results)) {
  if ("resp_tot" %in% names(model_results)) {
    model_results <- model_results %>% dplyr::rename(K = resp_tot)
  } else if ("resp_rate" %in% names(model_results)) {
    model_results <- model_results %>% dplyr::rename(K = resp_rate)
  } else {
    stop("No column named 'K', 'resp_tot', or 'resp_rate' in oxygen_results_with_R.csv")
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
if ("fit_ok" %in% names(model_results)) {
  model_results <- model_results %>% dplyr::filter(fit_ok %in% c(TRUE, "TRUE", 1))
}
model_results <- model_results %>%
  dplyr::filter(is.finite(K), is.finite(N0_cells_per_L), N0_cells_per_L > 0)
message("Rows after basic QC: ", nrow(model_results))

# ---- ensure O2_ref exists ---------------------------------------------------
need_o2ref <- (!("O2_ref" %in% names(model_results))) || any(!is.finite(model_results$O2_ref))
if (isTRUE(need_o2ref)) {
  message("O2_ref missing/NA -> computing from Oxygen_Data_Filtered.csv")
  oxygen_filtered <- readr::read_csv(FILTERED_CSV, show_col_types = FALSE) %>%
    dplyr::mutate(Taxon = as.character(Taxon), Replicate = as.character(Replicate),
                  Time = as.numeric(Time), Oxygen = as.numeric(Oxygen))

  O2_ref_series <- oxygen_filtered %>%
    dplyr::group_by(Taxon, Replicate) %>%
    dplyr::group_modify(~{
      df <- .x %>% dplyr::arrange(Time)
      if (nrow(df) < 5) return(tibble::tibble(O2_ref_calc = NA_real_))
      df <- df %>% dplyr::mutate(
        dO2 = c(NA_real_, diff(Oxygen)), dt = c(NA_real_, diff(Time)),
        dO2dt = dO2 / dt, sm = zoo::rollmean(dO2dt, 3, fill = NA_real_, align = "right"))
      idx <- which(df$sm < -1e-7)[1]
      idx <- if (is.na(idx) || idx > nrow(df) - 2) 1 else min(idx + 5, nrow(df) - 2)
      df_trim <- df[idx:nrow(df), , drop = FALSE]
      if (nrow(df_trim) < 3) return(tibble::tibble(O2_ref_calc = NA_real_))
      tibble::tibble(O2_ref_calc = as.numeric(mean(head(df_trim$Oxygen, 3), na.rm = TRUE)))
    }) %>%
    dplyr::ungroup()

  model_results <- model_results %>%
    dplyr::left_join(O2_ref_series, by = c("Taxon", "Replicate")) %>%
    dplyr::mutate(O2_ref = dplyr::coalesce(O2_ref, O2_ref_calc),
                  O2_ref = as.numeric(O2_ref)) %>%
    dplyr::select(-dplyr::any_of("O2_ref_calc"))
}
model_results <- model_results %>% dplyr::filter(is.finite(O2_ref), O2_ref > 0)
message("Rows with usable O2_ref: ", nrow(model_results))

# ---- vectorised N0 draw -----------------------------------------------------
draw_N0 <- function(mu, cv, model = c("lognormal", "truncnorm"), floor_frac = 0.2) {
  model <- match.arg(model)
  mu <- as.numeric(mu); cv <- as.numeric(cv)
  n <- length(mu); if (length(cv) == 1) cv <- rep(cv, n)
  out <- rep(NA_real_, n)
  ok  <- is.finite(mu) & mu > 0 & is.finite(cv) & cv >= 0
  if (!any(ok)) return(out)
  if (model == "lognormal") {
    sigma2  <- log1p(cv[ok]^2)
    sigma   <- sqrt(sigma2)
    meanlog <- log(mu[ok]) - 0.5 * sigma2
    out[ok] <- stats::rlnorm(sum(ok), meanlog = meanlog, sdlog = sigma)
    return(out)
  }
  sd <- cv[ok] * mu[ok]
  x  <- stats::rnorm(sum(ok), mean = mu[ok], sd = sd)
  out[ok] <- pmax(x, floor_frac * mu[ok])
  out
}

# ---- run scenarios ----------------------------------------------------------
all_taxon_summaries   <- list()
all_overall_summaries <- list()
all_series_summaries  <- list()

for (cv in CV_GRID) {
  tag <- paste0("cv", sprintf("%02d", round(100 * cv)), "_", DRAW_MODEL)
  message("Running scenario: ", tag)

  series_info <- model_results %>%
    dplyr::select(Taxon, Replicate, K, O2_ref, N0_cells_per_L) %>%
    dplyr::mutate(N0_mu = N0_cells_per_L, N0_cv = cv)

  mc_results <- series_info %>%
    tidyr::uncount(weights = N_MC, .id = "mc_id") %>%
    dplyr::mutate(
      N0_draw = draw_N0(N0_mu, N0_cv, model = DRAW_MODEL, floor_frac = TRUNC_FLOOR_FRAC),
      R_draw  = K * O2_ref / N0_draw, N0_CV = cv, draw_model = DRAW_MODEL
    ) %>%
    dplyr::filter(is.finite(N0_draw), N0_draw > 0, is.finite(R_draw), R_draw > 0)

  if (isTRUE(SAVE_RAW_DRAWS))
    readr::write_csv(mc_results, tbl(paste0("N0_MC_ourmodel_raw_R_draws_", tag, ".csv")))

  mc_summary <- mc_results %>%
    dplyr::group_by(Taxon, Replicate, N0_CV, draw_model) %>%
    dplyr::summarise(
      n_MC = dplyr::n(), N0_mu_used = dplyr::first(N0_mu),
      R_mean = mean(R_draw, na.rm = TRUE), R_sd = stats::sd(R_draw, na.rm = TRUE),
      R_q025 = stats::quantile(R_draw, 0.025, na.rm = TRUE),
      R_q975 = stats::quantile(R_draw, 0.975, na.rm = TRUE),
      R_rel_sd = stats::sd(R_draw, na.rm = TRUE) / mean(R_draw, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(N0_CV_pct = round(100 * N0_CV), scenario = tag)

  mc_summary_taxon <- mc_summary %>%
    dplyr::group_by(Taxon, N0_CV, draw_model) %>%
    dplyr::summarise(n_reps = dplyr::n(),
                     R_rel_sd_med = stats::median(R_rel_sd, na.rm = TRUE),
                     R_rel_sd_min = min(R_rel_sd, na.rm = TRUE),
                     R_rel_sd_max = max(R_rel_sd, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(N0_CV_pct = round(100 * N0_CV), scenario = tag)

  overall_stats <- mc_summary %>%
    dplyr::summarise(N0_CV = cv, draw_model = DRAW_MODEL, n_series = dplyr::n(),
                     R_rel_sd_med = stats::median(R_rel_sd, na.rm = TRUE),
                     R_rel_sd_min = min(R_rel_sd, na.rm = TRUE),
                     R_rel_sd_max = max(R_rel_sd, na.rm = TRUE)) %>%
    dplyr::mutate(N0_CV_pct = round(100 * N0_CV), scenario = tag)

  all_taxon_summaries[[tag]]   <- mc_summary_taxon
  all_overall_summaries[[tag]] <- overall_stats
  all_series_summaries[[tag]]  <- mc_summary

  if (isTRUE(SAVE_PER_SCENARIO_CSVS)) {
    readr::write_csv(mc_summary,       tbl(paste0("N0_MC_ourmodel_summary_per_series_", tag, ".csv")))
    readr::write_csv(mc_summary_taxon, tbl(paste0("N0_MC_ourmodel_taxon_summary_", tag, ".csv")))
  }
}

taxon_summary_all <- dplyr::bind_rows(all_taxon_summaries) %>%
  dplyr::arrange(draw_model, N0_CV, Taxon)
readr::write_csv(taxon_summary_all, tbl("N0_MC_ourmodel_taxon_summary_ALL.csv"))

overall_summary_all <- dplyr::bind_rows(all_overall_summaries) %>%
  dplyr::arrange(draw_model, N0_CV)
readr::write_csv(overall_summary_all, tbl("N0_MC_ourmodel_overall_summary_ALL.csv"))

series_summary_all <- dplyr::bind_rows(all_series_summaries) %>%
  dplyr::arrange(draw_model, N0_CV, Taxon, Replicate)
if (isTRUE(SAVE_SERIES_ALL_CSV))
  readr::write_csv(series_summary_all, tbl("N0_MC_ourmodel_summary_per_series_ALL.csv"))

# ---- figures ----------------------------------------------------------------
if (isTRUE(SAVE_PLOTS) && nrow(series_summary_all) > 0) {
  p_rel_sd <- ggplot2::ggplot(series_summary_all, ggplot2::aes(x = R_rel_sd)) +
    ggplot2::geom_density(alpha = 0.35) +
    ggplot2::scale_x_continuous(labels = scales::scientific) +
    ggplot2::facet_wrap(~ scenario, scales = "free_y") +
    ggplot2::labs(x = "Relative SD of R (SD/mean)", y = "Density",
                  title = "MC sensitivity of R to N0 uncertainty (all scenarios)") +
    ggplot2::theme_classic(base_size = 12)
  ggplot2::ggsave(fig("N0_MC_ourmodel_R_rel_sd_density_ALL.png"), p_rel_sd, width = 10, height = 6, dpi = 300)
  ggplot2::ggsave(fig("N0_MC_ourmodel_R_rel_sd_density_ALL.pdf"), p_rel_sd, width = 10, height = 6)

  p_box <- ggplot2::ggplot(series_summary_all, ggplot2::aes(x = Taxon, y = R_rel_sd)) +
    ggplot2::geom_boxplot(outlier.size = 0.7, fill = "grey80") +
    ggplot2::scale_y_continuous(labels = scales::scientific) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ scenario, scales = "free_x") +
    ggplot2::labs(x = "Taxon", y = "Relative SD of R (SD/mean)",
                  title = "Relative uncertainty in R by taxon (all scenarios)") +
    ggplot2::theme_classic(base_size = 11)
  ggplot2::ggsave(fig("N0_MC_ourmodel_R_rel_sd_by_taxon_ALL.png"), p_box, width = 12, height = 8, dpi = 300)
  ggplot2::ggsave(fig("N0_MC_ourmodel_R_rel_sd_by_taxon_ALL.pdf"), p_box, width = 12, height = 8)
}

message("09_montecarlo_N0: done.")
