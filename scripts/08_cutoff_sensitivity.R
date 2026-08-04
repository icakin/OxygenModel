# =============================================================================
# 08_cutoff_sensitivity.R - Reviewer sensitivity test: O2_norm >= 0.5
# =============================================================================
# Re-fits the normalised O2 model twice per (Taxon, Replicate): once on the full
# trimmed window and once using only points with Oxygen_norm >= 0.5, then
# compares r and K between the two. Uses the trimmed/filtered data from 02
# (it does NOT re-trim, so it never overwrites the shared tables).
#
#   Input : results/tables/Oxygen_Data_Filtered.csv          (from 02)
#   Output: results/tables/oxygen_model_results_main.csv
#           results/tables/oxygen_model_results_O2_ge_0.5.csv
#           results/tables/Table_S3_oxygen_model_comparison_O2_ge_0.5.csv
#           results/tables/oxygen_model_comparison_O2_ge_0.5_summary.csv
#   Figures: results/figures/oxygen_dynamics_all_models.pdf
#            results/figures/oxygen_dynamics_fullsize_per_page.pdf
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
  library(minpack.lm)
  library(patchwork)
})

O2_THRESHOLD <- 0.5   # reviewer: use only Oxygen_norm >= 0.5

oxygen_data <- readr::read_csv(FILTERED_CSV, show_col_types = FALSE) %>%
  dplyr::mutate(Taxon = as.character(Taxon), Replicate = as.character(Replicate))
oxygen_data <- drop_excluded_taxa(oxygen_data, "filtered")
oxygen_data <- apply_excluded_curves(oxygen_data, "filtered curves")

empty_results <- function() tibble::tibble(
  Taxon = character(), Replicate = character(),
  r_per_minute = numeric(), r_per_hour = numeric(), K = numeric(),
  O2_0 = numeric(), O2_ref = numeric(), AICc = numeric(),
  lnO2_change_per_min = numeric(), pseudo_R2 = numeric(), fit_ok = logical()
)
results_main      <- empty_results()
results_O2_ge_0.5 <- empty_results()
plots_list <- list()

combos <- dplyr::group_keys(dplyr::group_by(oxygen_data, Taxon, Replicate))

for (i in seq_len(nrow(combos))) {
  Tax <- combos$Taxon[i]
  Rep <- combos$Replicate[i]

  df <- oxygen_data %>%
    dplyr::filter(Taxon == Tax, Replicate == Rep) %>%
    dplyr::arrange(Time)
  if (nrow(df) < 5) next

  df <- df %>%
    dplyr::mutate(
      dO2   = c(NA, diff(Oxygen)),
      dt    = c(NA, diff(Time)),
      dO2dt = dO2 / dt,
      sm    = zoo::rollmean(dO2dt, 3, fill = NA, align = "right")
    )
  idx <- which(df$sm < -1e-7)[1]
  idx <- if (is.na(idx) || idx > nrow(df) - 2) 1 else idx + 5
  df  <- df[idx:nrow(df), ] %>% dplyr::mutate(Time0 = Time - min(Time, na.rm = TRUE))
  if (nrow(df) < 5) next

  O0 <- mean(head(df$Oxygen, 3), na.rm = TRUE)
  df <- df %>% dplyr::mutate(Oxygen_norm = Oxygen / O0)

  r_start <- {
    seg    <- head(df, max(3, floor(.3 * nrow(df))))
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

  pred    <- stats::predict(fit, df)
  keep_ix <- abs(df$Oxygen_norm - pred) < 2 * stats::sd(df$Oxygen_norm - pred, na.rm = TRUE)
  df_kept <- df[keep_ix, ]
  if (nrow(df_kept) < 5) next
  fit <- tryCatch(stats::update(fit, data = df_kept), error = function(e) fit)

  pars  <- summary(fit)$coefficients
  r_est <- pars["r", "Estimate"]
  K_est <- pars["K", "Estimate"]
  pseudo_R2 <- 1 - sum(stats::residuals(fit)^2) /
    sum((df_kept$Oxygen_norm - mean(df_kept$Oxygen_norm))^2)

  fit_ok_main <- all(
    abs(pars[, "Std. Error"] / pars[, "Estimate"]) < REL_SE_THRESHOLD,
    pars[, "Pr(>|t|)"] < PVAL_THRESHOLD,
    pseudo_R2 >= R2_THRESHOLD,
    diff(range(stats::residuals(fit), na.rm = TRUE)) < MAX_RESID_RANGE,
    is.finite(K_est), K_est > 0, K_est < 0.5,
    dplyr::between(r_est, 1e-4, 0.1),
    stats::AIC(stats::lm(Oxygen_norm ~ 1, data = df_kept)) - stats::AIC(fit) >= AIC_IMPROVEMENT,
    mean(abs(stats::residuals(fit) / df_kept$Oxygen_norm)) < MAPE_MAX
  )

  results_main <- results_main %>%
    dplyr::add_row(
      Taxon = Tax, Replicate = Rep,
      r_per_minute = r_est, r_per_hour = r_est * 60, K = K_est,
      O2_0 = coef(fit)["O2_0"], O2_ref = O0, AICc = stats::AIC(fit),
      lnO2_change_per_min = (log(dplyr::last(df_kept$Oxygen_norm)) -
                               log(dplyr::first(df_kept$Oxygen_norm))) /
        (dplyr::last(df_kept$Time0) - dplyr::first(df_kept$Time0)),
      pseudo_R2 = pseudo_R2, fit_ok = fit_ok_main
    )

  df_kept <- df_kept %>% dplyr::mutate(Pred = stats::predict(fit, df_kept))
  plot_key <- paste(Tax, Rep, sep = "_")
  plots_list[[plot_key]] <-
    ggplot2::ggplot(df_kept, ggplot2::aes(Time0, Oxygen_norm)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_line(ggplot2::aes(y = Pred), linewidth = .8, colour = "red") +
    ggplot2::labs(title = plot_key, x = expression(Time ~ "(min)"),
                  y = expression(Normalised ~ O[2])) +
    isme_theme()

  # ---- fit using only Oxygen_norm >= 0.5 ------------------------------------
  df_05 <- df %>% dplyr::filter(Oxygen_norm >= O2_THRESHOLD)
  if (nrow(df_05) >= 5) {
    r_start_05 <- {
      seg    <- head(df_05, max(3, floor(.3 * nrow(df_05))))
      slopes <- abs(diff(log(pmax(seg$Oxygen_norm, 1e-6))) / diff(seg$Time0))
      pmin(pmax(max(slopes, na.rm = TRUE), 1e-4), 5e-2)
    }
    K_start_05 <- {
      slope <- suppressWarnings(min(diff(df_05$Oxygen_norm) / diff(df_05$Time0), na.rm = TRUE))
      if (!is.finite(slope)) slope <- -1e-3
      pmin(pmax(abs(slope), 1e-5), 0.1)
    }
    fit_05 <- tryCatch(
      minpack.lm::nlsLM(Oxygen_norm ~ resp_model(r, K, Time0, O2_0), data = df_05,
                        start = list(r = r_start_05, K = K_start_05, O2_0 = 1),
                        lower = FIT_LOWER, upper = FIT_UPPER,
                        control = minpack.lm::nls.lm.control(maxiter = 300)),
      error = function(e) NULL)

    if (!is.null(fit_05)) {
      pred_05    <- stats::predict(fit_05, df_05)
      keep_ix_05 <- abs(df_05$Oxygen_norm - pred_05) < 2 * stats::sd(df_05$Oxygen_norm - pred_05, na.rm = TRUE)
      df_05_kept <- df_05[keep_ix_05, ]

      if (nrow(df_05_kept) >= 5) {
        fit_05 <- tryCatch(stats::update(fit_05, data = df_05_kept), error = function(e) fit_05)
        pars_05  <- summary(fit_05)$coefficients
        r_est_05 <- pars_05["r", "Estimate"]
        K_est_05 <- pars_05["K", "Estimate"]
        pseudo_R2_05 <- 1 - sum(stats::residuals(fit_05)^2) /
          sum((df_05_kept$Oxygen_norm - mean(df_05_kept$Oxygen_norm))^2)

        fit_ok_05 <- all(
          abs(pars_05[, "Std. Error"] / pars_05[, "Estimate"]) < REL_SE_THRESHOLD,
          pars_05[, "Pr(>|t|)"] < PVAL_THRESHOLD,
          pseudo_R2_05 >= R2_THRESHOLD,
          diff(range(stats::residuals(fit_05), na.rm = TRUE)) < MAX_RESID_RANGE,
          is.finite(K_est_05), K_est_05 > 0, K_est_05 < 0.5,
          dplyr::between(r_est_05, 1e-4, 0.1),
          stats::AIC(stats::lm(Oxygen_norm ~ 1, data = df_05_kept)) - stats::AIC(fit_05) >= AIC_IMPROVEMENT,
          mean(abs(stats::residuals(fit_05) / df_05_kept$Oxygen_norm)) < MAPE_MAX
        )

        results_O2_ge_0.5 <- results_O2_ge_0.5 %>%
          dplyr::add_row(
            Taxon = Tax, Replicate = Rep,
            r_per_minute = r_est_05, r_per_hour = r_est_05 * 60, K = K_est_05,
            O2_0 = coef(fit_05)["O2_0"], O2_ref = O0, AICc = stats::AIC(fit_05),
            lnO2_change_per_min = (log(dplyr::last(df_05_kept$Oxygen_norm)) -
                                     log(dplyr::first(df_05_kept$Oxygen_norm))) /
              (dplyr::last(df_05_kept$Time0) - dplyr::first(df_05_kept$Time0)),
            pseudo_R2 = pseudo_R2_05, fit_ok = fit_ok_05
          )
      }
    }
  }
}

readr::write_csv(results_main,      tbl("oxygen_model_results_main.csv"))
readr::write_csv(results_O2_ge_0.5, tbl("oxygen_model_results_O2_ge_0.5.csv"))

if (length(plots_list) > 0) {
  grDevices::pdf(fig("oxygen_dynamics_all_models.pdf"), width = 14, height = 10)
  print(patchwork::wrap_plots(plots_list))
  grDevices::dev.off()

  grDevices::pdf(fig("oxygen_dynamics_fullsize_per_page.pdf"), width = 8, height = 6)
  purrr::walk(plots_list, print)
  grDevices::dev.off()
}

# ---- comparison: full trimmed vs O2 >= 0.5 ----------------------------------
comparison_05 <- results_main %>%
  dplyr::filter(fit_ok) %>%
  dplyr::select(Taxon, Replicate, r_main = r_per_minute, K_main = K) %>%
  dplyr::inner_join(
    results_O2_ge_0.5 %>%
      dplyr::filter(fit_ok) %>%
      dplyr::select(Taxon, Replicate, r_05 = r_per_minute, K_05 = K),
    by = c("Taxon", "Replicate")
  ) %>%
  dplyr::mutate(r_ratio = r_05 / r_main, r_diff = r_05 - r_main,
                K_ratio = K_05 / K_main, K_diff = K_05 - K_main)

readr::write_csv(comparison_05, tbl("Table_S3_oxygen_model_comparison_O2_ge_0.5.csv"))

if (nrow(comparison_05) > 0) {
  comparison_05_summary <- comparison_05 %>%
    dplyr::summarise(
      n_paired_fits          = dplyr::n(),
      cor_r_main_r_05        = cor(r_main, r_05, use = "complete.obs"),
      median_r_ratio         = median(r_ratio, na.rm = TRUE),
      q25_r_ratio            = quantile(r_ratio, 0.25, na.rm = TRUE),
      q75_r_ratio            = quantile(r_ratio, 0.75, na.rm = TRUE),
      max_abs_r_ratio_minus1 = max(abs(r_ratio - 1), na.rm = TRUE),
      cor_K_main_K_05        = cor(K_main, K_05, use = "complete.obs"),
      median_K_ratio         = median(K_ratio, na.rm = TRUE),
      q25_K_ratio            = quantile(K_ratio, 0.25, na.rm = TRUE),
      q75_K_ratio            = quantile(K_ratio, 0.75, na.rm = TRUE),
      max_abs_K_ratio_minus1 = max(abs(K_ratio - 1), na.rm = TRUE)
    )
  readr::write_csv(comparison_05_summary, tbl("oxygen_model_comparison_O2_ge_0.5_summary.csv"))
}

message("07_cutoff_sensitivity: ", nrow(comparison_05), " paired fits compared.")
