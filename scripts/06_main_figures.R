# =============================================================================
# 06_main_figures.R - publication figures (Fig 2, 3, 4, 5, 6, S1, Supp 3)
# =============================================================================
# Builds the manuscript figures from the fitted curves and reconstructed
# results.
#   * Fig 2   - 4-panel normalised O2 dynamics + fit, with r and R labels
#   * Supp 3  - normalised O2 dynamics faceted by taxon, coloured by replicate
#   * Supp S1 - fit residual diagnostics (residual vs time, residual vs fitted)
#   * Fig 3/4/5 - cross-method growth comparison (OD600 / FC vs O2)
#   * Fig 6   - per-cell respiration vs per-cell growth (RIS mixed model)
#
#   Input : results/tables/oxygen_fit_curves.csv          (from 05)
#           results/tables/oxygen_results_with_R.csv       (from 05)
#   Figures (results/figures/):
#     Fig_2_oxygen_dynamics_facet4_no_overflow.tiff
#     supp_Fig_3_oxygen_all_replicates_NORMALISED_by_replicate.tiff
#     Supp_Fig_S1_residuals.tiff
#     Fig_3/4/5 cross-method growth comparison
#     Fig_6_COMBINED.tiff (+ Fig_6a, Fig_6b) and Fig 6 data tables
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
  library(patchwork); library(grid); library(lme4); library(lmerTest)
  library(multcomp); library(ggsignif); library(gridExtra); library(glue); library(stringr)
})

# ---- inputs -----------------------------------------------------------------
fit_curves <- readr::read_csv(FITCURVES_CSV, show_col_types = FALSE)
results    <- readr::read_csv(RESULTS_FINAL_CSV, show_col_types = FALSE)

# =============================================================================
# Supp Fig S1 - residual diagnostics
# =============================================================================
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
    dplyr::mutate(Taxon = factor(Taxon), Replicate = factor(Replicate))

  p_resid_time <- ggplot2::ggplot(resid_all, ggplot2::aes(Time0, Residual, colour = Replicate)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_point(alpha = 0.7, size = 0.8) +
    ggplot2::facet_wrap(~ Taxon, scales = "free_x") +
    ggplot2::labs(x = "Time (min, re-zeroed)", y = "Residual (normalised O₂)") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(legend.position = "bottom",
                   strip.text = ggplot2::element_text(face = "italic"),
                   axis.title = ggplot2::element_text(size = 11),
                   axis.text  = ggplot2::element_text(size = 9))

  p_resid_fit <- ggplot2::ggplot(resid_all, ggplot2::aes(Fitted, Residual, colour = Replicate)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_point(alpha = 0.7, size = 0.8) +
    ggplot2::facet_wrap(~ Taxon, scales = "free_x") +
    ggplot2::labs(x = "Fitted normalised O₂", y = "Residual (normalised O₂)") +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(legend.position = "none",
                   strip.text = ggplot2::element_text(face = "italic"),
                   axis.title = ggplot2::element_text(size = 11),
                   axis.text  = ggplot2::element_text(size = 9))

  S1_residuals <- p_resid_time / p_resid_fit +
    patchwork::plot_layout(ncol = 1, heights = c(1, 1)) +
    patchwork::plot_annotation(
      title = "Supplementary Fig. S1. Residual diagnostics for O₂-based growth–respiration fits",
      tag_levels = "A"
    )

  save_tiff(S1_residuals, fig("Supp_Fig_S1_residuals.tiff"), width = 12, height = 10)
}

# =============================================================================
# Supp Fig 3 - normalised O2 dynamics, colour by replicate
# =============================================================================
all_fits_df_norm <- fit_curves %>%
  dplyr::mutate(TaxonFull = as.character(Taxon), Replicate = factor(Replicate)) %>%
  dplyr::select(TaxonFull, Replicate, Time0_min, Oxygen_norm, Fit_norm) %>%
  dplyr::rename(Time = Time0_min, Oxygen_n = Oxygen_norm, Pred_n = Fit_norm) %>%
  dplyr::arrange(TaxonFull, Replicate, Time)

supp_plot_norm_rep <- ggplot2::ggplot(all_fits_df_norm, ggplot2::aes(Time, Oxygen_n)) +
  ggplot2::geom_point(colour = "grey60", size = 1, alpha = 0.55) +
  ggplot2::geom_line(ggplot2::aes(y = Pred_n, colour = Replicate, group = Replicate),
                     linewidth = 0.9, na.rm = TRUE) +
  ggplot2::facet_wrap(~ TaxonFull, scales = "free_y") +
  ggplot2::labs(
    title = "Supplementary: Normalised O₂ Dynamics (colour = Replicate)",
    x = "Time since fit start (minutes)",
    y = expression("Normalised Oxygen (" * O[2] * "/" * O[2][0] * ")"),
    colour = "Replicate"
  ) +
  ggplot2::scale_colour_brewer(palette = "Dark2", drop = FALSE) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(legend.position = "bottom",
                 strip.text = ggplot2::element_text(face = "bold"),
                 axis.title = ggplot2::element_text(size = 12),
                 axis.text  = ggplot2::element_text(size = 10))

save_tiff(supp_plot_norm_rep,
          fig("supp_Fig_3_oxygen_all_replicates_NORMALISED_by_replicate.tiff"),
          width = 14, height = 10)

# =============================================================================
# Fig 2 - 4-panel facet with r and R labels
# =============================================================================
TEXT_SIZE  <- 3
X_ANCHOR   <- 0.00001
X_MARGIN   <- 0.05
Y_GAP_FRAC <- 0.80
Y_MARGIN   <- 0.05
ONE_LINE   <- FALSE

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

selected_combos <- FIG2_SELECTED_COMBOS
letters_vec <- letters[seq_len(nrow(selected_combos))]

facet_data <- fit_curves %>%
  dplyr::mutate(Time = Time0_min, Predicted_O2 = Fit_norm, Oxygen = Oxygen_norm) %>%
  dplyr::inner_join(selected_combos, by = c("Taxon", "Replicate")) %>%
  dplyr::select(Taxon, Replicate, Time, Oxygen, Predicted_O2)

if (nrow(facet_data) > 0) {

  annot_info <- results %>%
    dplyr::inner_join(selected_combos, by = c("Taxon", "Replicate")) %>%
    dplyr::mutate(
      FacetLabel = paste0(
        "(", letters_vec[match(paste(Taxon, Replicate),
                               paste(selected_combos$Taxon, selected_combos$Replicate))],
        ")~italic('", Taxon, "')"
      ),
      label_text = if (ONE_LINE) {
        paste0("italic(r)==", sprintf("%.2f", r_per_hour), "~h^{-1}*','~~",
               "italic(R)==", sci_pm(R, digits = 2), "~mg~O[2]~cell^{-1}~min^{-1}")
      } else {
        paste0("atop(",
               "italic(r)==", sprintf("%.2f", r_per_hour), "~h^{-1},",
               "italic(R)==", sci_pm(R, digits = 2), "~mg~O[2]~cell^{-1}~min^{-1})")
      }
    ) %>%
    dplyr::select(Taxon, Replicate, FacetLabel, label_text)

  facet_data <- facet_data %>%
    dplyr::inner_join(annot_info, by = c("Taxon", "Replicate"))

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

      tibble::tibble(x_pos = x_pos, y_pos = y_pos)
    }) %>%
    dplyr::ungroup()

  annotations <- facet_data %>%
    dplyr::distinct(FacetLabel, label_text) %>%
    dplyr::left_join(label_positions, by = "FacetLabel")

  facet_plot <- ggplot2::ggplot(facet_data, ggplot2::aes(x = Time)) +
    ggplot2::geom_point(ggplot2::aes(y = Oxygen), size = 2.2, alpha = 0.85) +
    ggplot2::geom_line(ggplot2::aes(y = Predicted_O2), linewidth = 1.2) +
    ggplot2::facet_wrap(~ FacetLabel, nrow = 1, labeller = ggplot2::label_parsed) +
    ggplot2::geom_text(data = annotations,
                       ggplot2::aes(x = x_pos, y = y_pos, label = label_text),
                       inherit.aes = FALSE, parse = TRUE, hjust = 0, vjust = 1.05,
                       size = TEXT_SIZE) +
    ggplot2::coord_cartesian(clip = "on") +
    ggplot2::labs(x = "Time since fit start (minutes)",
                  y = expression("Normalised Oxygen (" * O[2] / O[2 * "," * 0] * ")")) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(size = 16, face = "italic"),
      axis.title       = ggplot2::element_text(size = 17),
      axis.text        = ggplot2::element_text(size = 14),
      panel.spacing    = grid::unit(1.2, "lines"),
      plot.margin      = ggplot2::margin(12, 20, 12, 20),
      axis.line        = ggplot2::element_line(linewidth = 0.7),
      axis.ticks       = ggplot2::element_line(linewidth = 0.6)
    )

  save_tiff(facet_plot, fig("Fig_2_oxygen_dynamics_facet4_no_overflow.tiff"),
            width = 14, height = 4.5)
} else {
  message("Fig 2: none of FIG2_SELECTED_COMBOS matched the fitted curves; skipped.")
}

# =============================================================================
# Fig 3 / 4 / 5 - cross-method growth comparison (OD600 / FC vs O2)
# =============================================================================
FIG3_GROWTH_TIF    <- fig("Fig_3_growth_rate_comparison.tiff")
METHOD_EFFECTS_TIF <- fig("method_effects_CI.tiff")
FIG4_REGRESS_TIF   <- fig("Fig_4_growth_rate_regression_normalized_combined.tiff")
FIG5_BA_TIF        <- fig("Fig_5_BlandAltman_AllReplicates_LME.tiff")
OD_FC_CSV          <- file.path(data_dir, "OD_r_FC_r.csv")

results <- readr::read_csv(RESULTS_FINAL_CSV, show_col_types = FALSE)
results <- drop_excluded_taxa(results, "results")
results <- apply_excluded_curves(results, "results curves")

# =============================================================================
# Fig 3 - boxplot + method-effect mixed model
# =============================================================================
od_fc_data <- readr::read_csv(OD_FC_CSV, show_col_types = FALSE) %>%
  dplyr::mutate(Taxon = as.character(Taxon), Replicate = as.character(Replicate), Duration = Time)

OD_FC <- od_fc_data %>%
  dplyr::mutate(r_OD600 = (log(OD_Final) - log(OD_Initial)) / Duration,
                r_FC    = (log(FC_Final) - log(FC_Initial)) / Duration) %>%
  dplyr::select(Taxon, Replicate, Duration, r_OD600, r_FC)

growth <- results %>%
  dplyr::filter(fit_ok) %>%
  dplyr::select(Taxon, Replicate, r_per_minute) %>%
  dplyr::rename(r_O2 = r_per_minute) %>%
  dplyr::left_join(OD_FC, by = c("Taxon", "Replicate")) %>%
  dplyr::select(Taxon, Replicate, r_O2, r_OD600, r_FC) %>%
  dplyr::mutate(doubling_time_O2 = log(2) / r_O2,
                doubling_time_OD600 = log(2) / r_OD600,
                doubling_time_FC = log(2) / r_FC)
readr::write_csv(growth, tbl("growth_rates_combined.csv"))

cb_colors <- c("Oxygen" = "#E69F00", "OD600" = "#56B4E9", "Flow Cytometry" = "#009E73")
growth_long <- growth %>%
  tidyr::pivot_longer(cols = c(r_O2, r_OD600, r_FC), names_to = "Method", values_to = "Growth_Rate") %>%
  dplyr::mutate(Method = factor(Method, levels = c("r_O2", "r_OD600", "r_FC"),
                                labels = c("Oxygen", "OD600", "Flow Cytometry")))

growth_comparison <- ggplot2::ggplot(growth_long, ggplot2::aes(Taxon, Growth_Rate, fill = Method)) +
  ggplot2::geom_boxplot(outlier.size = 1.2, outlier.shape = 16) +
  ggplot2::scale_fill_manual(values = cb_colors) + ggplot2::theme_classic(base_size = 16) +
  ggplot2::labs(x = "Taxon", y = expression(paste("Growth rate (min"^{-1}, ")")), fill = "Method") +
  ggplot2::theme(axis.title = ggplot2::element_text(size = 18),
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 14, face = "italic"),
                 axis.text.y = ggplot2::element_text(size = 14),
                 legend.title = ggplot2::element_text(size = 16), legend.text = ggplot2::element_text(size = 14),
                 legend.position = "top", strip.text = ggplot2::element_text(size = 16, face = "italic")) +
  ggplot2::facet_grid(. ~ "By Taxon")

global_comparison <- ggplot2::ggplot(growth_long, ggplot2::aes(Method, Growth_Rate, fill = Method)) +
  ggplot2::geom_boxplot(outlier.size = 1.2, outlier.shape = 16) +
  ggsignif::geom_signif(comparisons = list(c("Oxygen", "OD600"), c("Oxygen", "Flow Cytometry"),
                                           c("OD600", "Flow Cytometry")),
                        annotations = c("ns", "ns", "ns"), step_increase = 0.1, tip_length = 0.01) +
  ggplot2::scale_fill_manual(values = cb_colors) + ggplot2::theme_classic(base_size = 16) +
  ggplot2::labs(x = "Method", y = expression(paste("Growth rate (min"^{-1}, ")"))) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = 18),
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 14),
                 axis.text.y = ggplot2::element_text(size = 14), legend.position = "none") +
  ggplot2::facet_grid(. ~ "All Taxa")

combined_comparison <- growth_comparison + global_comparison + patchwork::plot_layout(widths = c(3, 1))
save_tiff(combined_comparison, FIG3_GROWTH_TIF, width = 14, height = 6)

mm_method   <- lmer(Growth_Rate ~ 1 + Method + (1 | Taxon), REML = FALSE, data = growth_long)
mm_method_1 <- lmer(Growth_Rate ~ 1 + (1 | Taxon), REML = FALSE, data = growth_long)
ci_method    <- confint(mm_method, method = "Wald")
coefficients <- fixef(mm_method)
mc      <- multcomp::glht(mm_method, linfct = multcomp::mcp(Method = "Tukey"),
                          test = multcomp::adjusted("bonferroni"))
summary_mc <- summary(mc)

method_effects <- data.frame(
  Method = c("Oxygen", "OD600", "Flow Cytometry"),
  Estimate = c(coefficients[1], coefficients[1] + coefficients["MethodOD600"],
               coefficients[1] + coefficients["MethodFlow Cytometry"]),
  CI_lower = c(ci_method["(Intercept)", 1], ci_method["(Intercept)", 1] + ci_method["MethodOD600", 1],
               ci_method["(Intercept)", 1] + ci_method["MethodFlow Cytometry", 1]),
  CI_upper = c(ci_method["(Intercept)", 2], ci_method["(Intercept)", 2] + ci_method["MethodOD600", 2],
               ci_method["(Intercept)", 2] + ci_method["MethodFlow Cytometry", 2]))

pairs <- data.frame(y.position = max(method_effects$CI_upper, na.rm = TRUE) + c(0.0005, 0.001, 0.0015),
                    xmin = c(1, 1, 2), xmax = c(2, 3, 3), p.value = summary_mc$test$pvalues) %>%
  dplyr::mutate(stars = dplyr::case_when(p.value > 0.05 ~ "ns", p.value <= 0.05 & p.value > 0.01 ~ "*",
                                         p.value <= 0.01 & p.value > 0.001 ~ "**", p.value <= 0.001 ~ "***"))

method_comparison <- ggplot2::ggplot(method_effects, ggplot2::aes(Method, Estimate)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = CI_lower, ymax = CI_upper), width = 0.2) +
  ggplot2::geom_segment(data = pairs, ggplot2::aes(x = xmin, xend = xmax, y = y.position, yend = y.position)) +
  ggplot2::geom_text(data = pairs, ggplot2::aes(x = (xmin + xmax) / 2, y = y.position, label = stars),
                     vjust = -0.5, size = 5) +
  ggplot2::labs(x = "Method", y = expression("Growth rate (" * min^{-1} * ")"),
                title = "Comparison of Method Effects (Mixed Model)") +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_tiff(method_comparison, METHOD_EFFECTS_TIF, width = 8, height = 6)

readr::write_csv(method_effects, tbl("method_effects_estimates.csv"))
readr::write_csv(pairs,          tbl("method_effects_significance.csv"))
capture.output(summary(mm_OD <- lmer(r_O2 ~ r_OD600 + (1 | Taxon), data = growth)),
               file = tbl("mixed_model_OD600_summary.txt"))
capture.output(summary(mm_FC <- lmer(r_O2 ~ r_FC + (1 | Taxon), data = growth)),
               file = tbl("mixed_model_FC_summary.txt"))

# =============================================================================
# Fig 4 - LME regressions (r_O2 random-effect-normalised vs r_OD600 / r_FC)
# =============================================================================
create_norm_data <- function(model) {
  ranef_taxon <- ranef(model)$Taxon; fixed_coef <- fixef(model)
  rand_effects <- data.frame(Taxon = rownames(ranef_taxon), rand_int = ranef_taxon[, 1])
  growth_norm <- growth %>% dplyr::left_join(rand_effects, by = "Taxon") %>%
    dplyr::mutate(r_O2_norm = r_O2 - rand_int)
  r2_mixed <- cor(fitted(model), growth$r_O2, use = "complete.obs")^2
  list(data = growth_norm, fixed_coef = fixed_coef,
       eq_text = sprintf("y = %.3f + %.3fx\nR² = %.3f", fixed_coef[1], fixed_coef[2], r2_mixed))
}
od_norm <- create_norm_data(mm_OD); fc_norm <- create_norm_data(mm_FC)
rate_range <- range(c(od_norm$data$r_O2_norm, od_norm$data$r_OD600,
                      fc_norm$data$r_O2_norm, fc_norm$data$r_FC), na.rm = TRUE)
okabe_ito_extended <- TAXON_PALETTE

p_od <- ggplot2::ggplot(od_norm$data, ggplot2::aes(r_OD600, r_O2_norm)) +
  ggplot2::geom_point(ggplot2::aes(color = Taxon), size = 3) +
  ggplot2::geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  ggplot2::geom_abline(intercept = od_norm$fixed_coef[1], slope = od_norm$fixed_coef[2], color = "black") +
  ggplot2::scale_color_manual(values = okabe_ito_extended) + ggplot2::theme_classic() +
  ggplot2::labs(x = expression(paste("Growth Rate - OD"[600], " (min"^-1, ")")),
                y = expression(paste("Growth Rate - O"[2], " (min"^-1, ")"))) +
  ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank()) +
  ggplot2::coord_cartesian(xlim = rate_range, ylim = rate_range) +
  ggplot2::annotate("text", x = min(rate_range) + 0.1 * diff(rate_range),
                    y = max(rate_range) - 0.1 * diff(rate_range), label = od_norm$eq_text,
                    hjust = 0, vjust = 1, size = 5) + ggplot2::ggtitle("(A)")
p_fc <- ggplot2::ggplot(fc_norm$data, ggplot2::aes(r_FC, r_O2_norm)) +
  ggplot2::geom_point(ggplot2::aes(color = Taxon), size = 3) +
  ggplot2::geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  ggplot2::geom_abline(intercept = fc_norm$fixed_coef[1], slope = fc_norm$fixed_coef[2], color = "black") +
  ggplot2::scale_color_manual(values = okabe_ito_extended) + ggplot2::theme_classic() +
  ggplot2::labs(x = expression(paste("Growth Rate - FC (min"^-1, ")")),
                y = expression(paste("Growth Rate - O"[2], " (min"^-1, ")"))) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::coord_cartesian(xlim = rate_range, ylim = rate_range) +
  ggplot2::annotate("text", x = min(rate_range) + 0.1 * diff(rate_range),
                    y = max(rate_range) - 0.1 * diff(rate_range), label = fc_norm$eq_text,
                    hjust = 0, vjust = 1, size = 5) + ggplot2::ggtitle("(B)")
save_tiff(p_od / p_fc, FIG4_REGRESS_TIF, width = 10, height = 16)

# =============================================================================
# Fig 5 - Bland-Altman (O2 vs OD600, O2 vs FC)
# =============================================================================
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && is.finite(a)) a else b
df_raw <- readr::read_csv(tbl("growth_rates_combined.csv"), show_col_types = FALSE) %>%
  dplyr::rename(Oxygen_r = r_O2, OD_r = r_OD600, FC_r = r_FC) %>%
  dplyr::mutate(Taxon = as.factor(Taxon))
get_lmer_slope_p <- function(fit) {
  tab <- coef(summary(fit)); if (is.null(tab)) return(NA_real_)
  rn <- rownames(tab); if (is.null(rn)) return(NA_real_)
  idx <- which(rn == "avg"); if (!length(idx)) return(NA_real_)
  val <- tab[idx, "Pr(>|t|)"]; if (length(val) == 0) NA_real_ else as.numeric(val)
}
bland_altman_plot_lme <- function(df, method1, method2, label1, label2, panel_letter) {
  df_ba <- df %>% dplyr::mutate(method1 = {{ method1 }}, method2 = {{ method2 }},
                                avg = (method1 + method2) / 2, diff = method1 - method2) %>%
    dplyr::filter(is.finite(avg), is.finite(diff))
  bias <- mean(df_ba$diff, na.rm = TRUE); sd_diff <- sd(df_ba$diff, na.rm = TRUE)
  upper <- bias + 1.96 * sd_diff; lower <- bias - 1.96 * sd_diff
  within_limits <- mean(df_ba$diff >= lower & df_ba$diff <= upper, na.rm = TRUE) * 100
  df_reg <- df_ba %>% dplyr::filter(diff >= lower, diff <= upper)
  use_lmer <- dplyr::n_distinct(df_reg$Taxon) >= 2 && nrow(df_reg) >= 5
  if (use_lmer) {
    mm_fit <- lmer(diff ~ avg + (1 | Taxon), data = df_reg, REML = TRUE); fe <- fixef(mm_fit)
    intercept <- unname(fe[1]) %||% NA_real_; slope <- unname(fe[2]) %||% NA_real_
  } else {
    mm_fit <- lm(diff ~ avg, data = df_reg); coefs <- coef(mm_fit)
    intercept <- unname(coefs[1]) %||% NA_real_; slope <- unname(coefs[2]) %||% NA_real_
  }
  title_wrapped <- stringr::str_wrap(glue::glue("({panel_letter}) {label1} vs {label2} — Bland–Altman (all replicates)"), width = 55)
  ggplot2::ggplot(df_ba, ggplot2::aes(avg, diff)) +
    ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper, fill = "#D55E00", alpha = 0.15) +
    ggplot2::geom_point(size = 2.5, alpha = 0.9, color = "black") +
    ggplot2::geom_hline(yintercept = bias, color = "black", linewidth = 0.7) +
    ggplot2::geom_hline(yintercept = upper, color = "firebrick", linewidth = 0.8, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = lower, color = "firebrick", linewidth = 0.8, linetype = "dashed") +
    ggplot2::annotate("text", x = min(df_ba$avg, na.rm = TRUE), y = upper, label = sprintf("Upper LoA = %.3f", upper),
                      hjust = 0, vjust = -0.8, size = 4, fontface = "bold") +
    ggplot2::annotate("text", x = min(df_ba$avg, na.rm = TRUE), y = lower, label = sprintf("Lower LoA = %.3f", lower),
                      hjust = 0, vjust = 1.8, size = 4, fontface = "bold") +
    ggplot2::annotate("text", x = min(df_ba$avg, na.rm = TRUE), y = bias, label = sprintf("Bias = %.3f", bias),
                      hjust = 0, vjust = -2.6, size = 4, fontface = "bold") +
    ggplot2::annotate("text", x = min(df_ba$avg, na.rm = TRUE), y = lower - 0.005,
                      label = sprintf("%% within limits: %.1f%%", within_limits), hjust = 0, size = 4, fontface = "italic") +
    ggplot2::labs(title = title_wrapped, x = "Mean growth rate (per replicate)", y = "Difference in growth rate (per replicate)") +
    ggplot2::theme_classic(base_size = 14) + ggplot2::theme(plot.title = ggplot2::element_text(lineheight = 1.05))
}
plot1 <- bland_altman_plot_lme(df_raw, Oxygen_r, OD_r, "Oxygen", "OD600", "A")
plot2 <- bland_altman_plot_lme(df_raw, Oxygen_r, FC_r, "Oxygen", "Flow cytometry", "B")
save_tiff(gridExtra::arrangeGrob(plot1, plot2, ncol = 2), FIG5_BA_TIF, width = 14, height = 6)

message("06_main_figures: Fig 2, 3, 4, 5, S1, Supp 3 written.")

# =============================================================================
# Fig 6 - per-cell respiration vs per-cell growth (RIS mixed model)
# Folded in from the former standalone 07_fig6_growth_resp.R (restored from the
# original OxygenModel.R). Uses the same `results` table already loaded above.
# =============================================================================
plots_dir <- figures_dir

results <- readr::read_csv(file.path(tables_dir, "oxygen_results_with_R.csv"), show_col_types = FALSE)
if (!"fit_ok" %in% names(results)) results$fit_ok <- TRUE

FIG6_RIS_MAIN_TIF  <- file.path(plots_dir, "Fig_6a_r_vs_R_RIS_collapsed_MAIN.tiff")
FIG6B_SLOPES_TIF   <- file.path(plots_dir, "Fig_6b_taxon_specific_slopes.tiff")
FIG6_COMBINED_TIF  <- file.path(plots_dir, "Fig_6_COMBINED.tiff")

rm_df <- results %>%
  dplyr::filter(fit_ok, is.finite(G_C_fg_cell_h), G_C_fg_cell_h > 0,
                is.finite(R_C_fg_cell_h), R_C_fg_cell_h > 0) %>%
  dplyr::mutate(Taxon = factor(Taxon), log_G = log10(G_C_fg_cell_h), log_R = log10(R_C_fg_cell_h))

if (nrow(rm_df) < 8 || dplyr::n_distinct(rm_df$Taxon) < 2) {
  warning("Not enough data / taxa for Fig 6 random-slope RIS model.")
} else {
  x0 <- mean(rm_df$log_R, na.rm = TRUE)
  rm_df2 <- rm_df %>% dplyr::mutate(log_Rc = log_R - x0)

  mm <- suppressWarnings(lmerTest::lmer(log_G ~ log_Rc + (1 + log_Rc | Taxon), data = rm_df2, REML = FALSE))
  utils::capture.output(summary(mm), file = file.path(tables_dir, "mixed_model_Fig6_RIS_centeredX.txt"))
  fe <- lme4::fixef(mm)

  re_tbl <- lme4::ranef(mm)$Taxon %>% as.data.frame() %>% tibble::rownames_to_column("Taxon") %>%
    dplyr::transmute(Taxon = factor(Taxon, levels = levels(rm_df2$Taxon)),
                     b0 = .data[["(Intercept)"]], b1 = .data[["log_Rc"]])

  rm_collapsed <- rm_df2 %>% dplyr::left_join(re_tbl, by = "Taxon") %>%
    dplyr::mutate(log_G_collapsed = log_G - b0)

  x_seq   <- seq(min(rm_df2$log_R, na.rm = TRUE), max(rm_df2$log_R, na.rm = TRUE), length.out = 250)
  x_seq_c <- x_seq - x0
  V <- as.matrix(vcov(mm)); X <- cbind(1, x_seq_c)
  mu <- as.numeric(X %*% fe); se <- sqrt(pmax(0, rowSums((X %*% V) * X)))
  pred_global <- tibble::tibble(log_R = x_seq, mu = mu, lo = mu - 1.96 * se, hi = mu + 1.96 * se)

  EXPAND_FRAC <- 0.25; MIN_SPAN <- 0.35
  x_global_min <- min(rm_df2$log_R, na.rm = TRUE); x_global_max <- max(rm_df2$log_R, na.rm = TRUE)
  taxon_ranges <- rm_df2 %>% dplyr::group_by(Taxon) %>%
    dplyr::summarise(xmin = min(log_R, na.rm = TRUE), xmax = max(log_R, na.rm = TRUE),
                     xmid = mean(log_R, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(span = pmax(xmax - xmin, MIN_SPAN),
                  xmin2 = pmax(x_global_min, xmid - 0.5 * span), xmax2 = pmin(x_global_max, xmid + 0.5 * span),
                  xmin3 = pmax(x_global_min, xmin2 - EXPAND_FRAC * span), xmax3 = pmin(x_global_max, xmax2 + EXPAND_FRAC * span))
  taxon_lines <- re_tbl %>% dplyr::left_join(taxon_ranges, by = "Taxon") %>%
    dplyr::mutate(slope_tax = unname(fe[2]) + b1, intercept_collapsed = unname(fe[1])) %>%
    dplyr::rowwise() %>% dplyr::mutate(log_R_seq = list(seq(xmin3, xmax3, length.out = 120))) %>%
    dplyr::ungroup() %>% tidyr::unnest(log_R_seq) %>% dplyr::mutate(log_Rc_seq = log_R_seq - x0) %>%
    dplyr::transmute(Taxon, log_R = log_R_seq, mu_tax_collapsed = intercept_collapsed + slope_tax * log_Rc_seq)

  leg_nrow <- 3; n_tax <- nlevels(rm_df2$Taxon); leg_ncol <- ceiling(n_tax / leg_nrow)
  okabe <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#000000",
             "#332288","#88CCEE","#117733","#999933","#CC6677","#AA4499","#661100","#BBBBBB")
  pal <- if (n_tax <= length(okabe)) okabe[seq_len(n_tax)] else rep(okabe, length.out = n_tax)
  pal_named <- stats::setNames(pal, levels(rm_df2$Taxon))

  isme_theme_pub <- function() ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 15), axis.text = ggplot2::element_text(size = 13),
                   axis.line = ggplot2::element_line(linewidth = 0.9), axis.ticks = ggplot2::element_line(linewidth = 0.7),
                   legend.position = "top", legend.direction = "horizontal", legend.justification = "center",
                   legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_text(face = "italic", size = 11),
                   legend.key.width = grid::unit(1.0, "lines"), legend.key.height = grid::unit(0.75, "lines"),
                   legend.spacing.x = grid::unit(0.30, "lines"), legend.box = "vertical",
                   legend.margin = ggplot2::margin(b = 6, unit = "pt"), plot.margin = grid::unit(c(10,42,10,22), "pt"))

  r2_mixed <- cor(stats::fitted(mm), rm_df2$log_G, use = "complete.obs")^2
  m_fixed <- unname(fe[2]); n_fixed <- unname(fe[1] - fe[2] * x0)
  ann_text <- sprintf("Fixed effect: y = %.2f x + %.2f\nFixed slope = %.2f   R² (fitted vs observed) = %.2f",
                      m_fixed, n_fixed, m_fixed, r2_mixed)

  p6 <- ggplot2::ggplot(rm_collapsed, ggplot2::aes(x = log_R, y = log_G_collapsed, colour = Taxon)) +
    ggplot2::geom_ribbon(data = pred_global, ggplot2::aes(x = log_R, ymin = lo, ymax = hi), inherit.aes = FALSE, alpha = 0.06, colour = NA) +
    ggplot2::geom_line(data = pred_global, ggplot2::aes(x = log_R, y = mu), inherit.aes = FALSE, linewidth = 1.15, colour = "black") +
    ggplot2::geom_point(size = 2.2, alpha = 0.55) +
    ggplot2::geom_line(data = taxon_lines, ggplot2::aes(x = log_R, y = mu_tax_collapsed, colour = Taxon, group = Taxon), inherit.aes = FALSE, linewidth = 1.55, alpha = 0.98) +
    ggplot2::scale_color_manual(values = pal_named, labels = function(x) gsub("_", " ", x)) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = leg_nrow, ncol = leg_ncol, byrow = TRUE, override.aes = list(size = 3.1, alpha = 1, linewidth = 1.2))) +
    ggplot2::labs(x = expression(log[10](R~"(fg C cell"^{-1}~" h"^{-1}*")")),
                  y = expression("Taxon-collapsed " * log[10](G~"(fg C cell"^{-1}~" h"^{-1}*")"))) +
    ggplot2::annotate("text", x = min(rm_collapsed$log_R, na.rm = TRUE) + 0.03 * diff(range(rm_collapsed$log_R, na.rm = TRUE)),
                      y = max(rm_collapsed$log_G_collapsed, na.rm = TRUE) - 0.05 * diff(range(rm_collapsed$log_G_collapsed, na.rm = TRUE)),
                      label = ann_text, hjust = 0, vjust = 1, size = 4.1, lineheight = 1.05) +
    isme_theme_pub()
  ggplot2::ggsave(FIG6_RIS_MAIN_TIF, p6, width = 9.0, height = 6.9 + 0.55 * (leg_nrow - 1), dpi = 600, device = "tiff", compression = "lzw")
  readr::write_csv(rm_collapsed, file.path(tables_dir, "data_Fig6_RIS_COLLAPSED_centeredX.csv"))
  readr::write_csv(taxon_lines,  file.path(tables_dir, "Fig6_taxon_specific_lines_COLLAPSED_centeredX_EXTENDED.csv"))

  beta1 <- lme4::fixef(mm)[["log_Rc"]]; se_beta1 <- sqrt(as.matrix(vcov(mm))[["log_Rc", "log_Rc"]])
  re_tax <- lme4::ranef(mm, condVar = TRUE)$Taxon; pv <- attr(re_tax, "postVar")
  b1 <- re_tax[, "log_Rc"]; se_b1 <- sqrt(pmax(0, pv[2, 2, ]))
  slopes_df <- tibble::tibble(Taxon = factor(rownames(re_tax), levels = levels(rm_df2$Taxon)),
                              slope_tax = as.numeric(beta1 + b1), se_tax = sqrt(se_beta1^2 + se_b1^2),
                              lo = slope_tax - 1.96 * se_tax, hi = slope_tax + 1.96 * se_tax) %>%
    dplyr::filter(!is.na(Taxon)) %>% dplyr::arrange(slope_tax) %>%
    dplyr::mutate(Taxon_ord = factor(as.character(Taxon), levels = as.character(Taxon)))
  readr::write_csv(slopes_df, file.path(tables_dir, "Fig6b_taxon_specific_slopes.csv"))

  p6b <- ggplot2::ggplot(slopes_df, ggplot2::aes(x = slope_tax, y = Taxon_ord, colour = Taxon)) +
    ggplot2::geom_vline(xintercept = beta1, linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0.0, linewidth = 0.9) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::scale_color_manual(values = pal_named, guide = "none") +
    ggplot2::labs(x = "Taxon-specific slope", y = NULL) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic", size = 11), axis.title.x = ggplot2::element_text(size = 15))
  ggplot2::ggsave(FIG6B_SLOPES_TIF, p6b, width = 7.8, height = max(4.8, 0.28 * nrow(slopes_df) + 1.2), dpi = 600, device = "tiff", compression = "lzw")

  fig6_combined <- (p6 + p6b) + patchwork::plot_layout(widths = c(2.2, 1)) +
    patchwork::plot_annotation(tag_levels = list(c("(A)", "(B)")))
  ggplot2::ggsave(FIG6_COMBINED_TIF, fig6_combined, width = 16, height = 7.2, dpi = 600, device = "tiff", compression = "lzw")
  message("Figure 6 written: Fig_6_COMBINED.tiff (fixed slope ", round(m_fixed, 3), ", intercept ", round(n_fixed, 3), ")")
}
