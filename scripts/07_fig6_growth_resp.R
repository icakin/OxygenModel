# =============================================================================
# 07_fig6_growth_resp.R - Figure 6: per-cell respiration vs per-cell growth
# =============================================================================
# Restored from the original publication code (OxygenModel.R). Fits the
# random-intercept-and-slope mixed model log10(G) ~ log10(R) with centered x,
# taxon-collapses the y-axis (removes taxon random intercepts) for display, and
# builds the combined A/B figure:
#   Panel A - taxon-collapsed log10(G) vs log10(R), taxon + fixed-effect lines
#   Panel B - taxon-specific slope estimates with 95% intervals
#   INPUT   results/tables/oxygen_results_with_R.csv
#   OUTPUT  results/figures/Fig_6_COMBINED.tiff (+ 6a, 6b) and data tables
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
  library(dplyr); library(tidyr); library(tibble); library(readr)
  library(ggplot2); library(patchwork); library(lme4); library(lmerTest); library(grid)
})
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
