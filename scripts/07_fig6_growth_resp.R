# =============================================================================
# 07_fig6_growth_resp.R - Figure 6: per-cell respiration vs per-cell growth
#                         (the r-vs-R relationship), reproducible from pipeline.
# =============================================================================
# Rebuilds Figure 6 and its data directly from oxygen_results_with_R.csv:
#   Panel A - log10 R vs log10 G, per-taxon regression lines + overall line
#   Panel B - per-taxon slope estimates with 95% confidence intervals
# so the figure has a generating script and a traceable source.
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
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork) })

res <- readr::read_csv(file.path(tables_dir, "oxygen_results_with_R.csv"), show_col_types = FALSE) |>
  dplyr::filter(is.finite(R_C_fg_cell_h), R_C_fg_cell_h > 0,
                is.finite(G_C_fg_cell_h), G_C_fg_cell_h > 0) |>
  dplyr::mutate(log_R = log10(R_C_fg_cell_h), log_G = log10(G_C_fg_cell_h))
readr::write_csv(res |> dplyr::select(Taxon, Replicate, log_R, log_G),
                 file.path(tables_dir, "data_Fig6_growth_vs_resp.csv"))

# overall slope: mixed-model fixed effect (matches the paper's approach), with a
# pooled-lm fallback if the RIS fails to fit.
fe <- tryCatch({
  suppressWarnings({
    ris <- lme4::lmer(log_R ~ log_G + (1 + log_G | Taxon), data = res,
                      control = lme4::lmerControl(optimizer = "bobyqa"))
    lme4::fixef(ris)
  })
}, error = function(e) coef(lm(log_R ~ log_G, data = res)))
fe_int <- as.numeric(fe[1]); fe_slp <- as.numeric(fe[2])

# per-taxon regression slopes + 95% CI (each taxon's own lm)
bytax <- res |> dplyr::group_by(Taxon) |> dplyr::group_modify(~{
  m <- lm(log_R ~ log_G, data = .x); ci <- confint(m)["log_G", ]
  tibble::tibble(int = coef(m)[1], slope = coef(m)[2], lo = ci[1], hi = ci[2])
}) |> dplyr::ungroup()
readr::write_csv(bytax, file.path(tables_dir, "Fig6_taxon_slopes.csv"))

# ---- Panel A ----------------------------------------------------------------
pA <- ggplot(res, aes(log_G, log_R, colour = Taxon)) +
  geom_point(size = 1.6, alpha = .8) +
  geom_abline(data = bytax, aes(intercept = int, slope = slope, colour = Taxon), linewidth = .4) +
  geom_abline(intercept = fe_int, slope = fe_slp, linewidth = 1.1, colour = "black") +
  annotate("text", x = min(res$log_G), y = max(res$log_R),
           label = sprintf("overall slope = %.2f", fe_slp), hjust = 0, size = 3.8) +
  labs(title = "A", x = expression(log[10]~"growth (fg C cell"^{-1}~"h"^{-1}*")"),
       y = expression(log[10]~"respiration (fg C cell"^{-1}~"h"^{-1}*")")) +
  theme_classic(base_size = 12) + theme(legend.position = "none")

# ---- Panel B ----------------------------------------------------------------
bt <- bytax |> dplyr::arrange(slope) |> dplyr::mutate(Taxon = factor(Taxon, levels = Taxon))
pB <- ggplot(bt, aes(slope, Taxon, colour = Taxon)) +
  geom_vline(xintercept = fe_slp, linetype = 2, colour = "grey55") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .7) +
  geom_point(size = 2.4) +
  labs(title = "B", x = "taxon-specific slope (95% CI)", y = NULL) +
  theme_classic(base_size = 12) + theme(legend.position = "none")

ggsave(file.path(figures_dir, "Fig_6_growth_vs_respiration.tiff"),
       pA + pB + patchwork::plot_layout(widths = c(1.4, 1)),
       width = 12, height = 6, dpi = 300, compression = "lzw")
message("Figure 6 written: Fig_6_growth_vs_respiration.tiff (overall slope ", round(fe_slp, 3), ")")
