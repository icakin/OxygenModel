# =============================================================================
# 07_supp_growth_resp.R - reproducible SUPPLEMENTARY figure: per-cell growth
#                          vs per-cell respiration (the r-vs-R relationship).
# =============================================================================
# Rebuilds the growth-vs-respiration figure and its data directly from the
# committed pipeline output oxygen_results_with_R.csv, so the figure has a
# generating script and a traceable source. Source config.R first for paths.
# (This is the figure the review flagged as having no producer.)
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

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(lme4) })

res <- readr::read_csv(file.path(tables_dir, "oxygen_results_with_R.csv"),
                       show_col_types = FALSE) |>
  dplyr::filter(is.finite(R_C_fg_cell_h), R_C_fg_cell_h > 0,
                is.finite(G_C_fg_cell_h), G_C_fg_cell_h > 0) |>
  dplyr::mutate(log_R = log(R_C_fg_cell_h), log_G = log(G_C_fg_cell_h))

readr::write_csv(res |> dplyr::select(Taxon, Replicate, log_R, log_G),
                 file.path(tables_dir, "data_growth_vs_resp_pooled.csv"))

pooled <- lm(log_R ~ log_G, data = res)                        # pooled slope
ris    <- lmer(log_R ~ log_G + (1 + log_G | Taxon), data = res) # per-taxon RIS

p <- ggplot(res, aes(log_G, log_R, colour = Taxon)) +
  geom_point(size = 1.6, alpha = .8) +
  geom_smooth(aes(group = Taxon), method = "lm", se = FALSE, linewidth = .4) +
  geom_abline(intercept = coef(pooled)[1], slope = coef(pooled)[2], linewidth = 1) +
  labs(x = "log growth (fg C cell⁻¹ h⁻¹)",
       y = "log respiration (fg C cell⁻¹ h⁻¹)",
       title = "Supplementary figure - per-cell respiration vs growth") +
  theme_bw()

ggsave(file.path(figures_dir, "Supp_Fig_growth_vs_respiration.tiff"),
       p, width = 7, height = 5, dpi = 300)
