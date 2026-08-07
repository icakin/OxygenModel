# =============================================================================
# 13_depletion_frac_sensitivity.R - Sensitivity of per-cell respiration (and its
#   between-taxon ordering) to the O2-depletion threshold DEPLETION_FRAC.
# =============================================================================
# N0 is anchored at the time each vial reaches DEPLETION_FRAC of its starting O2
# (default 0.10, i.e. 90% depleted). 13 of 75 curves never reach that threshold
# and fall back to the recording end as the growth-stop, so the choice of
# threshold could bias N0 -> per-cell R -> CUE. This re-derives N0 and per-cell
# respiration across a grid of thresholds and reports, for each: how many curves
# fail to deplete, the median per-cell R, its % shift from the 0.10 default, and
# whether the between-taxon ordering is preserved (Spearman rho vs default).
#
# Only t_depletion changes with the threshold; r, K, C_tot, T_end and the fit
# windows are fixed by 05, so this reuses the committed fits rather than refitting.
#
# INPUT  results/tables/oxygen_fit_results.csv  (from 05)
#        results/tables/Oxygen_All_Long.csv     (full trace -> depletion time)
#        data/OD_r_FC_r.csv                      (FC_Final)
# OUTPUT results/tables/depletion_frac_sensitivity.csv
#        results/figures/Fig_depletion_frac_sensitivity.png
# RUN    Rscript scripts/13_depletion_frac_sensitivity.R
# =============================================================================

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
# config may install a figure whitelist (FIG_KEEP); allow this figure through.
if (exists("FIG_KEEP")) FIG_KEEP <- unique(c(FIG_KEEP, "Fig_depletion_frac_sensitivity"))
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

stopifnot(file.exists(RESULTS_FIT_CSV), file.exists(LONG_CSV))
fits <- readr::read_csv(RESULTS_FIT_CSV, show_col_types = FALSE) %>%
  dplyr::filter(is.finite(r_per_minute), r_per_minute > 0,
                is.finite(C_tot_mg_per_L), C_tot_mg_per_L > 0,
                is.finite(T_end_min), T_end_min > 0)
fcf <- load_fc_final()
if (is.null(fcf)) stop("Need FC_Final (OD_r_FC_r.csv) for the depletion N0 route.")
raw <- readr::read_csv(LONG_CSV, show_col_types = FALSE)
raw <- raw[order(raw$Taxon, raw$Replicate, raw$Time), ]

# depletion time per curve at an arbitrary threshold (mirrors config load_depletion_table)
dep_at <- function(frac) {
  out <- by(raw, list(raw$Taxon, raw$Replicate), function(g) {
    if (nrow(g) < 3) return(NULL)
    o0  <- mean(head(g$Oxygen, 3), na.rm = TRUE)
    hit <- which(g$Oxygen <= frac * o0); reached <- length(hit) > 0
    data.frame(Taxon = as.character(g$Taxon[1]), Replicate = as.character(g$Replicate[1]),
               t_depletion_min = if (reached) g$Time[hit[1]] else max(g$Time, na.rm = TRUE),
               depleted = reached, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# per-cell respiration for every curve at a given threshold
R_at <- function(frac) {
  fits %>%
    dplyr::inner_join(dep_at(frac), by = c("Taxon", "Replicate")) %>%
    dplyr::inner_join(fcf,          by = c("Taxon", "Replicate")) %>%
    dplyr::mutate(
      N0  = dplyr::if_else(is.finite(FC_Final) & FC_Final > 0 &
                             (t_depletion_min - fit_start_min) > 0,
                           n0_depletion(FC_Final, r_per_minute, t_depletion_min, fit_start_min), NA_real_),
      bmi = N0 * (exp(r_per_minute * T_end_min) - 1) / r_per_minute,
      R_C_fg_cell_h = dplyr::if_else(is.finite(bmi) & bmi > 0,
                        (C_tot_mg_per_L / bmi) * O2_to_C_mass * MG_TO_FG * MIN_TO_H, NA_real_)) %>%
    dplyr::filter(is.finite(R_C_fg_cell_h), R_C_fg_cell_h > 0)
}

GRID    <- c(0.05, 0.075, 0.10, 0.125, 0.15, 0.20, 0.25)
DEFAULT <- 0.10

# reference (default) taxon-median R, for ordering + % shift
ref <- R_at(DEFAULT) %>% dplyr::group_by(Taxon) %>%
  dplyr::summarise(R_ref = median(R_C_fg_cell_h), .groups = "drop")

rows <- lapply(GRID, function(fr) {
  d  <- dep_at(fr)
  rr <- R_at(fr)
  tx <- rr %>% dplyr::group_by(Taxon) %>%
    dplyr::summarise(R = median(R_C_fg_cell_h), .groups = "drop") %>%
    dplyr::inner_join(ref, by = "Taxon")
  data.frame(
    depletion_frac              = fr,
    pct_O2_depleted             = round(100 * (1 - fr)),
    n_curves_used               = nrow(rr),
    n_not_depleted              = sum(!d$depleted, na.rm = TRUE),
    median_R_fgC_cell_h         = median(rr$R_C_fg_cell_h),
    median_pct_shift_vs_default = 100 * median(tx$R / tx$R_ref - 1),
    ordering_rho_vs_default     = suppressWarnings(cor(tx$R, tx$R_ref, method = "spearman")))
})
tab <- do.call(rbind, rows)
readr::write_csv(tab, tbl("depletion_frac_sensitivity.csv"))
cat("\n== DEPLETION_FRAC sensitivity (per-cell respiration) ==\n")
print(tab, row.names = FALSE, digits = 3)
cat("\nDefault threshold =", DEFAULT, "(", round(100 * (1 - DEFAULT)),
    "% depleted). Ordering rho ~1 means the between-taxon ranking is preserved.\n")

p <- ggplot(tab, aes(pct_O2_depleted, median_pct_shift_vs_default)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_vline(xintercept = 100 * (1 - DEFAULT), linetype = "22", colour = "grey40") +
  geom_line(colour = "firebrick", linewidth = 1) +
  geom_point(colour = "firebrick", size = 2.4) +
  labs(x = "% O2 depleted at the N0 anchor", y = "median % shift in per-cell R vs default (90%)",
       title = "Sensitivity of per-cell respiration to the depletion threshold",
       subtitle = sprintf("Between-taxon ordering preserved across the grid (Spearman rho %.3f-%.3f); dashed line = default",
                          min(tab$ordering_rho_vs_default), max(tab$ordering_rho_vs_default))) +
  theme_bw(base_size = 12)
ggsave(fig("Fig_depletion_frac_sensitivity.png"), p, width = 8, height = 4.6, dpi = 200, bg = "white")
cat("\nwrote depletion_frac_sensitivity.csv and Fig_depletion_frac_sensitivity.png\n")
