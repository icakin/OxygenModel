# =============================================================================
# 14_heldout_and_oxygen_per_cell.R
#   Two checks of the kinetic model that need no carbon assumptions.
# =============================================================================
# CONTEXT. The matched reconciliation FC_Initial vs N0(fit start) sits at ~1.53x.
# Neither RQ nor carbon-per-cell can move that number: they do not appear in
# N0 = FC_Final * FC_TO_CELLS_PER_L * exp(-r * (t_depletion - fit_start)). A
# constant SYBR overcount also cancels, because both sides of the ratio are SYBR
# counts. What remains is
#
#     FC_Initial / N0_hat = (s_i/s_f) * (1/h) * exp(r*dt - G)
#
#   s_i/s_f  event-to-cell mapping at 45 min vs end of run
#   h        cell-number change between depletion and final FC sampling
#   exp()    growth-trajectory error: the back-projection assumes constant r over
#            fit start -> depletion, but the fit window ends with a median ~44.5%
#            of starting O2 remaining, so that interval is largely EXTRAPOLATED.
#
# ln(1.53) = 0.425, i.e. 0.61 of a doubling over the whole extrapolation, so this
# is ordinary model-biology mismatch rather than a calibration error.
#
# CHECK A (held-out trajectory). The model is fitted only to the early window and
# passes every in-window diagnostic (white residuals; nested saturation and O2
# terms rejected at chance level). The region from the window end to depletion is
# therefore genuinely out-of-sample. If constant-r/constant-R over-projects, the
# OBSERVED O2 should sit systematically ABOVE the prediction there. Uses no flow
# cytometry and no carbon assumptions at all: the strongest check available.
#
# CHECK B (oxygen per cell). Under the model, dO2 = R*integral(N) = (R/r)*dN, so
#     Y_kin = R / r          (mg O2 per cell, from the fit)
#     Y_obs = dO2 / dN       (mg O2 per cell, from endpoints)
# should agree. RQ and carbon-per-cell drop out of BOTH sides, so unlike the
# CUE-vs-CUE comparison this does not share the calibrations it is testing.
# Caveat: within the fit window the identity is exact by construction, so the
# test is only meaningful over the FULL run, where the model is extrapolated.
#
# INPUT  results/tables/oxygen_results_with_R.csv, oxygen_fit_curves.csv,
#        Oxygen_All_Long.csv;  data/OD_r_FC_r.csv
# OUTPUT results/tables/heldout_extrapolation_percurve.csv
#        results/tables/oxygen_per_cell_check.csv
#        results/figures/Fig_heldout_extrapolation.png
# RUN    Rscript scripts/14_heldout_and_oxygen_per_cell.R
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
if (exists("FIG_KEEP")) FIG_KEEP <- unique(c(FIG_KEEP, "Fig_heldout_extrapolation"))
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

res  <- readr::read_csv(RESULTS_FINAL_CSV, show_col_types = FALSE)
long <- readr::read_csv(LONG_CSV,          show_col_types = FALSE) %>%
  mutate(Taxon = as.character(Taxon), Replicate = as.character(Replicate))
fcf  <- load_fc_final()
odfc <- readr::read_csv(OD_FC_CSV, show_col_types = FALSE)
names(odfc) <- trimws(names(odfc))

# =============================================================================
# CHECK A - held-out extrapolation
# =============================================================================
message("[A] held-out extrapolation past the fit window ...")

rowsA <- list()
for (i in seq_len(nrow(res))) {
  tx <- res$Taxon[i]; rp <- res$Replicate[i]
  r  <- res$r_per_minute[i]; K <- res$K[i]; O0 <- res$O2_0[i]
  Oref <- res$O2_ref[i]; fs <- res$fit_start_min[i]
  wend <- fs + res$T_end_min[i]          # end of the FITTED window (absolute min)
  tdep <- res$t_depletion_min[i]
  if (!all(is.finite(c(r, K, O0, Oref, fs, wend, tdep))) || tdep <= wend) next

  g <- long %>% filter(Taxon == tx, Replicate == rp, Time >= wend, Time <= tdep) %>%
    arrange(Time)
  if (nrow(g) < 5) next

  # model prediction on the normalised scale, extended past the window
  t0   <- g$Time - fs
  pred <- O0 + (K / r) * (1 - exp(r * t0))
  obs  <- g$Oxygen / Oref
  dev  <- obs - pred                      # >0 means observed ABOVE prediction

  # oxygen remaining at window end, as a fraction of start
  frac_end <- (res$O2_ref[i] * (O0 + (K/r)*(1 - exp(r*res$T_end_min[i])))) / res$O2_start[i]

  rowsA[[length(rowsA)+1]] <- data.frame(
    Taxon = tx, Replicate = rp,
    r_per_minute = r,
    extrap_min       = tdep - wend,
    frac_O2_at_wend  = frac_end,
    mean_dev_norm    = mean(dev, na.rm = TRUE),
    max_dev_norm     = max(dev,  na.rm = TRUE),
    frac_obs_above   = mean(dev > 0, na.rm = TRUE),
    depleted         = isTRUE(res$depleted[i]))
}
A <- bind_rows(rowsA)
readr::write_csv(A, tbl("heldout_extrapolation_percurve.csv"))

cat("\n== CHECK A: held-out region (fit-window end -> depletion) ==\n")
cat(sprintf("  curves tested                : %d\n", nrow(A)))
cat(sprintf("  median extrapolated interval : %.1f min\n", median(A$extrap_min)))
cat(sprintf("  median O2 remaining at window end: %.1f%% of start\n",
            100 * median(A$frac_O2_at_wend, na.rm = TRUE)))
cat(sprintf("  median fraction of held-out points with OBSERVED ABOVE PREDICTED: %.2f\n",
            median(A$frac_obs_above)))
cat(sprintf("  curves where observed is above prediction >75%% of the time: %d / %d\n",
            sum(A$frac_obs_above > 0.75), nrow(A)))
cat("  (systematic 'observed above predicted' = constant-r over-projects growth,\n")
cat("   which under-estimates N0 and inflates FC_Initial/N0 in the observed direction.)\n")

# =============================================================================
# CHECK B - oxygen per new cell, carbon-free
# =============================================================================
message("[B] oxygen per new cell (RQ- and carbon-free) ...")

# total O2 drawdown over the whole recording, per curve (mg/L)
draw <- long %>% group_by(Taxon, Replicate) %>%
  summarise(O2_first = mean(head(Oxygen, 3), na.rm = TRUE),
            O2_last  = mean(tail(Oxygen, 3), na.rm = TRUE),
            .groups = "drop") %>%
  mutate(dO2_mg_per_L = O2_first - O2_last)

B <- res %>%
  select(Taxon, Replicate, r_per_minute, R, N0_cells_per_L) %>%
  left_join(draw, by = c("Taxon", "Replicate")) %>%
  left_join(odfc %>% transmute(Taxon = as.character(Taxon),
                               Replicate = as.character(Replicate),
                               FC_Initial = as.numeric(FC_Initial),
                               FC_Final   = as.numeric(FC_Final)),
            by = c("Taxon", "Replicate")) %>%
  mutate(
    dN_cells_per_L = (FC_Final - FC_Initial) * FC_TO_CELLS_PER_L,
    Y_obs = dO2_mg_per_L / dN_cells_per_L,   # mg O2 per new cell, from endpoints
    Y_kin = R / r_per_minute,                # mg O2 per cell, from the fit
    ratio = Y_obs / Y_kin) %>%
  filter(is.finite(ratio), ratio > 0)

readr::write_csv(B, tbl("oxygen_per_cell_check.csv"))

cat("\n== CHECK B: oxygen per new cell (no RQ, no carbon per cell) ==\n")
cat(sprintf("  curves               : %d\n", nrow(B)))
cat(sprintf("  Y_obs / Y_kin        : median %.2fx  (IQR %.2f - %.2f)\n",
            median(B$ratio), quantile(B$ratio, .25), quantile(B$ratio, .75)))
cat(sprintf("  log2 ratio           : median %.2f\n", median(log2(B$ratio))))
cat("  RQ and carbon-per-cell cancel from both sides, so this tests the kinetic\n")
cat("  oxygen-biomass coupling itself rather than the carbon calibration.\n")

# =============================================================================
# figure
# =============================================================================
p <- ggplot(A, aes(extrap_min, mean_dev_norm)) +
  geom_hline(yintercept = 0, colour = "grey55") +
  geom_point(aes(colour = depleted), size = 2, alpha = .85) +
  scale_colour_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "grey45"),
                      name = "reached 90% depletion") +
  labs(x = "extrapolated interval past the fit window (min)",
       y = "mean (observed - predicted) normalised O2",
       title = "Held-out test: does the fitted model over-project past its window?",
       subtitle = "Points above zero = observed O2 higher than predicted = growth over-projected") +
  theme_bw(base_size = 12)
ggsave(fig("Fig_heldout_extrapolation.png"), p, width = 8.5, height = 5, dpi = 200, bg = "white")
cat("\nwrote heldout_extrapolation_percurve.csv, oxygen_per_cell_check.csv,",
    "Fig_heldout_extrapolation.png\n")
