# =============================================================================
# d2_partB_routes.R - PART B: forward vs backward N0
# =============================================================================
# Both routes estimate N0, the cell density at the start of the fitted oxygen
# window, which cannot be measured in a sealed vial:
#
#   FORWARD  (as published, config.R:157 N0_BACKPROJECT <- TRUE)
#       N0 = N_inoc * exp(r * delta)
#     N_inoc from data/Ninoc.csv, delta = inoculation -> fit-window start.
#     Assumes exponential growth at the oxygen-derived r through the ~45 min
#     stabilization period, an interval the paper itself says the oxygen signal
#     cannot be trusted in.
#
#   BACKWARD
#       N0 = N_end * exp(-r * (t_end - t_0))
#     N_end from the per-replicate endpoint flow-cytometry count of the SAME
#     SensorDish vial, t_0 = fit-window start, t_end = when the endpoint sample
#     was taken. Assumes exponential growth only over the interval the oxygen
#     model already assumes it over.
#
# Both share one assumption. It is tested FIRST, per replicate, against the
# independent flow-cytometry fold change.
#
# Read-only over data/ and runs/D1_baseline/. Writes to runs/D2_analysis/.
# Run:  Rscript reports/D2_n0_routes/scripts/d2_partB_routes.R
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse) })

BASE <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", "..", ".."),
  mustWork = FALSE)
if (!dir.exists(file.path(BASE, "scripts"))) BASE <- getwd()
OUT <- file.path(BASE, "runs", "D2_analysis")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

rd <- function(p) suppressWarnings(suppressMessages(
  readr::read_csv(p, show_col_types = FALSE, progress = FALSE)))

# Same dilution chain as PART A.
FC_DILUTION <- (500 / 490) * (500 / 50) * (502.5 / 500)   # 10.2551
UL_PER_L    <- 1e6

# Exponential-growth consistency window. A replicate is kept for the backward
# route if the OBSERVED flow-cytometry fold change is within this factor of the
# fold change the oxygen-derived r predicts over the same interval. Stated, not
# hidden; the stricter variant is reported as a sensitivity.
GROWTH_CHECK_FACTOR        <- 2.0     # primary   -> keep ratio in [0.5, 2]
GROWTH_CHECK_FACTOR_STRICT <- 1.5     # sensitivity -> keep ratio in [1/1.5, 1.5]

# Carbon conversion, identical to config.R so R_C is comparable to the pipeline.
O2_to_C_mass <- (12.011 / 31.998) * 1
MG_TO_FG <- 1e12
MIN_TO_H <- 60

# =============================================================================
# Load
# =============================================================================
od  <- rd(file.path(BASE, "data", "OD_r_FC_r.csv")); names(od)[1] <- "Taxon"
fit <- rd(file.path(BASE, "runs", "D1_baseline", "tables", "oxygen_results_with_R.csv"))
lng <- rd(file.path(BASE, "runs", "D1_baseline", "tables", "Oxygen_All_Long.csv"))

# t_end: the deposit gives two mutually inconsistent readings, so BOTH are
# computed and the data are used to decide between them.
#  (A) LAST OXYGEN TIMESTAMP (197-212 min). The paper says the endpoint OD600/FC
#      samples were taken "from the corresponding SensorDish vials at the end of
#      the run". The endpoint sample time itself is not recorded in the deposit.
#  (B) THE RECORDED INTERVAL (101-139 min). data/OD_r_FC_r.csv carries a `Time`
#      column that 06_main_figures.R uses as the interval between the initial and
#      final OD/FC samples, i.e. the endpoint sample was taken `Time` minutes
#      after inoculation.
# These cannot both be right. The tie-break is an INTERNAL CONSISTENCY test that
# needs no endpoint at all: project the MEASURED initial flow-cytometry count
# forward to the fit-window start,
#       N0_from_FCstart = N_start * exp(r * t_0),
# and ask which backward estimate it agrees with. Reading (B) agrees to a median
# of 0.99; reading (A) is low by a factor of ~10. (B) is therefore used as the
# PRIMARY backward estimate and (A) is reported as the sensitivity.
run_end <- lng %>%
  dplyr::group_by(Taxon, Replicate) %>%
  dplyr::summarise(t_last_oxygen_min = max(Time), n_oxygen_points = dplyr::n(),
                   .groups = "drop")

d <- od %>%
  dplyr::transmute(Taxon = as.character(Taxon), Replicate = as.character(Replicate),
                   FC_Initial, FC_Final, Duration_recorded_min = Time,
                   OD_Initial, OD_Final) %>%
  dplyr::left_join(run_end, by = c("Taxon", "Replicate")) %>%
  dplyr::left_join(
    fit %>% dplyr::select(Taxon, Replicate, r_per_minute, K, O2_ref, T_end_min,
                          fit_start_min, C_tot_mg_per_L, cell_carbon_fg,
                          N_inoc_cells_per_L, delta_Ninoc_to_N0_min,
                          N0_cells_per_L, R, R_C_fg_cell_h, G_C_fg_cell_h, fit_ok),
    by = c("Taxon", "Replicate")) %>%
  dplyr::mutate(
    t_0_min     = fit_start_min,
    fit_end_min = fit_start_min + T_end_min,
    N_end_cells_per_L   = FC_Final   * FC_DILUTION * UL_PER_L,
    N_start_cells_per_L = FC_Initial * FC_DILUTION * UL_PER_L
  )

# =============================================================================
# B1. The SHARED assumption: is exponential growth at the oxygen-derived r
#     consistent with the independent flow-cytometry fold change?
# =============================================================================
b1 <- d %>%
  dplyr::mutate(
    fold_observed  = FC_Final / FC_Initial,
    fold_predicted = exp(r_per_minute * Duration_recorded_min),
    growth_check_ratio = fold_observed / fold_predicted,
    r_FC_implied   = log(fold_observed) / Duration_recorded_min,
    r_ratio_FC_over_O2 = r_FC_implied / r_per_minute,
    keep_backward = growth_check_ratio >= 1 / GROWTH_CHECK_FACTOR &
                    growth_check_ratio <= GROWTH_CHECK_FACTOR,
    keep_backward_strict = growth_check_ratio >= 1 / GROWTH_CHECK_FACTOR_STRICT &
                           growth_check_ratio <= GROWTH_CHECK_FACTOR_STRICT,
    exclusion_reason = dplyr::case_when(
      growth_check_ratio < 1 / GROWTH_CHECK_FACTOR ~
        sprintf("observed fold change %.2fx is <1/%.1f of the exponential prediction - growth saturated or cells died",
                growth_check_ratio, GROWTH_CHECK_FACTOR),
      growth_check_ratio > GROWTH_CHECK_FACTOR ~
        sprintf("observed fold change %.2fx exceeds the exponential prediction by >%.1fx - r underestimated or count inflated",
                growth_check_ratio, GROWTH_CHECK_FACTOR),
      TRUE ~ NA_character_)
  )
readr::write_csv(
  b1 %>% dplyr::select(Taxon, Replicate, r_per_minute, Duration_recorded_min,
                       FC_Initial, FC_Final, fold_observed, fold_predicted,
                       growth_check_ratio, r_FC_implied, r_ratio_FC_over_O2,
                       keep_backward, keep_backward_strict, exclusion_reason),
  file.path(OUT, "B1_exponential_growth_check.csv"))

b1_summary <- tibble::tibble(
  statistic = c("n replicates", "ratio min", "ratio 5th pct", "ratio 25th pct",
                "ratio median", "ratio 75th pct", "ratio 95th pct", "ratio max",
                "n within a factor of 2", "n within a factor of 1.5",
                "n excluded (primary, factor 2)", "n excluded (strict, factor 1.5)",
                "median r_FC / r_O2"),
  value = c(nrow(b1), min(b1$growth_check_ratio),
            stats::quantile(b1$growth_check_ratio, .05),
            stats::quantile(b1$growth_check_ratio, .25),
            stats::median(b1$growth_check_ratio),
            stats::quantile(b1$growth_check_ratio, .75),
            stats::quantile(b1$growth_check_ratio, .95),
            max(b1$growth_check_ratio),
            sum(b1$keep_backward), sum(b1$keep_backward_strict),
            sum(!b1$keep_backward), sum(!b1$keep_backward_strict),
            stats::median(b1$r_ratio_FC_over_O2)))
readr::write_csv(b1_summary, file.path(OUT, "B1_growth_check_summary.csv"))

b1_taxon <- b1 %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(n = dplyr::n(),
                   fold_obs_med = stats::median(fold_observed),
                   fold_pred_med = stats::median(fold_predicted),
                   ratio_min = min(growth_check_ratio),
                   ratio_med = stats::median(growth_check_ratio),
                   ratio_max = max(growth_check_ratio),
                   n_excluded = sum(!keep_backward), .groups = "drop")
readr::write_csv(b1_taxon, file.path(OUT, "B1_growth_check_by_taxon.csv"))

# =============================================================================
# B2. Both N0 estimates, and their ratio
# =============================================================================
b2 <- b1 %>%
  dplyr::mutate(
    # FORWARD, exactly as 05_oxygen_fits.R computes it
    N0_forward = N_inoc_cells_per_L * exp(r_per_minute * delta_Ninoc_to_N0_min),
    # BACKWARD (PRIMARY), t_end = the recorded OD/FC sampling interval
    interval_min   = Duration_recorded_min - t_0_min,
    N0_backward    = N_end_cells_per_L * exp(-r_per_minute * interval_min),
    # BACKWARD (SENSITIVITY), t_end = last oxygen timestamp
    interval_alt_min = t_last_oxygen_min - t_0_min,
    N0_backward_alt  = N_end_cells_per_L * exp(-r_per_minute * interval_alt_min),
    # CROSS-CHECK that needs no endpoint: project the measured INITIAL count
    # forward to the fit-window start. Independent of t_end entirely.
    N0_from_FCstart = N_start_cells_per_L * exp(r_per_minute * t_0_min),
    ratio_back_over_forward     = N0_backward / N0_forward,
    ratio_backalt_over_forward  = N0_backward_alt / N0_forward,
    ratio_FCstart_over_forward  = N0_from_FCstart / N0_forward,
    check_back_over_FCstart     = N0_backward / N0_from_FCstart,
    check_backalt_over_FCstart  = N0_backward_alt / N0_from_FCstart,
    # how big is each correction itself?
    forward_amplification = exp(r_per_minute * delta_Ninoc_to_N0_min),
    backward_deflation    = exp(-r_per_minute * interval_min)
  )
readr::write_csv(
  b2 %>% dplyr::select(Taxon, Replicate, r_per_minute, delta_Ninoc_to_N0_min,
                       t_0_min, t_last_oxygen_min, Duration_recorded_min,
                       interval_min, interval_alt_min,
                       N_inoc_cells_per_L, N_start_cells_per_L, N_end_cells_per_L,
                       N0_forward, N0_backward, N0_backward_alt, N0_from_FCstart,
                       ratio_back_over_forward, ratio_backalt_over_forward,
                       ratio_FCstart_over_forward,
                       check_back_over_FCstart, check_backalt_over_FCstart,
                       forward_amplification, backward_deflation, keep_backward),
  file.path(OUT, "B2_N0_both_routes_per_replicate.csv"))

kept <- b2 %>% dplyr::filter(keep_backward)

b2_summary <- tibble::tibble(
  quantity = c("n replicates (all)", "n kept for the backward route",
               "N0_forward median (cells/L)",
               "N0_backward median (cells/L), t_end = recorded interval",
               "N0_backward_alt median (cells/L), t_end = last oxygen reading",
               "N0_from_FCstart median (cells/L), no endpoint used",
               "CHECK backward / from_FCstart (median, should be 1)",
               "CHECK backward_alt / from_FCstart (median, should be 1)",
               "ratio backward/forward min", "ratio backward/forward 25th pct",
               "ratio backward/forward MEDIAN", "ratio backward/forward 75th pct",
               "ratio backward/forward max", "ratio fold spread",
               "n with ratio > 1", "n with ratio > 10",
               "forward amplification exp(r*delta) median",
               "backward deflation exp(-r*(t_end-t_0)) median",
               "ratio with t_end = last oxygen reading (median, SENSITIVITY)"),
  value = c(nrow(b2), nrow(kept),
            stats::median(b2$N0_forward), stats::median(kept$N0_backward),
            stats::median(kept$N0_backward_alt), stats::median(kept$N0_from_FCstart),
            stats::median(kept$check_back_over_FCstart),
            stats::median(kept$check_backalt_over_FCstart),
            min(kept$ratio_back_over_forward),
            stats::quantile(kept$ratio_back_over_forward, .25),
            stats::median(kept$ratio_back_over_forward),
            stats::quantile(kept$ratio_back_over_forward, .75),
            max(kept$ratio_back_over_forward),
            max(kept$ratio_back_over_forward) / min(kept$ratio_back_over_forward),
            sum(kept$ratio_back_over_forward > 1),
            sum(kept$ratio_back_over_forward > 10),
            stats::median(b2$forward_amplification),
            stats::median(kept$backward_deflation),
            stats::median(kept$ratio_backalt_over_forward)))
readr::write_csv(b2_summary, file.path(OUT, "B2_N0_ratio_summary.csv"))

b2_taxon <- kept %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(n = dplyr::n(),
                   N0_forward_med = stats::median(N0_forward),
                   N0_backward_med = stats::median(N0_backward),
                   ratio_min = min(ratio_back_over_forward),
                   ratio_med = stats::median(ratio_back_over_forward),
                   ratio_max = max(ratio_back_over_forward),
                   r_med = stats::median(r_per_minute), .groups = "drop") %>%
  dplyr::arrange(ratio_med)
readr::write_csv(b2_taxon, file.path(OUT, "B2_N0_ratio_by_taxon.csv"))

# Does the ratio depend on r? (The forward route makes R an explicit decreasing
# function of r, so a trend here is the diagnostic that matters.)
rfit  <- stats::lm(log(ratio_back_over_forward) ~ r_per_minute, data = kept)
rfit2 <- stats::lm(log(ratio_back_over_forward) ~ log(r_per_minute), data = kept)
b2_vs_r <- tibble::tibble(
  model = c("log(ratio) ~ r", "log(ratio) ~ log(r)"),
  intercept = c(stats::coef(rfit)[1], stats::coef(rfit2)[1]),
  slope     = c(stats::coef(rfit)[2], stats::coef(rfit2)[2]),
  slope_se  = c(summary(rfit)$coefficients[2, 2], summary(rfit2)$coefficients[2, 2]),
  p_value   = c(summary(rfit)$coefficients[2, 4], summary(rfit2)$coefficients[2, 4]),
  r_squared = c(summary(rfit)$r.squared, summary(rfit2)$r.squared),
  spearman_rho = c(stats::cor(kept$ratio_back_over_forward, kept$r_per_minute,
                              method = "spearman"), NA_real_))
readr::write_csv(b2_vs_r, file.path(OUT, "B2_N0_ratio_vs_r.csv"))

# =============================================================================
# B3. Propagate to per-cell respiration R and to CUE
# =============================================================================
# R = C_tot / [N0 * (exp(r*T)-1)/r]  so R is exactly inversely proportional to N0.
b3 <- kept %>%
  dplyr::mutate(
    biomass_integral_fwd = N0_forward     * (exp(r_per_minute * T_end_min) - 1) / r_per_minute,
    biomass_integral_bwd = N0_backward    * (exp(r_per_minute * T_end_min) - 1) / r_per_minute,
    R_forward  = C_tot_mg_per_L / biomass_integral_fwd,
    R_backward = C_tot_mg_per_L / biomass_integral_bwd,
    R_C_fwd = R_forward  * O2_to_C_mass * MG_TO_FG * MIN_TO_H,
    R_C_bwd = R_backward * O2_to_C_mass * MG_TO_FG * MIN_TO_H,
    G_C     = G_C_fg_cell_h,          # independent of N0
    CUE_forward  = G_C / (G_C + R_C_fwd),
    CUE_backward = G_C / (G_C + R_C_bwd),
    R_ratio   = R_backward / R_forward,
    CUE_delta = CUE_backward - CUE_forward
  )
readr::write_csv(
  b3 %>% dplyr::select(Taxon, Replicate, r_per_minute, N0_forward, N0_backward,
                       R_forward, R_backward, R_ratio, R_C_fwd, R_C_bwd, G_C,
                       CUE_forward, CUE_backward, CUE_delta),
  file.path(OUT, "B3_R_and_CUE_per_replicate.csv"))

b3_taxon <- b3 %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(n = dplyr::n(),
                   R_C_fwd_med = stats::median(R_C_fwd),
                   R_C_bwd_med = stats::median(R_C_bwd),
                   R_ratio_med = stats::median(R_ratio),
                   G_C_med = stats::median(G_C),
                   CUE_fwd_med = stats::median(CUE_forward),
                   CUE_bwd_med = stats::median(CUE_backward),
                   CUE_delta_med = stats::median(CUE_delta), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(R_ratio_med))
readr::write_csv(b3_taxon, file.path(OUT, "B3_R_and_CUE_by_taxon.csv"))

# Set the N0-route effect against the window-choice sensitivity D1 measured.
# D1: choosing the manual window rather than 02's automatic one moved K by 23%,
# R by 23% and C_tot by 69% (runs/D1_baseline/comparisons/D2_per_quantity.csv).
D1_WINDOW_R_EFFECT <- 0.2289896   # max relative difference in R, D1 comparison D2
a2 <- rd(file.path(OUT, "A2_three_way_density_per_replicate.csv"))
b3_scale <- tibble::tibble(
  effect = c("fit-window choice (D1: max relative difference in R)",
             "forward back-projection exp(r*delta) alone (median)",
             "N_inoc vs the measured OD600 count (median, PART A)",
             "N_inoc vs the flow-cytometry count (median, PART A)",
             "N0 route: forward -> backward, effect on R (median)",
             "N0 route: forward -> backward, effect on R (max)"),
  fold_change_in_R = c(1 + D1_WINDOW_R_EFFECT,
                       stats::median(b2$forward_amplification),
                       stats::median(a2$ratio_measuredOD_over_Ninoc),
                       stats::median(a2$ratio_FC_over_Ninoc),
                       1 / stats::median(b3$R_ratio),
                       1 / min(b3$R_ratio)),
  direction = c("either", "lowers R", "lowers R", "lowers R", "lowers R", "lowers R"))
readr::write_csv(b3_scale, file.path(OUT, "B3_effect_sizes_vs_window_choice.csv"))

message("d2_partB: wrote B1-B3 artefacts to ", OUT)
