# =============================================================================
# D8 PART C - integrated CUE from mass balance, as an independent check
# =============================================================================
# THE DERIVATION, with every assumption stated.
#
#   biomass carbon produced, per litre
#     = (FC_Final - FC_Initial) [events/uL] * FC_TO_CELLS_PER_L [cells L^-1 per
#       event/uL] * cell_carbon_fg [fg C / cell]                    -> fg C / L
#
#   carbon respired, per litre
#     = total O2 drawdown over the WHOLE run [mg O2 / L]
#       * O2_to_C_mass [mg C per mg O2, = (12.011/31.998) * RQ]
#       * 1e12 [fg per mg]                                          -> fg C / L
#
#   integrated CUE = biomass C / (biomass C + respired C)
#
# VIAL VOLUME CANCELS. Both terms are per litre, so no volume is needed - which
# is fortunate, since none is recorded in the repository.
#
# THE TWO CUEs ARE DIFFERENT ESTIMANDS, and the comparison must say so:
#   kinetic CUE     instantaneous, over the FIT WINDOW only, from r and K
#   integrated CUE  over the WHOLE RUN, including the post-window decline where
#                   growth has stopped and maintenance dominates
# The integrated value should therefore be systematically LOWER, and by more
# where a larger share of the oxygen was consumed after the window closed. That
# is a prediction, and it is tested here rather than asserted.
#
# WHICH INPUTS EACH NEEDS: the kinetic CUE needs N0 (hence FC_TO_CELLS_PER_L and
# the depletion anchor) and the quota; the mass-balance CUE needs
# FC_TO_CELLS_PER_L and the quota but NOT N0 or the back-calculation. Their
# RATIO cancels the quota, so the ratio is the scale-free comparison.
#
# OUTPUT  runs/D8_analysis/C1_*.csv .. C4_*.csv
# =============================================================================

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d8_common.R"))

message("\n== D8 PART C: integrated CUE from mass balance ==")

res  <- readr::read_csv(RESULTS_FINAL_CSV, show_col_types = FALSE) |>
  dplyr::filter(fit_ok)
long <- readr::read_csv(LONG_CSV, show_col_types = FALSE)
odfc <- readr::read_csv(OD_FC_CSV, show_col_types = FALSE)
names(odfc) <- trimws(names(odfc))

# Per-curve oxygen accounting over the whole run, and the post-window share.
acct <- long |>
  dplyr::mutate(Taxon = as.character(Taxon), Replicate = as.character(Replicate)) |>
  dplyr::arrange(Taxon, Replicate, Time) |>
  dplyr::group_by(Taxon, Replicate) |>
  dplyr::summarise(
    t_start_run = min(Time), t_end_run = max(Time), n_points = dplyr::n(),
    O2_start_run = mean(utils::head(Oxygen, 3)),
    O2_end_run   = mean(utils::tail(Oxygen, 3)),
    .groups = "drop") |>
  dplyr::mutate(drawdown_total_mgL = O2_start_run - O2_end_run)

# Oxygen consumed AFTER the fit window closes.
win <- res |> dplyr::transmute(Taxon = as.character(Taxon),
                               Replicate = as.character(Replicate),
                               fit_start_min, T_end_min,
                               win_end_abs = fit_start_min + T_end_min)
post <- long |>
  dplyr::mutate(Taxon = as.character(Taxon), Replicate = as.character(Replicate)) |>
  dplyr::inner_join(win, by = c("Taxon", "Replicate")) |>
  dplyr::group_by(Taxon, Replicate) |>
  dplyr::summarise(
    O2_at_window_end = {
      i <- which.min(abs(Time - dplyr::first(win_end_abs))); Oxygen[i]
    },
    O2_end_post = mean(utils::tail(Oxygen[order(Time)], 3)),
    run_end = max(Time), win_end = dplyr::first(win_end_abs),
    .groups = "drop") |>
  dplyr::mutate(drawdown_post_mgL = pmax(O2_at_window_end - O2_end_post, 0),
                minutes_post_window = pmax(run_end - win_end, 0))

mb <- res |>
  dplyr::select(-dplyr::any_of("FC_Final")) |>   # re-joined from OD_r_FC_r below
  dplyr::mutate(Taxon = as.character(Taxon), Replicate = as.character(Replicate)) |>
  dplyr::inner_join(acct, by = c("Taxon", "Replicate")) |>
  dplyr::inner_join(post,  by = c("Taxon", "Replicate")) |>
  dplyr::inner_join(dplyr::transmute(odfc, Taxon = as.character(Taxon),
                                     Replicate = as.character(Replicate),
                                     FC_Initial = as.numeric(FC_Initial),
                                     FC_Final = as.numeric(FC_Final)),
                    by = c("Taxon", "Replicate")) |>
  dplyr::mutate(
    biomass_C_fg_per_L  = (FC_Final - FC_Initial) * FC_TO_CELLS_PER_L * cell_carbon_fg,
    respired_C_fg_per_L = drawdown_total_mgL * O2_to_C_mass * MG_TO_FG,
    CUE_massbalance = biomass_C_fg_per_L / (biomass_C_fg_per_L + respired_C_fg_per_L),
    CUE_kinetic     = G_C_fg_cell_h / (G_C_fg_cell_h + R_C_fg_cell_h),
    ratio_mb_over_kin = CUE_massbalance / CUE_kinetic,
    frac_O2_post_window = drawdown_post_mgL / drawdown_total_mgL,
    frac_time_post_window = minutes_post_window / (t_end_run - t_start_run),
    physically_impossible = !is.finite(CUE_massbalance) |
      CUE_massbalance < 0 | CUE_massbalance > 1 | (FC_Final - FC_Initial) <= 0)
d8_write(mb |> dplyr::select(Taxon, Replicate, FC_Initial, FC_Final,
                             drawdown_total_mgL, drawdown_post_mgL,
                             biomass_C_fg_per_L, respired_C_fg_per_L,
                             CUE_massbalance, CUE_kinetic, ratio_mb_over_kin,
                             frac_O2_post_window, frac_time_post_window,
                             physically_impossible),
         "C1_massbalance_per_replicate.csv")

C2 <- mb |> dplyr::group_by(Taxon) |>
  dplyr::summarise(n = dplyr::n(),
                   CUE_massbalance = stats::median(CUE_massbalance),
                   CUE_kinetic = stats::median(CUE_kinetic),
                   ratio = stats::median(ratio_mb_over_kin),
                   frac_O2_post_window = stats::median(frac_O2_post_window),
                   n_impossible = sum(physically_impossible),
                   .groups = "drop") |>
  dplyr::arrange(ratio)
d8_write(C2, "C2_massbalance_by_taxon.csv")

ct <- function(x, y) {
  o <- suppressWarnings(stats::cor.test(x, y))
  c(r = unname(o$estimate), p = unname(o$p.value))
}
g1 <- ct(mb$frac_O2_post_window, mb$ratio_mb_over_kin)
g2 <- ct(mb$frac_time_post_window, mb$ratio_mb_over_kin)
C3 <- tibble::tibble(
  quantity = c("median CUE, mass balance (whole run)",
               "median CUE, kinetic (fit window)",
               "median ratio massbalance / kinetic",
               "IQR of that ratio",
               "replicates with CUE outside [0,1] or FC_Final <= FC_Initial",
               "median fraction of O2 consumed AFTER the window closed",
               "corr(ratio, fraction of O2 post-window)",
               "  its p-value",
               "corr(ratio, fraction of TIME post-window)",
               "  its p-value"),
  value = c(stats::median(mb$CUE_massbalance), stats::median(mb$CUE_kinetic),
            stats::median(mb$ratio_mb_over_kin), stats::IQR(mb$ratio_mb_over_kin),
            sum(mb$physically_impossible),
            stats::median(mb$frac_O2_post_window),
            g1[["r"]], g1[["p"]], g2[["r"]], g2[["p"]]))
d8_write(C3, "C3_massbalance_headline.csv")

# ---- CORRECTING THE BRIEF: the ratio does NOT cancel the carbon quota --------
# The brief states that the ratio of the two CUEs cancels the quota. It does not.
# Both are of the form 1/(1 + A/q), so
#   CUE_mb / CUE_kin = (1 + B/q) / (1 + A/q)
# which depends on q. The dependence is tested over a plausible range rather than
# assumed away: the quota is rescaled and the ratio recomputed.
q_scale <- c(0.25, 0.5, 1, 2, 4)
C5 <- dplyr::bind_rows(lapply(q_scale, function(f) {
  B  <- mb$biomass_C_fg_per_L * f
  Cm <- B / (B + mb$respired_C_fg_per_L)
  Gk <- mb$G_C_fg_cell_h * f
  Ck <- Gk / (Gk + mb$R_C_fg_cell_h)
  tibble::tibble(quota_scale = f,
                 median_CUE_massbalance = stats::median(Cm),
                 median_CUE_kinetic = stats::median(Ck),
                 median_ratio = stats::median(Cm / Ck),
                 IQR_ratio = stats::IQR(Cm / Ck))
}))
d8_write(C5, "C5_ratio_sensitivity_to_quota.csv")

C4 <- mb |> dplyr::filter(physically_impossible) |>
  dplyr::select(Taxon, Replicate, FC_Initial, FC_Final,
                biomass_C_fg_per_L, respired_C_fg_per_L, CUE_massbalance)
d8_write(C4, "C4_physically_impossible.csv")

message("\n  HEADLINE:")
print(as.data.frame(C3), row.names = FALSE, digits = 4)
message("\n  BY TAXON:")
print(as.data.frame(C2), row.names = FALSE, digits = 4)
message("\n  SENSITIVITY OF THE RATIO TO THE CARBON QUOTA:")
print(as.data.frame(C5), row.names = FALSE, digits = 4)
if (nrow(C4)) { message("\n  PHYSICALLY IMPOSSIBLE REPLICATES:")
  print(as.data.frame(C4), row.names = FALSE, digits = 4) }
message("\nD8 PART C done.")
