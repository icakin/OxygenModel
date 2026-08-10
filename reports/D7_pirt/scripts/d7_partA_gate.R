# =============================================================================
# D7 PART A - THE GATE: is a regression of R on G measuring physiology, or the
#             shared dependence of both axes on r?
# =============================================================================
# THE ALGEBRA, from scripts/05_oxygen_fits.R under N0_METHOD = "depletion":
#
#   G_C_fg_cell_h = r * 60 * q          q = cell_carbon_fg, a per-TAXON constant
#     => G is EXACTLY proportional to r. Within a taxon, G carries no other
#        information at all.
#
#   C_tot            = (K/r)(e^{rT}-1) O2_ref
#   biomass_integral = N0 (e^{rT}-1)/r
#   R = C_tot / biomass_integral = K * O2_ref / N0      <- the (e^{rT}-1)/r cancel
#   N0 = FC_Final * FC_TO_CELLS_PER_L * exp(-r (t_dep - fit_start))
#     => R = [K O2_ref / (FC c)] * exp(+r * dt),  dt = t_dep - fit_start
#
# So BOTH axes carry r: G linearly, R exponentially through the depletion anchor.
# A positive R-on-G slope could therefore be physiology (Pirt) or simply that
# exp(r dt) and r both increase with r. D6 found an analogous case where an
# apparent predictor turned out to be r in disguise, so this is tested, not
# assumed.
#
# THE TEST: recompute R with N0 held FIXED per taxon - no r dependence in N0 at
# all - and refit. If the Pirt parameters survive, the relationship is not an
# artefact of shared r. If they collapse, there is no decomposition to be had.
#
# OUTPUT  runs/D7_analysis/A1_*.csv .. A3_*.csv
# =============================================================================

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d7_common.R"))

message("\n== D7 PART A: the shared-r gate ==")

res <- readr::read_csv(RESULTS_FINAL_CSV, show_col_types = FALSE) |>
  dplyr::filter(fit_ok, is.finite(G_C_fg_cell_h), G_C_fg_cell_h > 0,
                is.finite(R_C_fg_cell_h), R_C_fg_cell_h > 0)

# Confirm the algebra numerically rather than asserting it.
chk <- res |>
  dplyr::mutate(
    G_predicted = r_per_minute * MIN_TO_H * cell_carbon_fg,
    R_reconstructed = (K * O2_ref / N0_cells_per_L) * O2_to_C_mass * MG_TO_FG * MIN_TO_H,
    dt_min = t_depletion_min - fit_start_min)
A1 <- tibble::tibble(
  check = c("max |G - r*60*cell_carbon| / G",
            "max |R - K*O2ref/N0 (converted)| / R",
            "corr(G, r) across curves",
            "median dt = t_depletion - fit_start (min)",
            "corr(log R, log r) under the default depletion N0",
            "corr(G, R) as the pipeline reports them"),
  value = c(max(abs(chk$G_predicted - chk$G_C_fg_cell_h) / chk$G_C_fg_cell_h),
            max(abs(chk$R_reconstructed - chk$R_C_fg_cell_h) / chk$R_C_fg_cell_h),
            stats::cor(chk$G_C_fg_cell_h, chk$r_per_minute),
            stats::median(chk$dt_min),
            stats::cor(log(chk$R_C_fg_cell_h), log(chk$r_per_minute)),
            stats::cor(chk$G_C_fg_cell_h, chk$R_C_fg_cell_h)))
d7_write(A1, "A1_algebra_check.csv")

# ---- the gate ----------------------------------------------------------------
# Three N0 treatments, all on the SAME K and r, so only the r-dependence of N0
# changes:
#   default   N0 = FC c exp(-r dt)                       (r-dependent)
#   fixed_tax N0 = taxon-mean of the default N0          (no r dependence)
#   fixed_all N0 = grand mean of the default N0          (no r dependence)
gate <- res |>
  dplyr::group_by(Taxon) |>
  dplyr::mutate(N0_fixed_taxon = mean(N0_cells_per_L)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    N0_fixed_all = mean(N0_cells_per_L),
    R_default   = R_C_fg_cell_h,
    R_fixed_tax = (K * O2_ref / N0_fixed_taxon) * O2_to_C_mass * MG_TO_FG * MIN_TO_H,
    R_fixed_all = (K * O2_ref / N0_fixed_all)   * O2_to_C_mass * MG_TO_FG * MIN_TO_H)
d7_write(gate |> dplyr::select(Taxon, Replicate, r_per_minute, K, G_C_fg_cell_h,
                               N0_cells_per_L, N0_fixed_taxon,
                               R_default, R_fixed_tax, R_fixed_all),
         "A2_gate_per_curve.csv")

A3 <- dplyr::bind_rows(
  pirt_fit(gate$G_C_fg_cell_h, gate$R_default,   "N0 default (r-dependent)"),
  pirt_fit(gate$G_C_fg_cell_h, gate$R_fixed_tax, "N0 fixed per taxon"),
  pirt_fit(gate$G_C_fg_cell_h, gate$R_fixed_all, "N0 fixed globally"))
d7_write(A3, "A3_gate_pirt_fits.csv")

# ---- the same gate for the Pseudomonas temperature series --------------------
# 11_temperature_cue.R has its OWN INOC_DELAY_MIN <- 45 (line 40), overriding
# config.R's 0, so its N0 = 550 * exp(r * 45) carries r as well - the shared-r
# worry applies to B2 too and is tested identically.
#
#   growth_C_fg_per_hr = r * 60 * N0_cells_per_mL * cell_C     (per mL)
#   resp_C_fg_per_hr   = K * O0_mol_per_mL * 60 * 12e15        (per mL, no N0)
#
# Pirt needs SPECIFIC rates, so both are divided by biomass:
#   mu = r * 60                                   (N0 cancels exactly)
#   q  = resp_per_mL / (N0 * cell_C)  ->  proportional to K / exp(45 r)
# Holding N0 fixed removes that exp(45r), exactly as the taxon test above does.
cue <- readr::read_csv(file.path(tables_dir, "oxygen_model_results_good_only_NEWformula.csv"),
                       show_col_types = FALSE) |>
  dplyr::filter(fit_ok, is.finite(r_per_minute), r_per_minute > 0,
                is.finite(resp_C_fg_per_hr), resp_C_fg_per_hr > 0)

INOC_CELLS_PER_uL_11 <- 550   # as set in 11_temperature_cue.R
INOC_DELAY_MIN_11    <- 45    # ditto - NOT config.R's 0
cue <- cue |>
  dplyr::mutate(
    biomass_default = N0_cells_per_mL * cell_carbon_of(Taxon),
    N0_fixed_mL     = INOC_CELLS_PER_uL_11 * 1e3,
    biomass_fixed   = N0_fixed_mL * cell_carbon_of(Taxon),
    mu_per_hr       = r_per_hour,
    q_default       = resp_C_fg_per_hr / biomass_default,
    q_fixed         = resp_C_fg_per_hr / biomass_fixed)
d7_write(cue |> dplyr::select(Taxon, Temperature, Replicate, r_per_hour, K,
                              N0_cells_per_mL, mu_per_hr, q_default, q_fixed),
         "A4_gate_pseudomonas_per_curve.csv")

A5 <- dplyr::bind_rows(
  pirt_fit(cue$mu_per_hr, cue$q_default, "Pseudomonas, N0 default (exp(45r))"),
  pirt_fit(cue$mu_per_hr, cue$q_fixed,   "Pseudomonas, N0 fixed (no r)"))
d7_write(A5, "A5_gate_pseudomonas_fits.csv")

message("\n  ALGEBRA:")
print(as.data.frame(A1), row.names = FALSE, digits = 4)
message("\n  THE GATE - Pirt fit under each N0 treatment:")
print(as.data.frame(A3[, c("dataset", "n", "m_intercept", "slope_inv_Y",
                           "slope_p", "Y", "r_squared", "corr_G_R")]),
      row.names = FALSE, digits = 4)
message("\n  THE GATE - Pseudomonas (within-taxon, across temperature):")
print(as.data.frame(A5[, c("dataset", "n", "m_intercept", "m_lo", "m_hi",
                           "slope_inv_Y", "slope_p", "Y", "r_squared", "corr_G_R")]),
      row.names = FALSE, digits = 4)
message("\nD7 PART A done.")
