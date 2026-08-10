# =============================================================================
# D8 PARTS A + B - the scale-free balance quantities, and how far the N0
#                  temperature dependence contaminates them
# =============================================================================
# PART A works from the FITTED parameters r and K only. No conversion to
# per-cell carbon anywhere.
#
# THE ALGEBRA. In 11_temperature_cue.R,
#   growth_C_per_mL  = r * 60 * N0 * q          (q = cell carbon)
#   resp_C_per_mL    = K * O2_ref_mol * 60 * 12e15
#   CUE = growth / (growth + resp) = 1 / (1 + (c/q) * (K * O2_ref) / r)
# so the argmax of CUE is the ARGMIN of (K * O2_ref)/r, and the constant c/q -
# which contains FC_TO_CELLS_PER_L, the carbon quota, RQ and the absolute N0
# scale - drops out entirely.
#
# ONE CORRECTION TO THE BRIEF, and it matters. The brief writes the scale-free
# ratio as K/r. But K is defined on the NORMALISED oxygen trace, so the physical
# consumption rate is K * O2_ref, and O2_ref is oxygen SOLUBILITY, which falls
# with temperature: median 7.21 mg/L at 20 C to 5.31 at 40 C, a factor of 1.36
# across the range. The scale-free quantity is therefore (K * O2_ref)/r. Both are
# computed and reported so the difference is visible.
#
# THE CONDITION, tested not assumed: this holds only if N0 is TEMPERATURE-
# INDEPENDENT. 11_temperature_cue.R sets N0 = 550e3 * exp(r * 45), which is not.
# PART B quantifies exactly what that does.
#
# OUTPUT  runs/D8_analysis/A1_*.csv .. B5_*.csv
# =============================================================================

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d8_common.R"))

message("\n== D8 PARTS A + B ==")

cue <- readr::read_csv(file.path(tables_dir,
        "oxygen_model_results_good_only_NEWformula.csv"), show_col_types = FALSE) |>
  dplyr::filter(fit_ok, is.finite(r_per_minute), r_per_minute > 0,
                is.finite(K), K > 0) |>
  dplyr::mutate(
    r_hr        = r_per_hour,
    K_O2ref     = K * O2_ref,               # physical consumption rate scale
    ratio_K_r   = K / r_per_minute,         # the brief's version
    ratio_phys  = (K * O2_ref) / r_per_minute)  # solubility-corrected
d8_write(cue |> dplyr::select(Taxon, Temperature, Replicate, r_per_minute, r_hr,
                              K, O2_ref, K_O2ref, ratio_K_r, ratio_phys),
         "A1_pseudomonas_scale_free.csv")

# ---- A: activation energies --------------------------------------------------
A2 <- dplyr::bind_rows(
  ba_fit(cue$Temperature, cue$r_hr,    "E_r  (growth rate r)"),
  ba_fit(cue$Temperature, cue$K,       "E_K  (K, normalised)"),
  ba_fit(cue$Temperature, cue$K_O2ref, "E_K  (K * O2_ref, physical)"))
# E_K - E_r with an interval, by bootstrapping curves
set.seed(20260808)
dE <- replicate(2000, {
  i <- sample(nrow(cue), replace = TRUE); d <- cue[i, ]
  if (length(unique(d$Temperature)) < 3) return(NA_real_)
  a <- ba_fit(d$Temperature, d$r_hr, rising_only = TRUE)$E_eV
  b <- ba_fit(d$Temperature, d$K_O2ref, rising_only = TRUE)$E_eV
  b - a
})
dE <- dE[is.finite(dE)]
A2 <- dplyr::bind_rows(A2, tibble::tibble(
  quantity = "E_K - E_r  (physical K)", n = nrow(cue), T_peak_C = NA_real_,
  E_eV = A2$E_eV[3] - A2$E_eV[1],
  E_lo = unname(stats::quantile(dE, .025)), E_hi = unname(stats::quantile(dE, .975)),
  p = NA_real_, r_squared = NA_real_))
d8_write(A2, "A2_activation_energies.csv")

# ---- A: argmin of the balance ratio -----------------------------------------
am_Kr   <- argmin_quad(cue$Temperature, cue$ratio_K_r)
am_phys <- argmin_quad(cue$Temperature, cue$ratio_phys)

tpc <- readr::read_csv(file.path(tables_dir,
        "SharpeSchoolfield_Temperature_Params_NEWformula.csv"), show_col_types = FALSE)
Topt_pipeline <- tpc$Topt_C[tpc$response == "carbon_use_efficiency"]

rng <- cue |> dplyr::group_by(Temperature) |>
  dplyr::summarise(m_Kr = stats::median(ratio_K_r),
                   m_phys = stats::median(ratio_phys), .groups = "drop")

A3 <- tibble::tibble(
  quantity = c("argmin (K/r)  [brief's version]",
               "argmin (K*O2_ref/r)  [solubility-corrected]",
               "pipeline T_opt(CUE) from 11_temperature_cue.R",
               "fold-change in K/r across 20-40 C",
               "fold-change in K*O2_ref/r across 20-40 C",
               "O2_ref at 20 C / O2_ref at 40 C"),
  value = c(am_Kr[["argmin"]], am_phys[["argmin"]], Topt_pipeline,
            max(rng$m_Kr) / min(rng$m_Kr),
            max(rng$m_phys) / min(rng$m_phys),
            stats::median(cue$O2_ref[cue$Temperature == 20]) /
              stats::median(cue$O2_ref[cue$Temperature == 40])),
  lo = c(am_Kr[["lo"]], am_phys[["lo"]], NA, NA, NA, NA),
  hi = c(am_Kr[["hi"]], am_phys[["hi"]], NA, NA, NA, NA))
d8_write(A3, "A3_argmin_and_foldchange.csv")
d8_write(rng, "A4_ratio_by_temperature.csv")

# ---- A: the 15-taxon set at 30 C, ranking only -------------------------------
res <- readr::read_csv(RESULTS_FINAL_CSV, show_col_types = FALSE) |>
  dplyr::filter(fit_ok, is.finite(r_per_minute), r_per_minute > 0, is.finite(K), K > 0)
A5 <- res |>
  dplyr::mutate(ratio_phys = (K * O2_ref) / r_per_minute) |>
  dplyr::group_by(Taxon) |>
  dplyr::summarise(n = dplyr::n(),
                   median_r = stats::median(r_per_minute),
                   median_K_O2ref = stats::median(K * O2_ref),
                   median_ratio_phys = stats::median(ratio_phys),
                   .groups = "drop") |>
  dplyr::arrange(median_ratio_phys) |>
  dplyr::mutate(rank_best_balance = dplyr::row_number())
d8_write(A5, "A5_taxon_ranking_30C.csv")

# =============================================================================
# PART B - the N0 temperature dependence
# =============================================================================
# 11_temperature_cue.R: N0 = 550e3 * exp(r * 45). The 15-taxon pipeline instead
# uses the depletion anchor N0 = FC c exp(-r dt). Both carry r; the question is
# how much that moves the quantities above.
B1 <- cue |>
  dplyr::mutate(
    delta_t_11   = INOC_DELAY_MIN_11,
    expfac_11    = exp(r_per_minute * INOC_DELAY_MIN_11),
    expfac_cfg   = exp(r_per_minute * INOC_DELAY_MIN_CFG),
    expfac_10min = exp(r_per_minute * 10)) |>
  dplyr::select(Temperature, Replicate, r_per_minute, expfac_11, expfac_cfg, expfac_10min)
d8_write(B1, "B1_N0_expfactor_per_curve.csv")

B2 <- B1 |> dplyr::group_by(Temperature) |>
  dplyr::summarise(median_r = stats::median(r_per_minute),
                   median_expfac_45min = stats::median(expfac_11),
                   median_expfac_10min = stats::median(expfac_10min),
                   .groups = "drop")
d8_write(B2, "B2_N0_expfactor_by_temperature.csv")

# CUE under three N0 treatments. Only N0 changes; r, K and O2_ref are the fitted
# values throughout.
mgL_to_mol_per_mL <- function(mg_per_L) (mg_per_L * 1e-3) / 32 / 1000
qC <- cell_carbon_of("Pseudomonas")[1]
cue_variants <- cue |>
  dplyr::mutate(
    N0_45  = INOC_CELLS_PER_uL_11 * 1e3 * exp(r_per_minute * INOC_DELAY_MIN_11),
    N0_0   = INOC_CELLS_PER_uL_11 * 1e3,
    N0_10  = INOC_CELLS_PER_uL_11 * 1e3 * exp(r_per_minute * 10),
    resp_C = K * mgL_to_mol_per_mL(O2_ref) * 60 * 12e15,
    g45 = r_per_minute * 60 * N0_45 * qC,
    g0  = r_per_minute * 60 * N0_0  * qC,
    g10 = r_per_minute * 60 * N0_10 * qC,
    CUE_45 = g45 / (g45 + resp_C),
    CUE_0  = g0  / (g0  + resp_C),
    CUE_10 = g10 / (g10 + resp_C))
d8_write(cue_variants |> dplyr::select(Temperature, Replicate, r_per_minute,
                                       CUE_45, CUE_0, CUE_10),
         "B3_CUE_under_N0_treatments.csv")

argmax_quad <- function(Tc, v) {
  m <- stats::lm(log(v) ~ poly(Tc, 2, raw = TRUE)); b <- unname(stats::coef(m))
  if (b[3] >= 0) return(NA_real_)
  -b[2] / (2 * b[3])
}
B4 <- tibble::tibble(
  treatment = c("N0 = 550e3 * exp(45 r)   [as 11_temperature_cue.R computes it]",
                "N0 = 550e3               [config.R's INOC_DELAY_MIN = 0]",
                "N0 = 550e3 * exp(10 r)   [pre-equilibrated, delta_t = 10 min]"),
  T_opt_CUE = c(argmax_quad(cue_variants$Temperature, cue_variants$CUE_45),
                argmax_quad(cue_variants$Temperature, cue_variants$CUE_0),
                argmax_quad(cue_variants$Temperature, cue_variants$CUE_10)),
  E_K_physical = A2$E_eV[3],
  argmin_ratio_phys = am_phys[["argmin"]],
  E_r = A2$E_eV[1])
d8_write(B4, "B4_Topt_under_N0_treatments.csv")

# The pre-equilibration projection, coldest and warmest.
cold <- min(cue$Temperature); warm <- max(cue$Temperature)
pick <- function(Tv, col) stats::median(B1[[col]][B1$Temperature == Tv])
B5 <- tibble::tibble(
  temperature_C = c(cold, warm, cold, warm),
  delta_t_min = c(45, 45, 10, 10),
  median_expfactor = c(pick(cold, "expfac_11"), pick(warm, "expfac_11"),
                       pick(cold, "expfac_10min"), pick(warm, "expfac_10min"))) |>
  dplyr::mutate(pct_inflation = 100 * (median_expfactor - 1))
# the cold-warm differential is what actually contaminates a thermal comparison
B5 <- dplyr::bind_rows(B5, tibble::tibble(
  temperature_C = NA_real_, delta_t_min = c(45, 10),
  median_expfactor = c(pick(warm, "expfac_11") / pick(cold, "expfac_11"),
                       pick(warm, "expfac_10min") / pick(cold, "expfac_10min")),
  pct_inflation = NA_real_))
d8_write(B5, "B5_preequilibration_projection.csv")

# ---- decomposing the gap between 25.3 C and the published 32.44 C ----------
# Validated first: the CUE reconstructed here with N0 = exp(45r) reproduces the
# pipeline's own carbon_use_efficiency column to 1.3e-7, so the three numbers
# below differ ONLY in the two respects named.
tb_cue <- readr::read_csv(file.path(tables_dir,
        "oxygen_model_results_good_only_NEWformula.csv"), show_col_types = FALSE) |>
  dplyr::filter(fit_ok) |>
  dplyr::select(Temperature, Replicate, cue_pipeline = carbon_use_efficiency) |>
  dplyr::inner_join(dplyr::select(cue_variants, Temperature, Replicate, CUE_45, CUE_0),
                    by = c("Temperature", "Replicate"))

B6 <- tibble::tibble(
  step = c("CUE with N0 temperature-INDEPENDENT, quadratic-in-log fit",
           "CUE as 11_temperature_cue.R computes it (N0 = exp(45 r)), same fit",
           "pipeline's published Sharpe-Schoolfield T_opt(CUE)",
           "  -> shift attributable to the N0 exp(45 r) back-projection",
           "  -> shift attributable to the FUNCTIONAL FORM (SS vs quadratic)",
           "  -> total gap",
           "validation: max |reconstructed CUE_45 - pipeline CUE|"),
  value_C = c(argmax_quad(cue_variants$Temperature, cue_variants$CUE_0),
              argmax_quad(cue_variants$Temperature, cue_variants$CUE_45),
              Topt_pipeline,
              argmax_quad(cue_variants$Temperature, cue_variants$CUE_45) -
                argmax_quad(cue_variants$Temperature, cue_variants$CUE_0),
              Topt_pipeline - argmax_quad(cue_variants$Temperature, cue_variants$CUE_45),
              Topt_pipeline - argmax_quad(cue_variants$Temperature, cue_variants$CUE_0),
              max(abs(tb_cue$CUE_45 - tb_cue$cue_pipeline))))
d8_write(B6, "B6_Topt_gap_decomposition.csv")

message("\n  DECOMPOSING THE T_opt GAP:")
print(as.data.frame(B6), row.names = FALSE, digits = 5)

message("\n  ACTIVATION ENERGIES:")
print(as.data.frame(A2), row.names = FALSE, digits = 4)
message("\n  ARGMIN / FOLD-CHANGE:")
print(as.data.frame(A3), row.names = FALSE, digits = 5)
message("\n  T_opt(CUE) UNDER EACH N0 TREATMENT:")
print(as.data.frame(B4), row.names = FALSE, digits = 5)
message("\n  exp(r*delta_t) BY TEMPERATURE:")
print(as.data.frame(B2), row.names = FALSE, digits = 4)
message("\n  PRE-EQUILIBRATION PROJECTION (last two rows = warm/cold differential):")
print(as.data.frame(B5), row.names = FALSE, digits = 4)
message("\nD8 PARTS A+B done.")
