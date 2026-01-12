# Oxygen-based growth & respiration — analysis scripts (updated)

All scripts assume this folder structure:

- `data/`   — input CSVs (raw/trimmed O₂, OD/FC growth, inoculation + timing tables, temperature data)
- `Tables/` — model outputs, summaries, comparison tables
- `plots/`  — figures for main text + Supplementary (multi-page diagnostics stay as PDF; main figs now also saved as TIFF)

---

## Model used throughout (normalised)

All fitting scripts use the same normalised “OUR model” (with N₀ fixed to 1 inside the fit):

O2_norm(t) = O2_0 + (K / r) * (1 - exp(r * t))

Where:

- `r` = growth rate (min⁻¹)
- `K` = total scaling parameter in the normalised model
- `O2_0` = fitted intercept in normalised O₂
- `O2_ref` = initial O₂ used for normalisation (typically mean of first few O₂ points after onset trimming)

Per-cell respiration is reconstructed after fitting using:

R = (K * O2_ref) / N0

(`R` in mg O₂ cell⁻¹ min⁻¹; `N0` in cells L⁻¹)

---

## 1) `OxygenModel.R`

Purpose: Full pipeline for the multi-taxon experiment (trimming → fitting → reconstruct N₀ & R → comparisons + figures).

What it does:

- Trims raw SensorDish® O₂ time series to the exponential-decline window (peak → second inflection), logs skipped series, and writes a trimming diagnostics PDF.
- Fits the normalised model per Taxon × Replicate, with QC thresholds (SE, p-values, pseudo-R², residual range, AIC gain vs intercept-only, MAPE, parameter bounds).
- Exports fitted curves (per timepoint) and a multi-page PDF of per-series fit plots.
- Reads `data/Ninoc_and_deltaTime_to_N0.csv` and reconstructs N₀ at O₂ start from inoculation density + delay:
  - N0 = N_inoc * exp(r * delta_t)
  - then computes per-cell respiration R and additional derived quantities.
- Produces cross-method growth-rate comparisons (OD₆₀₀ vs O₂; FC vs O₂), mixed models, regressions, Bland–Altman, and summary figures.

Key inputs (in `data/`):

- `Oxygen_Data_Long.csv` (raw O₂ time series; must include `Taxon, Replicate, Time, Oxygen`)
- `Ninoc_and_deltaTime_to_N0.csv` (inoculation density + delay to O₂ start; expected columns include `Taxon, Replicate, N_inoculation_cells_per_L, delta_Ninoc_to_N0_min`)
- (Plus your OD/FC growth inputs used in the comparison sections)

Key outputs:

Tables:
- `Tables/Oxygen_Data_Trimmed.csv`
- `Tables/Oxygen_Data_Filtered.csv`
- `Tables/Skipped_Series_Log.csv`
- `Tables/oxygen_fit_results.csv` (fit parameters + QC)
- `Tables/oxygen_results_with_R.csv` (adds reconstructed `N0` + per-cell `R`)
- `Tables/oxygen_fit_curves.csv` (time-resolved fitted curves)
- Additional figure-specific tables (e.g., taxon-specific slope tables)

Plots:
- `plots/oxygen_trimming_diagnostics.pdf` (multi-page)
- `plots/oxygen_model_fit_curves.pdf` (multi-page, per-series fits)
- Main figures saved as TIFF, including:
  - `plots/Fig_2_oxygen_dynamics_facet4_no_overflow.tiff`
  - `plots/Fig_3_growth_rate_comparison.tiff`
  - `plots/Fig_4_growth_rate_regression_normalized_combined.tiff`
  - `plots/Fig_5_BlandAltman_AllReplicates_LME.tiff`
  - `plots/Fig_6_r_vs_R_RIS_derand_full_MAIN.tiff`
  - plus supplementary TIFFs (residuals, all-replicate normalised plots, etc.)

Notes on figure order: RIS is now Fig 6 and is placed after Fig 5; TIFF export is used for main/supp figures, while multi-page diagnostics remain PDF.

---

## 2) `MC.R`

Purpose: Monte Carlo sensitivity of reconstructed R to plausible N₀ uncertainty (minimal-output version).

What it does:

- Takes `Tables/oxygen_results_with_R.csv` and, for each Taxon × Replicate:
  - uses that replicate’s fitted N0 as the mean (`N0_mu`)
  - samples N0 from a chosen distribution (default: lognormal) across a grid of CV values
  - propagates uncertainty into `R = (K * O2_ref) / N0` and summarises relative SD / quantiles.
- If `O2_ref` is missing in `oxygen_results_with_R.csv`, it can compute it from `Tables/Oxygen_Data_Filtered.csv`.

Inputs:
- `Tables/oxygen_results_with_R.csv` (must contain `Taxon, Replicate, N0_cells_per_L` and `K`; `O2_ref` optional)
- `Tables/Oxygen_Data_Filtered.csv` (only used if `O2_ref` needs reconstructing)

Default outputs (only 2 CSVs):
- `Tables/N0_MC_ourmodel_taxon_summary_ALL.csv`
- `Tables/N0_MC_ourmodel_overall_summary_ALL.csv`

Optional (off by default via flags in the script):
- combined per-series summaries and a couple of diagnostic plots

---

## 3) `Simulation.R`

Purpose: Synthetic-data validation (parameter recovery) for the same normalised O₂ model.

What it does:

- Simulates noisy O₂ trajectories from `O2_norm(t)` across grids of:
  - true r, true per-cell R, sampling interval, duration, noise SD
- Mimics trimming by stopping when O₂ reaches ~40% of its starting value.
- Re-fits the model via `nlsLM()` and quantifies bias, RMSE and relative error in r and R.

Outputs:
- `Tables/synthetic_parameter_recovery_results.csv`
- `Tables/synthetic_parameter_recovery_summary.csv`
- Parameter-recovery plots in `plots/` (PDF/PNG as defined in the script)

---

## 4) `05CUTOFF.R`

Purpose: Reviewer sensitivity check — re-fit using only O₂_norm ≥ 0.5.

What it does:

1) Trims raw O₂ time series (same general trimming logic as the main pipeline) and writes diagnostics.
2) Fits the model twice per Taxon × Replicate:
   - (i) full trimmed window
   - (ii) restricted window where `Oxygen_norm ≥ 0.5`
3) Compares paired estimates and exports ratios/differences for r and K.

Inputs (in `data/`):
- `Oxygen_Data_Long.csv`

Outputs:

Tables:
- `Tables/Oxygen_Data_Trimmed.csv`
- `Tables/Oxygen_Data_Filtered.csv`
- `Tables/Skipped_Series_Log.csv`
- `Tables/oxygen_model_results_main.csv`
- `Tables/oxygen_model_results_O2_ge_0.5.csv`
- `Tables/Table_S3_oxygen_model_comparison_O2_ge_0.5.csv`
- `Tables/oxygen_model_comparison_O2_ge_0.5_summary.csv`

Plots:
- `plots/oxygen_trimming_diagnostics.pdf`
- `plots/oxygen_dynamics_all_models.pdf`
- `plots/oxygen_dynamics_fullsize_per_page.pdf`

---

## 5) `CUE.R`

Purpose: Temperature-gradient analysis (Pseudomonas): growth, respiration and CUE vs temperature.

What it does:

- Fits the same normalised O₂ model per Temperature × Replicate.
- Reconstructs N0 at O₂ start from a fixed inoculation density + delay and fitted r.
- Converts growth + respiration to carbon units and computes respiration:growth and CUE.
- Fits temperature-response curves (Sharpe–Schoolfield / Arrhenius), selects by AIC, and extracts parameters with ~95% CIs.

Input (in `data/`):
- `Oxygen_Data_Filtered_CUE.csv`

Outputs:

Tables:
- `Tables/oxygen_model_results_good_only_NEWformula.csv`
- `Tables/SharpeSchoolfield_Temperature_Params_NEWformula.csv`

Plots:
- `plots/oxygen_dynamics_all_models_NEWformula.pdf`
- `plots/oxygen_dynamics_fullsize_per_page_NEWformula.pdf`
- `plots/Fig_7_SharpeSchoolfield_Temperature_Fits_NEWformula.pdf`
