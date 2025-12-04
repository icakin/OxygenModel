# Oxygen-based growth & respiration – analysis code

This repository contains the R scripts used to analyse dissolved-oxygen time series with a unified growth–respiration model, compare O₂-based growth rates to OD₆₀₀ and flow cytometry, run synthetic parameter-recovery tests, and quantify temperature effects on growth, respiration, and carbon use efficiency (CUE).

All scripts follow the same basic layout and write their outputs to:

- `data/`   – input CSV files (raw / trimmed oxygen, OD/FC growth rates, inoculation densities, etc.)
- `Tables/` – processed data tables, model coefficients, summaries
- `plots/`  – figures used in the main text and Supplementary Information

You can run the scripts independently, but they are loosely ordered from “core pipeline” to “validation and sensitivity analyses”.

---

## 1. Main O₂ pipeline: trimming, growth & respiration, N₀ from inoculation, cross-method comparison

**Script (example name):** `01_O2_trimming_growth_resp_inocN0.R`  
**Purpose:** Core pipeline for the main manuscript.

**What it does:**

- Automatically trims raw SensorDish® oxygen time series to the steep decline phase (peak to second inflection), dropping noisy stabilisation and late plateau regions.
- Normalises oxygen, fits the mechanistic growth–respiration model  
  `O2_norm(t) = O2_0 + (resp_tot / r) * (1 - exp(r * t))`  
  to each Taxon × Replicate series using `nlsLM()`, with N₀ fixed to 1 in the model.
- Reconstructs N₀ at the start of the O₂ trajectory from inoculation densities and the fitted r, then converts the total scaling parameter (`resp_tot`) into per-cell respiration (mg O₂ cell⁻¹ min⁻¹).
- Produces residual diagnostics (Supp. Fig. S1), 4-panel example fits (Fig. 2), supplementary fit panels, and the cross-method growth comparison (O₂ vs OD₆₀₀ vs FC) including mixed-effects models and Bland–Altman plots.

**Key inputs (in `data/`):**

- `Oxygen_Data_Long.csv` – raw long-format O₂ time series  
- `Inoculation_Density.csv` – inoculation densities (cells/µL)  
- `OD_r_FC_r.csv` – OD₆₀₀ and flow-cytometry growth-rate data

**Key outputs:**

- `Tables/Oxygen_Data_Trimmed.csv`, `Tables/Oxygen_Data_Filtered.csv`, `Tables/Skipped_Series_Log.csv`
- `Tables/oxygen_model_results.csv` (+ intermediate “raw” and QC tables)
- `Tables/growth_rates_combined.csv`, `Tables/method_effects_estimates.csv`, etc.
- `plots/oxygen_trimming_diagnostics.pdf`
- `plots/oxygen_dynamics_all_models.pdf`, `plots/oxygen_dynamics_fullsize_per_page.pdf`
- `plots/Fig_2_oxygen_dynamics_facet4_no_overflow.pdf`
- `plots/Fig_3_growth_rate_comparison.pdf`
- `plots/Fig_4_growth_rate_regression_normalized_combined.pdf`
- `plots/Fig_5_BlandAltman_AllReplicates_LME.pdf`

---

## 2. Synthetic parameter-recovery test

**Script (example name):** `02_synthetic_parameter_recovery.R`  
**Purpose:** Validates the oxygen model by testing how well r and R can be recovered from noisy synthetic DO time series.

**What it does:**

- Uses the same normalised model  
  `O2_norm(t) = O2_0 + (resp_tot / r) * (1 - exp(r * t))`  
  with fixed “true” values of r, per-cell R, N₀ and O₂₀.
- Simulates DO trajectories across combinations of:
  - growth rates (r), respiration rates (R),
  - sampling intervals (dt),
  - total duration,
  - noise levels (Gaussian noise on O₂_norm).
- Mimics the trimming logic (stopping when O₂ drops to ~40% of its starting value).
- Refits the model with `nlsLM()` and summarises bias, RMSE, and relative error for r and R across scenarios.

**Key inputs:**

- None beyond the script itself (fully synthetic).

**Key outputs:**

- `Tables/synthetic_parameter_recovery_results.csv`  
- `Tables/synthetic_parameter_recovery_summary.csv`  
- `plots/Fig_S2_synthetic_param_recovery_combined.pdf`  
- `plots/Fig_S2_synthetic_param_recovery_combined.png`

---

## 3. Reviewer sensitivity test: refitting using only O₂ ≥ 0.5

**Script (example name):** `03_O2_trimming_newmodel_O2_ge_0.5.R`  
**Purpose:** Sensitivity analysis for the reviewer: do r and R change if we only fit the model to oxygen values ≥ 0.5 (normalised)?

**What it does:**

- Repeats the trimming and normalisation of the DO time series (peak → second inflection).
- Fits the same O₂ model (N₀ fixed to 1) to:
  1. The full trimmed series (main analysis).
  2. A restricted window where `Oxygen_norm ≥ 0.5`.
- Applies identical QC criteria to both fits.
- Compares r and `resp_tot` between the two analyses for each Taxon × Replicate pair (ratios, differences, correlations) and summarises the global effect.

**Key inputs:**

- `data/Oxygen_Data_Long.csv`

**Key outputs:**

- `Tables/Oxygen_Data_Trimmed.csv`, `Tables/Oxygen_Data_Filtered.csv`, `Tables/Skipped_Series_Log.csv`
- `Tables/oxygen_model_results_main.csv` (full trimmed window)  
- `Tables/oxygen_model_results_O2_ge_0.5.csv` (O₂_norm ≥ 0.5 only)  
- `Tables/Table_S3_oxygen_model_comparison_O2_ge_0.5.csv`  
- `Tables/oxygen_model_comparison_O2_ge_0.5_summary.csv`
- `plots/oxygen_trimming_diagnostics.pdf`
- `plots/oxygen_dynamics_all_models.pdf`, `plots/oxygen_dynamics_fullsize_per_page.pdf`

---

## 4. Temperature-gradient analysis: growth, respiration & CUE vs temperature

**Script (example name):** `04_temperature_gradient_CUE_TPC.R`  
**Purpose:** Applies the O₂ model to a temperature gradient experiment to quantify how growth, respiration and CUE change with temperature, and to fit Sharpe–Schoolfield / Arrhenius temperature-response curves.

**What it does:**

- Fits the same normalised O₂ model (N₀ fixed to 1) to each Temperature × Replicate time series.
- Reconstructs N₀ at the start of O₂ measurements from a fixed inoculation density and r (using the inoculation delay).
- Converts:
  - growth to fg C mL⁻¹ h⁻¹ and fg C mL⁻¹ min⁻¹,  
  - respiration to fg C mL⁻¹ h⁻¹ and fg C mL⁻¹ min⁻¹,
  using rod-shaped geometry and an assumed carbon density.
- Computes respiration:growth ratios and CUE.
- Fits temperature-response curves (two-sided Sharpe–Schoolfield, one-sided SS, Arrhenius), picks the best model by AIC, and extracts Topt, E and Eh with approximate 95% CIs.
- Produces a 3-panel temperature figure (growth, respiration, and CUE vs temperature with fitted curves).

**Key inputs:**

- `data/Oxygen_Data_Filtered_CUE.csv`  
  (trimmed O₂ data for the temperature experiment; includes Taxon, Temperature, Replicate, Time, Oxygen)

**Key outputs:**

- `Tables/oxygen_model_results_good_only_NEWformula.csv`
- `Tables/SharpeSchoolfield_Temperature_Params_NEWformula.csv`
- `plots/oxygen_dynamics_all_models_NEWformula.pdf`
- `plots/oxygen_dynamics_fullsize_per_page_NEWformula.pdf`
- `plots/SharpeSchoolfield_Temperature_Fits_NEWformula.pdf`

---

## 5. (Optional) Monte Carlo propagation of N₀ uncertainty

**Script (example name):** `05_N0_uncertainty_MonteCarlo.R`  
**Purpose:** Propagates uncertainty in the initial cell density (N₀) into the fitted r and R estimates.

**Typical behaviour (if used):**

- Uses replicate flow-cytometry measurements of inoculation density to estimate mean and SD of N₀ per taxon.
- For each Taxon × Replicate DO series, repeatedly samples N₀ from the empirical normal distribution and refits the oxygen model with N₀ fixed at each draw.
- Quantifies the relative SD of r and R across Monte Carlo runs and summarises how much N₀ uncertainty contributes to parameter uncertainty.

**Expected outputs (convention):**

- `Tables/N0_MonteCarlo_raw.csv` – per-run estimates of r and R  
- `Tables/N0_MonteCarlo_summary.csv` – relative uncertainty per Taxon × Replicate  
- Optional plots in `plots/` showing distributions of r and R across runs

---

## Typical usage order

1. Run the **main O₂ pipeline** to generate core growth/respiration estimates and cross-method comparisons.  
2. Run the **synthetic validation** script to check parameter recovery (optional, but good for methods/appendix).  
3. Run the **O₂ ≥ 0.5 sensitivity** script to reproduce the reviewer analysis.  
4. Run the **temperature-gradient** script to generate CUE and temperature-response results.  
5. (Optional) Run the **N₀ Monte Carlo** script if you want to explicitly quantify how N₀ uncertainty affects r and R.

All scripts assume the `data/`, `Tables/`, and `plots/` folders exist and will create them if needed.
