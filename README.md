# Oxygen-based growth & respiration – analysis scripts

All scripts assume this folder structure:

- `data/`   – input CSVs (raw / trimmed O₂, OD/FC growth, inoculation, temperature data)
- `Tables/` – model outputs, summaries, comparison tables
- `plots/`  – figures for main text and Supplementary

---

## 1. OxygenModel.R

**Purpose:** Main analysis for the multi-taxon experiment.

**What it does:**

- Trims raw SensorDish® O₂ time series to the exponential-decline window (peak → second inflection).
- Normalises oxygen and fits the mechanistic model  
  `O2_norm(t) = O2_0 + (resp_tot / r) * (1 - exp(r * t))`  
  to each Taxon × Replicate (with N₀ fixed to 1 in the model).
- Reconstructs N₀ at O₂ start from inoculation densities and r; converts `resp_tot` to per-cell R (mg O₂ cell⁻¹ min⁻¹).
- Compares O₂-based r with OD₆₀₀ and flow-cytometry (mixed models, boxplots, regressions, Bland–Altman plots).

**Key outputs:**

- `Tables/Oxygen_Data_Trimmed.csv`, `Tables/Oxygen_Data_Filtered.csv`, `Tables/Skipped_Series_Log.csv`
- `Tables/oxygen_model_results.csv`, `Tables/growth_rates_combined.csv`
- `plots/oxygen_trimming_diagnostics.pdf`
- `plots/Fig_2_oxygen_dynamics_facet4_no_overflow.pdf`
- `plots/Fig_3_growth_rate_comparison.pdf`
- `plots/Fig_4_growth_rate_regression_normalized_combined.pdf`
- `plots/Fig_5_BlandAltman_AllReplicates_LME.pdf`

---

## 2. MC.R

**Purpose:** Monte Carlo propagation of N₀ uncertainty into r and R.

**What it does (conceptually):**

- Uses replicate N₀ measurements (e.g. flow cytometry) to estimate mean ± SD per taxon.
- For each oxygen time series, repeatedly samples N₀ from that distribution and re-fits the oxygen model with N₀ fixed.
- Quantifies how much realistic N₀ uncertainty changes r and R (relative SD, bias ranges).
- Summarises the impact on growth and respiration precision.

**Typical outputs:**

- `Tables/N0_MonteCarlo_raw.csv`
- `Tables/N0_MonteCarlo_summary.csv`
- Optional diagnostics in `plots/` (distributions of r and R across MC runs).

---

## 3. Simulation.R

**Purpose:** Synthetic-data validation of the oxygen growth–respiration model.

**What it does:**

- Simulates noisy O₂ trajectories from the same normalised model for known (r, R) across a grid of:
  - growth rates, respiration rates, sampling intervals, durations, noise levels.
- Mimics trimming (stops when O₂ reaches ~40% of its starting value).
- Re-fits the model with `nlsLM()` and quantifies bias, RMSE and relative error in r and R.
- Produces a main figure showing true vs estimated r and R under different noise levels.

**Outputs:**

- `Tables/synthetic_parameter_recovery_results.csv`
- `Tables/synthetic_parameter_recovery_summary.csv`
- `plots/Fig_S2_synthetic_param_recovery_combined.pdf`
- `plots/Fig_S2_synthetic_param_recovery_combined.png`

---

## 4. 05CUTOFF.R

**Purpose:** Reviewer sensitivity check – refitting the model using only O₂ ≥ 0.5 (normalised).

**What it does:**

- Repeats trimming and normalisation as in OxygenModel.R.
- Fits the model twice for each Taxon × Replicate:
  - (i) full trimmed series;  
  - (ii) restricted window where `Oxygen_norm ≥ 0.5`.
- Applies the same QC filters to both fits.
- Computes paired comparisons of r and `resp_tot` (ratios, differences, correlations) and summarises them globally.

**Outputs:**

- `Tables/oxygen_model_results_main.csv`
- `Tables/oxygen_model_results_O2_ge_0.5.csv`
- `Tables/Table_S3_oxygen_model_comparison_O2_ge_0.5.csv`
- `Tables/oxygen_model_comparison_O2_ge_0.5_summary.csv`
- `plots/oxygen_trimming_diagnostics.pdf`
- `plots/oxygen_dynamics_all_models.pdf`
- `plots/oxygen_dynamics_fullsize_per_page.pdf`

---

## 5. CUE.R

**Purpose:** Temperature-gradient analysis for Pseudomonas – growth, respiration and CUE vs temperature.

**What it does:**

- Fits the same normalised oxygen model (N₀ fixed to 1) for each Temperature × Replicate.
- Reconstructs N₀ at O₂ start from a fixed inoculation density and r.
- Converts growth and respiration to carbon units (fg C mL⁻¹ h⁻¹ / min⁻¹) using rod-shaped cell volume and assumed C density.
- Calculates respiration:growth ratios and carbon use efficiency (CUE).
- Fits temperature response curves (two-sided Sharpe–Schoolfield, one-sided SS, Arrhenius), selects by AIC, and extracts Topt, E and Eh with ~95% CIs.
- Produces the 3-panel temperature figure (growth, respiration and CUE vs temperature).

**Outputs:**

- `Tables/oxygen_model_results_good_only_NEWformula.csv`
- `Tables/SharpeSchoolfield_Temperature_Params_NEWformula.csv`
- `plots/oxygen_dynamics_all_models_NEWformula.pdf`
- `plots/oxygen_dynamics_fullsize_per_page_NEWformula.pdf`
- `plots/SharpeSchoolfield_Temperature_Fits_NEWformula.pdf`
