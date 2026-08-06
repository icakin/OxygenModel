# Oxygen-based bacterial growth & respiration — analysis pipeline

Estimating per-cell growth rate and respiration for 15 bacterial taxa from
single dissolved-oxygen time series, plus a temperature-gradient CUE experiment.

## Folder structure

- `data/`            — input CSVs (raw O₂ series, OD/flow-cytometry counts, cell sizes, temperature data)
- `scripts/`         — the numbered pipeline (below) plus `config.R` and `run_all.R`
- `results/tables/`  — model outputs, summaries, comparison tables
- `results/figures/` — main-text and supplementary figures (TIFF for main figs, PDF for multi-page diagnostics)
- `results/rds/`     — serialised model objects (e.g. the joint estimator fit)

## The model

Each Taxon × Replicate O₂ trace is fit with the normalised model (N₀ fixed to 1 inside the fit):

    O2_norm(t) = O2_0 + (K / r) * (1 - exp(r * t))

- `r`      growth rate (min⁻¹), from the curvature
- `K`      total scaling parameter in the normalised model
- `O2_0`   fitted normalised intercept
- `O2_ref` initial O₂ used for normalisation

## Per-cell respiration and N0

Per-cell respiration is `R = C_tot / (biomass integral)`, where the biomass
integral is `N0 * (exp(r*T_end) - 1) / r` over the fit window.

N0 (cell density at the fit start) is reconstructed from the per-replicate
inoculation densities in `data/Ninoc.csv`:

    N0 = N_inoculation_cells_per_L * exp(r * delta_Ninoc_to_N0_min)

where `delta_Ninoc_to_N0_min` is the inoculation-to-onset delay for that curve.
`r` and the fit windows are independent of N0, so the growth-rate results and the
OD600 / flow-cytometry validation do not depend on it; only per-cell respiration
and CUE do. The absolute respiration/CUE scale rides on the inoculation densities
and the carbon assumptions, not on the slope of the respiration-growth relationship.

## Running the pipeline

From the repo root:

    Rscript scripts/run_all.R

This runs the non-interactive steps in order (01, 02, 05–10). Two optional Shiny
apps (03, 04) let you hand-tune fit windows and per-taxon inputs between 02 and 05.
The temperature analysis (11) and the joint estimator (12) are run separately.

### Scripts

- `01_longdata.R`          raw `Oxygen_Data_Long.csv` → tidy long table (`Oxygen_All_Long.csv`)
- `02_trimming.R`          onset/end trimming, per-curve metadata, diagnostics PDF
- `03_trim_selector.R`     (Shiny app, optional) review/override each curve's fit window + exclusions
- `04_experiment_inputs.R` (Shiny app, optional) per-taxon inputs
- `05_oxygen_fits.R`       nlsLM fits (honours 03 windows) + N0 from Ninoc.csv and per-cell R
- `06_main_figures.R`      Fig 2 (O₂ dynamics), 3/4/5 (cross-method growth), 6 (growth vs respiration, RIS mixed model), Supp Fig 3, Supp Fig S1
- `07_cutoff_sensitivity.R` re-fit using only O₂_norm ≥ 0.5 (reviewer sensitivity)
- `08_window_sensitivity.R` fit-window robustness of K, R, r
- `09_montecarlo_N0.R`     Monte-Carlo of R vs N0 uncertainty
- `10_simulation_recovery.R` synthetic parameter recovery (Fig S2)
- `11_temperature_cue.R`   temperature-gradient growth / respiration / CUE (Fig 7); uses `data/Oxygen_Data_Filtered_CUE.csv`
- `12_joint_rK_estimator.R` joint hierarchical fit propagating the within-curve r–K covariance (requires a working `rstan` / Stan toolchain)

## Key data inputs (`data/`)

- `Oxygen_Data_Long.csv`         raw O₂ series (`Taxon, Replicate, Time, Oxygen`)
- `Ninoc.csv`                    per-replicate inoculation densities and delays (`N_inoculation_cells_per_L`, `delta_Ninoc_to_N0_min`); anchors N0
- `OD_r_FC_r.csv`                per-replicate OD and flow-cytometry counts (`FC_Initial`, `FC_Final`); used for the OD600/FC growth validation (Figs 3–5)
- `taxon_cell_sizes.csv`         per-taxon cell volume → carbon (100 fg C µm⁻³)
- `Cell_Counts.csv`              per-taxon mean of the at-inoculation counts
- `Oxygen_Data_Filtered_CUE.csv` temperature-gradient O₂ series (Pseudomonas)

## Main outputs

- `results/tables/oxygen_results_with_R.csv`  fits + reconstructed N0 + per-cell R and G (carbon units)
- `results/figures/Fig_2 … Fig_7`, `Fig_6_COMBINED.tiff`, and the supplementary figures
- sensitivity and Monte-Carlo tables/figures from 07–10
