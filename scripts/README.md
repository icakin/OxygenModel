# Refined Oxygen model — bacterial O₂ respiration pipeline (figures only)

Analyses dissolved-oxygen time series for several **bacterial taxa** (e.g.
*Bacillus*, *Burkholderia*, *Arthrobacter*, *Yersinia*). Each series is a
PreSens-style O₂ trajectory; the pipeline trims it to the exponential-decline
window, fits a normalised respiration model, and produces the manuscript
**figures** (plus the intermediate tables the figures are built from).

Organised like the companion Candida pipeline: a single auto-detecting
`config.R`, numbered scripts that each `source()` it, two interactive Shiny apps
(trim-selector + experiment-inputs), a `run_all.R`, and a `data/` + `results/{tables,
figures,rds}` layout.

> Scope: figures-only. The cross-method growth-rate comparison (OD₆₀₀ /
> flow-cytometry vs O₂ — old Fig 3/4/5 and the method-effect stats) is **not**
> ported. The N₀ / per-cell-respiration reconstruction is kept (Fig 6 + the N₀
> Monte-Carlo depend on it).

## The model

Every fit uses the same normalised model (N₀ fixed to 1 inside the fit):

```
O2_norm(t) = O2_0 + (K / r) * (1 - exp(r * t))
```

Per-cell respiration is reconstructed afterwards as `R = K · O2_ref / N0`, with
`N0 = N_inoc · exp(r · delta)` (delta = the inoculation→O₂-start delay that
`02_trimming` records per curve; `N_inoc` and cell carbon come from `config.R`
and the apps).

## Script map

| Script | What it does | Key inputs | Key outputs |
|---|---|---|---|
| `config.R` | Shared paths, model, fit/QC constants, **N_inoc + cell-size inputs**, manual-window/exclusion loader, taxon palette, theme, save helpers | — | sourced by all |
| `01_longdata.R` | Load & tidy the raw long table | `data/Oxygen_Data_Long.csv` | `Oxygen_All_Long.csv` |
| `02_trimming.R` | Hybrid spline start/end trimming (+ per-curve manual overrides); records the N₀ delay | `Oxygen_All_Long.csv` | `Oxygen_Data_Filtered.csv`, `Oxygen_Trimmed_Series_Metadata.csv`, `Oxygen_Curve_Code_Key.csv`, `oxygen_trimming_diagnostics.pdf` |
| `03_trim_selector.R` **(App)** | Click each curve's fit start/end, discard curves; auto-saves the two override CSVs | `Oxygen_All_Long.csv`, `...Series_Metadata.csv` | `manual_fit_windows.csv`, `plot_exclude_points.csv` |
| `04_experiment_inputs.R` **(App)** | Enter per-taxon **Ninoc** (inoculation) + cell size | `Oxygen_All_Long.csv`, `data/Cell_Counts.csv` | `data/Ninoc.csv`, `taxon_cell_sizes.csv` |
| `05_oxygen_fits.R` | nlsLM fits (honours 03 windows/exclusions), N₀/R reconstruction (config N_inoc + per-taxon cell carbon) | `Oxygen_Data_Filtered.csv`, `...Series_Metadata.csv` | `oxygen_fit_curves.csv`, `oxygen_results_with_R.csv`, `oxygen_model_fit_curves.pdf` |
| `06_main_figures.R` | Fig 2, Supp Fig 3, Supp S1, Fig 6 (r vs R RIS) | `oxygen_fit_curves.csv`, `oxygen_results_with_R.csv` | `Fig_2_*.tiff`, `supp_Fig_3_*.tiff`, `Supp_Fig_S1_residuals.tiff`, `Fig_6_*.tiff` |
| `07_cutoff_sensitivity.R` | O₂_norm ≥ 0.5 sensitivity re-fit | `Oxygen_Data_Filtered.csv` | `oxygen_dynamics_all_models.pdf`, `..._fullsize_per_page.pdf` |
| `08_montecarlo_N0.R` | Monte-Carlo of R vs N₀ uncertainty | `oxygen_results_with_R.csv` | `N0_MC_ourmodel_R_rel_sd_*.png/pdf` |
| `09_simulation_recovery.R` | Synthetic parameter recovery | — (synthetic) | `Fig_S2_synthetic_param_recovery_combined.pdf/png` |
| `10_temperature_cue.R` | Temperature-response growth/resp/CUE | `data/Oxygen_Data_Filtered_CUE.csv` | `Fig_7_*_NEWformula.pdf`, `oxygen_dynamics_*_NEWformula.pdf` |

`scripts/original_scripts/` keeps the original flat scripts for reference.

## Interactive trimming / inputs (the "new parts")

- **Trimming** (`02`) uses the hybrid spline detector (oxygen tip → deterministic
  end). Tune the knobs at the top of `02_trimming.R` (`TRIM_SPAR`,
  `START_SEARCH_*`, `END_DROP_FRAC`, `MAX_TIME_MIN`, guards) to your recording
  time-scale.
- **Trim-selector app** (`03`, Run App): click each curve's **Start**/**End**,
  tick *Don't include this sample* to discard, watch the blue model-guide curve.
  It auto-saves `manual_fit_windows.csv` + `plot_exclude_points.csv`; `config.R`
  loads both when `USE_APP_TRIM_FILES <- TRUE`, and `05` fits exactly those
  windows and drops the excluded curves.
- **Inoculation amount (Ninoc)**: set the scalar `N_inoculation_cells_per_L` in
  `config.R`, OR provide a per-curve `data/Ninoc.csv` (columns `Taxon, Replicate,
  N_inoculation_cells_per_L`), OR type it in the **experiment-inputs app**
  (`04`, Run App) — enter each taxon's initial cell count + a conversion factor
  (default 1e6 = counts as cells/µL, the CUE.R convention) and it writes
  `data/Ninoc.csv`. The CSV overrides the scalar per curve.
- **Cell size**: set the global rod dimensions in `config.R`, or use the
  **experiment-inputs app** (`04`, Run App) to give each taxon its own carbon
  per cell (`taxon_cell_sizes.csv`).

## Expected inputs (put these in `data/`)

- `Oxygen_Data_Long.csv` — columns `Taxon, Replicate, Time, Oxygen` (01/02/05/07).
- `Oxygen_Data_Filtered_CUE.csv` — columns `Taxon, Temperature, Replicate, Time,
  Oxygen` (10 only).

Optional inputs:
- `data/Ninoc.csv` — columns `Taxon, Replicate, N_inoculation_cells_per_L` (per-curve
  inoculation density; overrides the config scalar). Can be typed in the `04` app.
- `data/Cell_Counts.csv` — columns `Taxon, Replicate, FC_Initial, ...`; the `04`
  app prefills the Ninoc counts from `FC_Initial`.

(The old `Ninoc_and_deltaTime_to_N0.csv` is no longer required — the N₀ delay now
comes from `02_trimming`, and the density from `Ninoc.csv` / config / the app.)

## How to run

Fully automatic:

```r
source("scripts/run_all.R")               # 01, 02, 05–09
source("scripts/10_temperature_cue.R")    # temperature experiment, separate
```

With hand-tuning:

```r
source("scripts/01_longdata.R")
source("scripts/02_trimming.R")
# RStudio: open 03_trim_selector.R -> Run App   (set windows / exclusions)
# RStudio: open 04_experiment_inputs.R    -> Run App   (optional per-taxon sizes)
source("scripts/05_oxygen_fits.R")
source("scripts/06_main_figures.R")
```

Required packages: `tidyverse`, `zoo`, `minpack.lm`, `patchwork`, `lme4`,
`lmerTest`, `scales`, `stringr`, and (for the apps) `shiny`.
