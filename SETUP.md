# SETUP — getting this repository running on a new machine

Everything below assumes a machine with **R** and **git** and nothing else.
Verified on macOS 26.2 (Tahoe), R 4.5.2, aarch64-apple-darwin20.

---

## 1. Prerequisites

| Tool | Minimum | Why | Install |
|---|---|---|---|
| **R** | 4.5.x | the entire analysis | <https://cloud.r-project.org> |
| **git** | any | to clone, and to record provenance in `env/versions.json` | system package manager |
| **A C++ toolchain** | any | **stage 12 only** — `rstan` compiles a Stan model | macOS: `xcode-select --install`; Linux: `build-essential`; Windows: Rtools |

Stages 01–11 and 13 need no compiler. If you have none, skip stage 12 (see §4).

On Linux you will also need the usual headers that the tidyverse / `lme4` /
`stringi` binaries build against if CRAN has no binary for your distro:

```bash
sudo apt-get install -y build-essential libcurl4-openssl-dev libssl-dev \
                        libxml2-dev libfontconfig1-dev libharfbuzz-dev \
                        libfribidi-dev libfreetype6-dev libpng-dev \
                        libtiff5-dev libjpeg-dev
```

macOS and Windows get CRAN binaries and need none of this.

---

## 2. Clone

```bash
git clone <repo-url> OxygenModel
cd OxygenModel
```

---

## 3. Install — one command

```bash
Rscript scripts/00_install.R
```

This is **idempotent**; re-run it whenever you like. It:

1. installs `renv` if needed, then `renv::restore()`s the project library from
   **`renv.lock`** (150 packages, pinned to exact versions, R 4.5.2);
2. loads every package the pipeline actually uses and stops if any will not load;
3. checks for the C++ toolchain stage 12 needs — a **warning**, not a failure;
4. writes **`env/versions.json`** recording R, platform, `renv` and the git commit.

Expected first-run time: **5–20 min**, almost all of it downloading and (on
Linux) compiling packages. Subsequent runs finish in **~15 s** because
`renv::restore()` finds the library already synchronised.

> **Do not bump `RcppParallel` past 5.1.10.** 6.x ships a TBB release that
> dropped `tbb::task_scheduler_init`, which the StanHeaders code compiled into
> every Stan model still calls; stage 12 then dies at `stan_model()`. The pin is
> deliberate and commented in `scripts/00_install.R`.

---

## 4. Run the pipeline

```bash
bash scripts/run_all.sh
```

Fully unattended — no prompts, no browser windows. Output goes to `results/`
and a timestamped log tree to `logs/<run-id>/`.

To skip the one slow stage (12, the Stan fit) — useful on a machine with no
C++ toolchain:

```bash
OXYMODEL_SKIP_STAN=1 bash scripts/run_all.sh
```

`scripts/run_all.R` is the RStudio-oriented alternative: it `source()`s the
non-app stages into a single session. `run_all.sh` differs in three ways — it
gives each stage its own R process, it covers stages 11–13 (which `run_all.R`
does not list), and it records per-stage wall time and exit status to
`logs/<run-id>/stage_timings.csv`.

---

## 5. Expected runtime and disk

Measured on an Apple M-series laptop, R 4.5.2, warm package cache, with the
Stan model already compiled in the rstan cache.

| Stage | Script | Mode | Wall time |
|---|---|---|---|
| 01 | `01_longdata.R` | run | 2.3 s |
| 02 | `02_trimming.R` | run | 8.0 s |
| 03 | `03_trim_selector.R` | source-only (app **not** launched) | 5.6 s |
| 04 | `04_experiment_inputs.R` | source-only (app **not** launched) | 1.6 s |
| 05 | `05_oxygen_fits.R` | run | 7.1 s |
| 06 | `06_main_figures.R` | run | 10.1 s |
| 07 | `07_cutoff_sensitivity.R` | run | 10.4 s |
| 08 | `08_window_sensitivity.R` | run | 6.1 s |
| 09 | `09_montecarlo_N0.R` | run | 3.4 s |
| 10 | `10_simulation_recovery.R` | run | 16.3 s |
| 11 | `11_temperature_cue.R` | run | 5.4 s |
| 12 | `12_joint_rK_estimator.R` | run (needs rstan) | see note |
| 13 | `13_depletion_frac_sensitivity.R` | run | 2.5 s |
| | **total** | | **≈ 2–4 min** |

Stage 12 dominates. On a **first** run it also has to compile the Stan model,
which adds **40–60 s** on top of sampling; afterwards `rstan_options(auto_write
= TRUE)` caches the compiled model and only the sampling cost remains.
`OXYMODEL_SKIP_STAN=1` removes the stage entirely.

There is no genome-scale model here: 15 taxa × 5 replicates = 75 series.
Expect **a few minutes**, not hours.

Disk:

| Path | Size |
|---|---|
| repo without outputs (`data/` + `scripts/`) | ≈ 1 MB |
| `renv/library` (symlinks into the renv cache) | ≈ 5 MB |
| the renv **cache** itself, on a cold machine | ≈ 1.5 GB (rstan and its Stan headers dominate) |
| one full run (`results/`) | ≈ 30 MB (600-dpi TIFFs dominate) |
| `results/rds/joint_rK_estimator.rds` (stage 12) | ≈ 30 MB |
| `logs/<run-id>/` | ≈ 100 KB |

Budget **≈ 2 GB** on a machine that has neither the renv cache nor a compiler.

---

## 6. Inputs you must have

All live in `data/` and are already in the repository:

| File | Used by |
|---|---|
| `Oxygen_Data_Long.csv` | 01 → 02 → 05 / 07 |
| `Oxygen_Data_Filtered_CUE.csv` | 11 only |
| `Ninoc.csv` | 05, via `config::load_ninoc_table()` |
| `taxon_cell_sizes.csv` | `config::cell_carbon_of()` → 05 |
| `OD_r_FC_r.csv` | 05 (depletion-anchored N0), 06, 13 |
| `Cell_Counts.csv` | the 04 app only — **not** read by any automated stage |

Two further files are **manual curation** and are treated as inputs, never
regenerated by a run:

- `results/tables/manual_fit_windows.csv` — the per-curve fit window set in the
  03 trim-selector app (all 75 curves have one).
- `results/tables/plot_exclude_points.csv` — discarded curves (currently empty).

`scripts/config.R` reads those two from `results/tables/`.

---

## 7. What to check when it finishes

1. `logs/<run-id>/stage_timings.csv` — every stage `ok` (or `skipped` for 12 if
   you set `OXYMODEL_SKIP_STAN=1`).
2. The end of `logs/<run-id>/run_all.log` — the automatic grep for
   `Error` / `Warning` / `skipped` should show only the expected
   `02_trimming: kept 75 curves; skipped 0` line, the `nls.lm`
   `resetting maxiter` warnings, the `boundary (singular) fit` notes from the
   stage 06 mixed models, and the `geom_errorbarh()` deprecation notice.
3. `results/tables/oxygen_fit_results.csv` — 75 rows, `fit_ok` TRUE for all 75.
4. `results/tables/Skipped_Series_Log.csv` — empty (header only).
5. Stage 12 prints `max R-hat` — it should be ≤ 1.01.

---

## 8. Troubleshooting

**A package will not load after `renv::restore()`** — you are missing a system
library; install the headers listed in §1 and re-run step 3.

**Stage 12 fails with `symbol not found in flat namespace
'__ZN3tbb19task_scheduler_init...'`** — an RcppParallel newer than 5.1.10 got
into the library. Put it back:

```bash
R -e 'renv::install("RcppParallel@5.1.10"); renv::restore(prompt = FALSE)'
```

**Stage 12 fails to compile at all** — you have no C++ toolchain. Install one
(§1) or run `OXYMODEL_SKIP_STAN=1 bash scripts/run_all.sh`.

**A run seems to hang** — it should not. The two Shiny apps are guarded and will
not launch in a non-interactive session. If you *want* to open one:

```bash
OXYMODEL_LAUNCH_APPS=1 Rscript scripts/03_trim_selector.R
```
