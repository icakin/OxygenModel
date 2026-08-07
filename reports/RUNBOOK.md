# RUNBOOK — OxygenModel, from a clean machine to a full run

Definitive, from-scratch instructions. Reference material, deliberately plain
markdown (see `reports/README.md` for the rule).

If you only want to *install and run*, `SETUP.md` at the project root is
shorter. This runbook goes further: it tells you what to check at each stage and
how to regenerate the D1/D2 comparison artefacts.

> **Status, 2026-08-07.** §2 and §3 describe `scripts/00_install.R`,
> `scripts/run_all.sh`, `renv.lock` and `SETUP.md`. Those arrive with the
> **packaging PR**, which is separate from the change that added this runbook to
> `main` and may not have merged yet. Check with
> `ls scripts/run_all.sh renv.lock`. Until it lands, install packages yourself
> and drive the pipeline with `Rscript scripts/run_all.R` (§3.1). Everything
> else here is accurate against `main` today.

---

## 0. Prerequisites

R ≥ 4.5 and git. A **C++ toolchain** is needed by stage 12 only
(`xcode-select --install` on macOS, `build-essential` on Linux, Rtools on
Windows); stages 01–11 and 13 do not need one.

Quarto ≥ 1.8 and a LaTeX engine are needed **only** if you intend to re-render a
report — see §7, and note that neither report re-renders from a clean clone.

Budget **~2 GB** of disk on a machine with neither an renv cache nor TeX.

Linux only: install build headers first (see `SETUP.md` §1).

---

## 1. Clone

```bash
git clone <repo-url> OxygenModel
cd OxygenModel
```

Everything below is run from that directory.

**Check:** `ls data/` shows six CSVs and a README; `ls results/tables/` shows 34
files; `ls results/figures/` shows 25.

---

## 2. Install

```bash
Rscript scripts/00_install.R
```

Idempotent. Restores the pinned R library from `renv.lock`, verifies every
package the pipeline loads, checks the C++ toolchain stage 12 needs (a warning,
not a failure), and writes `env/versions.json`.

- First run: **5–20 min** (package downloads / compiles).
- Later runs: **~15 s**.

**Check:** the script ends with `environment is pinned and verified.` and
`env/versions.json` exists. If it stops, it prints what to do — it is designed
to fail loudly rather than let you find out later.

> **Do not bump `RcppParallel` past 5.1.10.** 6.x ships a TBB release that
> dropped `tbb::task_scheduler_init`, which the StanHeaders code compiled into
> every Stan model still calls; stage 12 then dies at `stan_model()` with
> `symbol not found in flat namespace '__ZN3tbb19task_scheduler_init...'`.

---

## 3. Run the pipeline

```bash
bash scripts/run_all.sh
```

Zero interaction. **~93 s** total. Output goes to `results/`; logs to
`logs/<run-id>/`.

Stage order — thirteen stages:

```
01 → 02 → (03, 04 sourced only) → 05 → 06 → 07 → 08
   → 09 → 10 → 11 → 12 → 13
```

| | Stage | Notes |
|---|---|---|
| 01 | `01_longdata.R` | |
| 02 | `02_trimming.R` | |
| 03 | `03_trim_selector.R` | Shiny app — **sourced, never launched** |
| 04 | `04_experiment_inputs.R` | Shiny app — **sourced, never launched** |
| 05 | `05_oxygen_fits.R` | fits + depletion-anchored N₀ |
| 06 | `06_main_figures.R` | |
| 07 | `07_cutoff_sensitivity.R` | |
| 08 | `08_window_sensitivity.R` | |
| 09 | `09_montecarlo_N0.R` | |
| 10 | `10_simulation_recovery.R` | |
| 11 | `11_temperature_cue.R` | |
| 12 | `12_joint_rK_estimator.R` | needs rstan + a compiler; slowest stage |
| 13 | `13_depletion_frac_sensitivity.R` | |

To skip the one slow stage on a machine with no compiler:

```bash
OXYMODEL_SKIP_STAN=1 bash scripts/run_all.sh
```

### 3.1 Without `run_all.sh`

`scripts/run_all.R` is the RStudio-oriented alternative — it `source()`s the
non-app stages into one session. It covers **01, 02, 05–10 only**; run 11, 12
and 13 yourself:

```bash
Rscript scripts/run_all.R
Rscript scripts/11_temperature_cue.R
Rscript scripts/12_joint_rK_estimator.R
Rscript scripts/13_depletion_frac_sensitivity.R
```

### 3.2 There is no output-redirection switch

Earlier drafts of this runbook used `OXYMODEL_RESULTS=<dir>` to send a run
somewhere other than `results/`. **`scripts/config.R` on `main` does not read
that variable** — paths are hard-wired to `results/`. A run therefore always
overwrites `results/`. To compare against the committed tree, copy it first:

```bash
cp -R results /tmp/results_committed     # then run, then diff
git status --short results/              # what the run changed
git checkout -- results/                 # put it back
```

**Check, in order:**

1. `logs/<run-id>/stage_timings.csv` — thirteen rows, every `status` = `ok`
   (or `skipped` for 12 if you set `OXYMODEL_SKIP_STAN=1`).
2. `logs/<run-id>/run_all.log`, last section — the automatic grep should show
   only `02_trimming: kept 75 curves; skipped 0`, the `nls.lm`
   `resetting 'maxiter' to 1024` warnings, the `boundary (singular) fit` notes
   from the stage 06 mixed models, and a `geom_errorbarh()` deprecation notice.
   Anything else is new and worth reading.
3. `results/tables/oxygen_fit_results.csv` — 75 rows, `fit_ok` TRUE on all 75.
4. `results/tables/Skipped_Series_Log.csv` — header only, no rows.
5. `results/tables/SharpeSchoolfield_Temperature_Params_NEWformula.csv` —
   3 rows. Missing means stage 11 did not run.
6. Stage 12 prints `max R-hat`; it should be ≤ 1.01.

### Things that are inputs, not outputs

Do **not** regenerate these; the pipeline reads them:

| File | What it is |
|---|---|
| `results/tables/manual_fit_windows.csv` | hand-set fit window for all 75 curves (the `03` app) |
| `results/tables/plot_exclude_points.csv` | discarded curves — currently empty |
| `data/taxon_cell_sizes.csv` | per-taxon cell carbon (the `04` app) |
| `data/Ninoc.csv` | per-replicate inoculation density and delay |
| `data/OD_r_FC_r.csv` | OD and flow-cytometry counts; drives depletion-anchored N₀ |

To edit the first three, open `03`/`04` in RStudio and use *Run App*, or:

```bash
OXYMODEL_LAUNCH_APPS=1 Rscript scripts/03_trim_selector.R
```

---

## 4. Verify the run touched nothing it should not

```bash
git status --short results/ data/
```

`data/` must be empty. `results/` will list what the run rewrote — that is
expected, since §3.2 explains a run always writes there. `git diff` on any file
shows exactly what moved; `git checkout -- results/` restores the committed
tree.

> A previous version of this runbook used
> `shasum -a 256 -c env/pre_run_checksums.txt`. That file is not tracked on
> `main`; git is the checksum record here.

---

## 5. Run the pre-restructure scripts (only for the D1 comparison)

```bash
bash reports/D1_baseline/scripts/run_originals.sh
```

**~49 s.** Copies `scripts/original_scripts/*.R` into `runs/D1_originals/src/`
and runs them from there with their own `data/`, `Tables/`, `plots/`.
`scripts/original_scripts/` is never modified or executed in place — it is the
provenance record of the published analysis, and upstream has since added a
README documenting exactly that.

The scratch tree supplies `data/Ninoc.csv` **also** under the legacy name
`data/Ninoc_and_deltaTime_to_N0.csv`, which `OxygenModel.R` hard-codes. Same
file, two names; no script is edited.

**Check:** `runs/D1_originals/logs/timings.csv` — five rows, all `ok`.

---

## 6. Regenerate the D1 comparison artefacts

```bash
Rscript reports/D1_baseline/scripts/d1_compare.R    # ~25 s
Rscript reports/D1_baseline/scripts/d1_audit.R      # ~15 s
```

Both write only to `runs/D1_baseline/comparisons/`.

These scripts were written when the pipeline could be pointed at
`runs/D1_baseline` via `OXYMODEL_RESULTS`. Without that switch (§3.2) you must
put a run's output there yourself — `cp -R results runs/D1_baseline` after a
run — before `d1_compare.R` will find a modular tree to compare.

**Check:** `runs/D1_baseline/comparisons/` holds 29 CSVs, including
`D1_modular_vs_committed_tables.csv`, `D2_per_quantity.csv`,
`D3_growth_rate_agreement.csv` and `E1`–`E4`.

---

## 7. Re-rendering a report

```bash
quarto render reports/D1_baseline/D1_baseline.qmd --to pdf   # ~30 s
```

**This does not work from a clean clone**, and that is by design. Each `.qmd`
reads comparison CSVs from `runs/` (`runs/D1_baseline/comparisons/`,
`runs/D2_analysis/`) plus `env/versions.json`. `runs/` is **deliberately not
tracked** — 58 MB of non-canonical output trees — so you must regenerate those
artefacts first (§5–§6) before a render will succeed. The committed PDF is the
deliverable; the `.qmd` is committed so the derivation is inspectable, not
because the PDF is expected to be rebuilt on demand.

The `.qmd` recomputes nothing. It renders CSVs, so the numbers in the PDF and
the numbers on disk cannot drift apart.

**Check:** the PDF exists and

```bash
pdftotext reports/D1_baseline/D1_baseline.pdf - | grep -c '??'
```

returns `0` (no unresolved cross-references).

---

## 8. Full sequence, copy-paste

```bash
git clone <repo-url> OxygenModel && cd OxygenModel
Rscript scripts/00_install.R
bash scripts/run_all.sh
git status --short data/            # must be empty
```

Total: **~2 min** after the first install. Add §5–§7 only if you need the D1
comparison or a re-rendered report.

---

## 9. Known issues you will hit

| Symptom | Cause | Fix |
|---|---|---|
| Stage 12: `symbol not found … __ZN3tbb19task_scheduler_init` | RcppParallel newer than 5.1.10 | `R -e 'renv::install("RcppParallel@5.1.10")'` |
| Stage 12 will not compile at all | no C++ toolchain | install one (§0), or `OXYMODEL_SKIP_STAN=1` |
| A package will not load after restore | missing system headers | install them (`SETUP.md` §1), re-run step 2 |
| Terminal appears to hang | should not happen — the apps are guarded | check you did not set `OXYMODEL_LAUNCH_APPS=1` |
| `quarto: command not found` | Quarto not installed | only needed for §7 |
| Render fails on a missing `.sty` | incomplete TeX | `R -e 'tinytex::tlmgr_install("<pkg>")'` |
| Render fails: cannot find `runs/D1_baseline/comparisons/…` | `runs/` is not tracked | regenerate it (§5–§6) |
| `Could not find Oxygen_All_Long.csv` from `03` | `results/tables/` is empty | run step 3 once, or restore the tree from git |

---

## 10. What this repository cannot do

Documented properly in `reports/D1_baseline/D1_baseline.pdf` §6.2, in short:

- **Figure 6 cannot be regenerated.** Six tables, three text files and five
  TIFFs in `results/` have no producing script anywhere in the repository.
- **`results/tables/per_curve_rK_covariance.csv` is committed but no script
  writes it.** Same class of problem, found while packaging.
- **The committed synthetic-recovery tables cannot be reproduced exactly** —
  they were made with a different R / `minpack.lm` version. Both current
  implementations agree with each other exactly and differ from the committed
  files by roughly 10⁻⁵.
- **`results/tables/method_effects_significance.csv` is not reproducible at
  all**, by design: its p-values come from `multcomp::glht`, which is not
  seeded, and move in the fourth decimal place between consecutive runs of
  identical code.
- **`data/Cell_Counts.csv` cannot be reconciled with `data/Ninoc.csv`** from
  inside the repository. See the D2 report and the "Subsequent developments"
  note in `reports/README.md` — that discrepancy is still live.
