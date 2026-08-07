# RUNBOOK — OxygenModel, from a clean machine to the D1 report

Definitive, from-scratch instructions. Reference material, deliberately plain
markdown (see `reports/README.md` for the rule).

If you only want to *install and run*, `SETUP.md` at the project root is
shorter. This runbook goes further: it also regenerates the D1 comparisons and
the report, and tells you what to check at each stage.

---

## 0. Prerequisites

R ≥ 4.5, Quarto ≥ 1.8, git. Nothing else — the install step handles the R
library and, if needed, the LaTeX engine. Budget **~2 GB** of disk on a machine
with neither an renv cache nor a TeX installation.

Linux only: install build headers first (see `SETUP.md` §1).

---

## 1. Clone

```bash
git clone <repo-url> OxygenModel
cd OxygenModel
```

Everything below is run from that directory.

**Check:** `ls data/` shows six CSVs; `ls results/tables/` shows 35 files.

---

## 2. Install

```bash
Rscript scripts/00_install.R
```

Idempotent. Restores the pinned R library from `renv.lock`, verifies every
package the pipeline loads, checks Quarto and a LaTeX engine (installing TinyTeX
if absent), renders a trivial `.qmd` to PDF as a smoke test, and writes
`env/versions.json`.

- First run: **5–20 min** (package downloads / compiles).
- Later runs: **~15 s**.

**Check:** the script ends with `environment is pinned and verified.` and
`env/versions.json` exists. If it stops, it prints what to do — it is designed
to fail loudly rather than let you find out at report-render time.

> The smoke test deliberately uses the real report preamble
> (`reports/_shared/header.tex`), a code block and a table. A hello-world `.qmd`
> renders on LaTeX installations the reports would still fail on — this was an
> actual failure during D1, so the test was strengthened.

---

## 3. Run the pipeline

Writing to `results/` (the default):

```bash
bash scripts/run_all.sh
```

Writing somewhere else, leaving `results/` untouched — **do this if you intend
to compare against the committed outputs**:

```bash
OXYMODEL_RESULTS=runs/D1_baseline bash scripts/run_all.sh
```

Zero interaction. **~63 s** total. Logs land in `logs/<run-id>/`.

Stage order: `01 → 02 → (03, 04 sourced only) → 05 → 06 → 07 → 08 → 09 → 10`.

**Check, in order:**

1. `logs/<run-id>/stage_timings.csv` — ten rows, every `status` = `ok`.
2. `logs/<run-id>/run_all.log`, last section — the automatic grep should show
   only `02_trimming: kept 75 curves; skipped 0` and two `nls.lm`
   `resetting 'maxiter' to 1024` warnings from stage 10. Anything else is new
   and worth reading.
3. `<output>/tables/oxygen_fit_results.csv` — 75 rows, `fit_ok` TRUE on all 75.
4. `<output>/tables/Skipped_Series_Log.csv` — header only, no rows.
5. `<output>/tables/SharpeSchoolfield_Temperature_Params_NEWformula.csv` —
   3 rows (`growth_r_per_hr`, `resp_C_fg_per_hr`, `carbon_use_efficiency`).
   If this file is missing, stage 10 did not run.
6. `<output>/figures/` — 20 files.

### Things that are inputs, not outputs

Do **not** regenerate these; the pipeline reads them:

| File | What it is |
|---|---|
| `results/tables/manual_fit_windows.csv` | hand-set fit window for all 75 curves (the `03` app) |
| `results/tables/plot_exclude_points.csv` | discarded curves — currently empty |
| `data/taxon_cell_sizes.csv` | per-taxon cell carbon (the `04` app) |
| `data/Ninoc.csv` | per-replicate inoculation density and delay |

`scripts/config.R` reads the first two from `results/tables/` even when
`OXYMODEL_RESULTS` points elsewhere. To edit them, open `03`/`04` in RStudio and
use *Run App*, or:

```bash
OXYMODEL_LAUNCH_APPS=1 Rscript scripts/03_trim_selector.R
```

---

## 4. Verify the run touched nothing it should not

```bash
shasum -a 256 -c env/pre_run_checksums.txt
```

Every line should read `OK`. This is a checksum, not an assertion — it is the
only way to know `results/` is untouched.

To re-baseline the checksums after a deliberate change to `results/`:

```bash
find results data -type f -not -name '.DS_Store' | sort | xargs shasum -a 256 > env/pre_run_checksums.txt
```

---

## 5. Run the pre-restructure scripts (only needed for comparisons)

```bash
bash reports/D1_baseline/scripts/run_originals.sh
```

**~49 s.** Copies `scripts/original_scripts/*.R` into `runs/D1_originals/src/`
and runs them from there with their own `data/`, `Tables/`, `plots/`.
`scripts/original_scripts/` is never modified or executed in place — it is the
provenance record of the published analysis.

The scratch tree supplies `data/Ninoc.csv` **also** under the legacy name
`data/Ninoc_and_deltaTime_to_N0.csv`, which `OxygenModel.R` hard-codes. Same
file, two names; no script is edited.

**Check:** `runs/D1_originals/logs/timings.csv` — five rows, all `ok`.

---

## 6. Regenerate the D1 comparisons

```bash
Rscript reports/D1_baseline/scripts/d1_compare.R                                  # ~25 s
OXYMODEL_RESULTS=runs/D1_baseline Rscript reports/D1_baseline/scripts/d1_audit.R  # ~15 s
```

Both write only to `runs/D1_baseline/comparisons/`. The second needs
`OXYMODEL_RESULTS` set so it reads the same tree the run wrote.

**Check:** `runs/D1_baseline/comparisons/` holds 29 CSVs, including
`D1_modular_vs_committed_tables.csv`, `D2_per_quantity.csv`,
`D3_growth_rate_agreement.csv` and `E1`–`E4`.

If you pinned a specific run for the report, refresh the manifest it reads:

```bash
RUNID=$(ls logs | tail -1)
cp "logs/$RUNID/stage_timings.csv" runs/D1_baseline/comparisons/D0_stage_timings.csv
cp runs/D1_originals/logs/timings.csv runs/D1_baseline/comparisons/D0_original_timings.csv
printf 'key,value\nrun_id,%s\nlog_dir,logs/%s\noutput_root,runs/D1_baseline\noriginals_scratch,runs/D1_originals\n' \
  "$RUNID" "$RUNID" > runs/D1_baseline/comparisons/D0_run_manifest.csv
```

---

## 7. Render the report

```bash
quarto render reports/D1_baseline/D1_baseline.qmd --to pdf   # ~30 s
```

**Check:** `reports/D1_baseline/D1_baseline.pdf` exists, is 18 pages, and

```bash
pdftotext reports/D1_baseline/D1_baseline.pdf - | grep -c '??'
```

returns `0` (no unresolved cross-references).

The `.qmd` recomputes nothing. It reads the CSVs from step 6 and
`env/versions.json`, so the numbers in the PDF and the numbers on disk cannot
drift apart.

---

## 8. Full sequence, copy-paste

```bash
git clone <repo-url> OxygenModel && cd OxygenModel
Rscript scripts/00_install.R
OXYMODEL_RESULTS=runs/D1_baseline bash scripts/run_all.sh
shasum -a 256 -c env/pre_run_checksums.txt
bash reports/D1_baseline/scripts/run_originals.sh
Rscript reports/D1_baseline/scripts/d1_compare.R
OXYMODEL_RESULTS=runs/D1_baseline Rscript reports/D1_baseline/scripts/d1_audit.R
quarto render reports/D1_baseline/D1_baseline.qmd --to pdf
```

Total: **~4 min** after the first install.

---

## 9. Known issues you will hit

| Symptom | Cause | Fix |
|---|---|---|
| `quarto: command not found` | Quarto not installed | `brew install --cask quarto`, then step 2 |
| Render fails, `finding package for fvextra.sty` … `Local TeX Live (2025) is older than remote` | the bundled TinyTeX is a release behind and cannot update across releases | `quarto install tinytex`, then step 2 |
| Render fails on a missing `.sty` | incomplete TeX | `R -e 'tinytex::tlmgr_install("<pkg>")'` |
| A package will not load after restore | missing system headers | install them (`SETUP.md` §1), re-run step 2 |
| Terminal appears to hang | should not happen — the apps are guarded | check you did not set `OXYMODEL_LAUNCH_APPS=1` |
| `Could not find Oxygen_All_Long.csv` from `03` | `results/tables/` is empty | `03` always reads the committed `results/tables/`; run step 3 once without `OXYMODEL_RESULTS`, or restore the tree from git |

## 10. What this repository cannot do

Documented properly in `reports/D1_baseline/D1_baseline.pdf` §6.2, in short:

- **Figure 6 cannot be regenerated.** Six tables, three text files and five
  TIFFs in `results/` have no producing script anywhere in the repository.
- **The committed synthetic-recovery tables cannot be reproduced exactly** —
  they were made with a different R / `minpack.lm` version. Both current
  implementations agree with each other exactly and differ from the committed
  files by roughly 10⁻⁵.
- **`data/Cell_Counts.csv` cannot be reconciled with `data/Ninoc.csv`** from
  inside the repository.
