# runs/ — manifest

Output trees produced by the D-series prompts. **None of these is canonical.**
`results/` is the committed record of the published analysis and is never
written to by anything here; every tree below is a parallel regeneration used as
evidence in a report.

## Tracking policy

**Exactly the tables the two committed reports read at render time are tracked.
Nothing else is.** That is a narrower rule than "all tables": it is the minimum
set that lets `reports/D1_baseline/D1_baseline.qmd` and
`reports/D2_n0_routes/D2_n0_routes.qmd` rebuild from a clean clone, which is the
only reason any of `runs/` is in the repository at all.

The set was determined by parsing every `rd()` and `file.path(CMP, ...)` call in
both sources — not by taking whole directories. Figures (600-dpi TIFFs,
multi-page diagnostic PDFs) and `.rds` model objects are ~47 MB and are always
excluded; a full run regenerates them. See `.gitignore`, section `runs/`.

| | on disk | tracked |
|---|---|---|
| all of `runs/` | ~47 MB | **133 KB** (48 files) |

### What is tracked, per tree

| Subtree | tracked | excluded | why excluded |
|---|---:|---:|---|
| `D1_baseline/comparisons/` | **24** | 5 | not read by `D1_baseline.qmd` |
| `D1_baseline/tables/` | 0 | 24 | a run's output; the report reads `comparisons/`, not this |
| `D1_baseline/figures/`, `rds/` | 0 | all | bulk |
| `D1_originals/Tables/` | 0 | 25 | consumed by `d1_compare.R`, not by any `.qmd` |
| `D1_originals/plots/`, `data/`, `src/`, `logs/` | 0 | all | bulk, or copies of tracked repo content |
| `D2_analysis/` | **24** | 4 | not read by `D2_n0_routes.qmd` |
| `D2_forward/`, `D2_backward/`, `D2_none/` | 0 | 2 each | their numbers reach the report via `D2_analysis/C1`–`C4` |

The five excluded `D1_baseline/comparisons/` files are
`D1_modular_vs_committed_columns.csv`, `D2_modular_vs_original_columns.csv`,
`D3_optima.csv`, `E3_qc_gate_per_series.csv` and
`E4_backprojection_per_series.csv`. The four excluded `D2_analysis/` files are
`A3_Ninoc_within_taxon_correlations.csv`, `A4_bland_altman_by_taxon.csv`,
`A4_growth_rates_perrep_vs_published.csv` and `B3_R_and_CUE_per_replicate.csv`.

**Consequence to be aware of:** the reports rebuild, but the *comparison
scripts* (`reports/D1_baseline/scripts/d1_compare.R` and the D2 `scripts/`) do
not — they need the excluded run outputs. Re-deriving the tables from scratch
means re-running the pipeline first; `reports/RUNBOOK.md` §5–§6 gives the
sequence. Rendering the reports needs none of that.

Per-stage timing logs are *not* tracked (`logs/` is ignored), but the two the D1
report cites were copied into `runs/D1_baseline/comparisons/` and those are
tracked.

> The `OXYMODEL_RESULTS=...` invocations recorded below are how these trees were
> produced **at the time**, on the D1/D2 work branch. `scripts/config.R` on
> `main` does not read that variable — see `reports/RUNBOOK.md` §3.2.

## Environment

Every tree below was produced under the same pinned environment, recorded in
`env/versions.json`:

R 4.5.2 (2025-10-31) · aarch64-apple-darwin20 · macOS Tahoe 26.2 ·
`renv.lock` with 135 packages · Quarto 1.8.27 · lualatex (TeX Live 2026).

---

## `D1_baseline/` — 24 files tracked (`comparisons/` only)

| | |
|---|---|
| **Produced by** | `OXYMODEL_RESULTS=runs/D1_baseline bash scripts/run_all.sh`, run id `20260726_211928` |
| **Prompt** | `prompts/D1_baseline_reproduction_prompt.md` (PART C) |
| **Commit** | run under `2c4315d`, recorded in `3b01180` |
| **Date** | 2026-07-26 |
| **Tracked** | `comparisons/` — the 24 files `D1_baseline.qmd` reads |
| **Excluded** | `tables/` (24), 5 unread `comparisons/` files, `figures/` (20 files, 21 MB), `rds/` |

**Underpins:** the whole of `reports/D1_baseline/D1_baseline.pdf`. Specifically
that the modular pipeline reproduces 19 of the 22 regenerable committed tables
as IDENTICAL, 10 of them at exactly 0.0 (`comparisons/D1_modular_vs_committed_tables.csv`);
that all 75 series pass all eight QC criteria (`comparisons/E3_qc_gate_by_criterion.csv`);
and the paper's headline numbers — R² 0.944 vs OD600 and 0.968 vs flow
cytometry, growth *T*<sub>opt</sub> 33.98 °C, CUE *T*<sub>opt</sub> 32.44 °C
(`comparisons/D3_*.csv`).

---

## `D1_originals/` — nothing tracked

| | |
|---|---|
| **Produced by** | `bash reports/D1_baseline/scripts/run_originals.sh` |
| **Prompt** | `prompts/D1_baseline_reproduction_prompt.md` (PART D, comparison D2) |
| **Commit** | `9084beb` |
| **Date** | 2026-07-26 |
| **Tracked** | none — no `.qmd` reads this tree directly |
| **Excluded** | `Tables/` (25, consumed by `d1_compare.R`), `plots/` (24 MB), `data/` and `src/` (copies of tracked repo content), `logs/` |

A scratch tree in which **copies** of `scripts/original_scripts/*.R` were run.
The originals themselves were never modified or executed in place — they are the
provenance record of the published analysis.

**Underpins:** that all five pre-restructure scripts still run unmodified; that
the temperature/CUE and synthetic-recovery ports are exact to 0.0 while the
bacterial O₂ fits diverge (*r* by 0.63%, *K* by 23%, *C*<sub>tot</sub> by 69%);
and — via `D1_baseline/comparisons/D2b_original_vs_committed_tables.csv` — that
`results/tables/` was produced by the **modular** pipeline rather than by these
scripts.

---

## `D2_analysis/` — 24 of 28 files tracked

| | |
|---|---|
| **Produced by** | `reports/D2_n0_routes/scripts/{d2_partA_reconcile,d2_partB_routes,d2_partC_temperature}.R` |
| **Prompt** | `prompts/D2_n0_reconciliation_and_dual_route_prompt.md` |
| **Commit** | `990ab09` (A), `3fbee63` (B), `9f2ca10` (C) |
| **Date** | 2026-07-26 |
| **Tracked** | the 24 files `D2_n0_routes.qmd` reads |
| **Excluded** | 4 unread files (listed under Tracking policy) |

**Underpins:** the whole of `reports/D2_n0_routes/D2_n0_routes.pdf`. The
three-way biomass reconciliation (`A2_*`), the finding that
`N_inoculation_cells_per_L`'s provenance cannot be recovered (`A3_*`), the
exponential-growth check (`B1_*`), the forward-versus-backward *N*₀ ratio of 57
(`B2_*`), the effect on *R* and CUE (`B3_*`), and the three-route temperature
comparison (`C1`–`C4`).

---

## `D2_forward/`, `D2_backward/`, `D2_none/` — nothing tracked

| | |
|---|---|
| **Produced by** | `OXYMODEL_RESULTS=runs/D2_<route> OXYMODEL_N0_ROUTE=<route> Rscript scripts/10_temperature_cue.R` |
| **Prompt** | `prompts/D2_n0_reconciliation_and_dual_route_prompt.md` (PART C) |
| **Commit** | `9f2ca10` |
| **Date** | 2026-07-26 |
| **Tracked** | none — these routes reach the report through `D2_analysis/C1`–`C4` |
| **Excluded** | `tables/` (2 files each), `figures/` (388 KB each) |

`backward` used `OXYMODEL_N0_SCALE=67.056258`, the *Pseudomonas*
backward/forward ratio measured in the bacterial dataset
(`D2_analysis/C0_pseudomonas_scale.txt`). `none` is a **lower bound**, retained
for contrast and not a candidate estimate.

**Underpins:** that growth *T*<sub>opt</sub> (33.98 °C) and respiration
*T*<sub>opt</sub> (37.53 °C) are identical to eight significant figures across
all three routes and the CUE optimum moves only 0.07 °C — so the published
*Pseudomonas* temperature conclusion does not depend on the *N*₀ route.
`D2_forward/` also serves as the regression check that the `OXYMODEL_N0_ROUTE`
switch changed no default: its two tables are byte-identical to
`D1_baseline/tables/`.

---

## Regenerating any of these

See `reports/RUNBOOK.md` §3–§6. The full sequence takes about four minutes after
the first install.
