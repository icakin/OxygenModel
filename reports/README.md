# reports/

Reports for the OxygenModel D-series.

## The rule

**Every prompt that performs a SUBSTANTIVE ANALYSIS writes a self-contained
Quarto report rendered to PDF, in its own `reports/<ID>_<name>/` folder, with
the `.qmd` committed alongside the PDF so it is regenerable.**

**Prompts that only produce reference material — runbooks, scoping notes,
environment records — write plain markdown instead. A rendered PDF is for
findings, not for documentation.**

Concretely:

| Kind of output | Format | Lives in |
|---|---|---|
| Findings: comparisons, results, verdicts, anything with numbers that someone might cite | `.qmd` **and** rendered `.pdf`, both committed | `reports/<ID>_<name>/` |
| Reference: how to run something, what a schema is, what we plan to do | plain `.md` | `reports/` |

A findings report is *self-contained*: readable by someone who has not followed
the work. It must open with a **verdict paragraph**, carry **actual numbers**
(never "matches" without a figure), list its **blockers** explicitly, and close
with a **provenance footer** — the commit, the environment from
`env/versions.json`, the run date, and a table mapping every number in the
report to the file it came from.

Reports read artefacts written by scripts. They do not recompute analyses
inside the `.qmd`; if a number needs computing, a script computes it, writes it
to a CSV, and the report renders that CSV. This keeps the PDF a *view* of the
run rather than a second, divergent implementation of it.

## Layout

```
reports/
├── _shared/                  style assets used by every report
│   ├── header.tex            LaTeX preamble (tables, verdict box, floats)
│   ├── references.bib        shared bibliography
│   └── _template.qmd         skeleton: front matter, section order, provenance
├── <ID>_<name>/
│   ├── <ID>_<name>.qmd       committed
│   ├── <ID>_<name>.pdf       committed, rendered
│   ├── figures/              figures the report generates itself
│   └── scripts/              scripts that produce this report's artefacts
├── RUNBOOK.md                reference: from-scratch instructions
├── MULTIDATASET_SCOPING.md   reference: multi-dataset scoping
└── README.md                 this file
```

## Starting a new report

```bash
mkdir -p reports/D2_example/figures reports/D2_example/scripts
cp reports/_shared/_template.qmd reports/D2_example/D2_example.qmd
# fill it in, then:
quarto render reports/D2_example/D2_example.qmd
git add reports/D2_example/D2_example.qmd reports/D2_example/D2_example.pdf
```

Rendering needs Quarto and a LaTeX engine; `Rscript scripts/00_install.R`
verifies both and fails loudly if either is missing.

## Index

| ID | Title | Date | Concludes |
|---|---|---|---|
| D1 | [Baseline reproduction](D1_baseline/D1_baseline.pdf) | 2026-07-26 | Yes, a third party can regenerate the published results from the deposited data — the modular pipeline reproduces 19 of the 22 regenerable committed tables as IDENTICAL, 10 of them at exactly 0.0. Two caveats: Figure 6 (6 tables, 3 text files, 5 TIFFs) is produced by no script in the repository, and `results/tables/` turns out to have come from the modular pipeline, not from `scripts/original_scripts/`. |
| D2 | [N₀ reconciliation and dual route](D2_n0_routes/D2_n0_routes.pdf) | 2026-07-26 | The published *Pseudomonas* temperature result does **not** need correcting: growth and respiration T<sub>opt</sub> are structurally independent of N₀ and the CUE optimum moves 0.07 °C. But the two independent measurements of cell density agree with each other and both sit 49–66× above `N_inoculation_cells_per_L`, so absolute per-cell respiration is high by a median factor of 57, and the provenance of `N_inoc` could not be recovered from the deposit. ⚠ **See "Subsequent developments" below.** |
| D5 | [Estimator uncertainty](D5_estimator_uncertainty/D5_estimator_uncertainty.pdf) | 2026-08-07 | The two-stage pipeline is **not** over-confident. The proposed mechanism — autocorrelated residuals inflating confidence ~10× — is absent: median AR(1) ρ = −0.03, and 0 of 75 curves are in the hypothesised regime, because `05` fits the raw oxygen series, not the spline. The *r*–*K* ridge is real (−0.985) but points along the direction per-cell *R* is insensitive to, so discarding it is conservative. What limits the design is replicate variation, 37× the within-curve SE; inflating the measurement covariance 100× changes nothing. The measured case is for **misspecification**: a 15-min lag biases *R* by +19%, a saturating late phase by −78%. D6 should be rescoped around that, not covariance propagation. |
| D6 | [Model structure](D6_model_structure/D6_model_structure.pdf) | 2026-08-07 | Neither misspecification D5 flagged is present in the real curves. A nested F-test rejects growth saturation and O2 limitation at chance level (3/75 each); every extension is *worse* by AIC by almost exactly the per-parameter penalty, and the extra parameters sit at the boundary where the model collapses back to M0. Not a power problem — the same models recover a planted 15-min lag in 100% of simulated curves at exactly 15.00 min. The reason is the trimming: the median window still holds 44.5% of its oxygen at the fit end. The window sensitivity that motivated the hypothesis is **arithmetic**, tracking fitted *r* at R²=0.977, and no model term reduces it. Between-taxon ordering identical (Spearman 1.000) under all six models. **No change to the default estimator follows.** |

## Subsequent developments

Reports are a record of what was found **at the time**. Where the repository has
moved since, it is noted here rather than by editing the report.

### D2 — N₀ reconciliation and dual route · noted 2026-08-07

D2 was written against `ed5e1a6`. Three things changed on `main` afterwards.

**1. The pipeline adopted D2's backward route as the default.** `dc0951c`
(2026-08-06) introduced depletion-anchored N₀ and `config.R` now ships
`N0_METHOD <- "depletion"`, i.e.

```
N0 = FC_Final · FC_TO_CELLS_PER_L · exp(−r · (t_depletion − fit_start))
```

That is the route D2 argued for — it assumes only what the O₂ model already
assumes. The forward route (`N_inoc · exp(r·δ)`) survives as the `"initial"`
fallback. D2's recommendation is therefore **implemented**, not pending.

**2. `FC_TO_CELLS_PER_L` was re-derived.** `1cb16fb` (2026-08-07) replaced
`909916` — a placeholder tuned to a target median R, implying an impossible
sub-1× dilution — with `10255100`, derived from the stated dilution chain
`(500/490)·(500/50)·(502.5/500) × 1e6`. This sets the **absolute** respiration
and CUE scale only; slopes and relative ordering are unaffected. **D2's absolute
per-cell respiration numbers predate this change and should not be quoted as
current.** Its relative conclusions are unaffected.

**3. `data/Ninoc.csv` was rewritten — but its values did not materially move.**
`5ed17b0` (2026-08-03, "Correct Ninoc (drop 3x clamp)") rewrites all 75 rows,
so the file differs textually from the one D2 analysed. The numbers do not:

| | D2 baseline (`ed5e1a6`) | current `main` |
|---|---|---|
| median `N_inoculation_cells_per_L` | 24,605,191 cells/L | 24,788,014 cells/L |
| median `N_inoc` / `FC_Initial` | 155,299 | 155,232 |

Median row-wise ratio new/old = **1.000×**. So **D2's central discrepancy is not
resolved by that rewrite and stands as current**: the flow-cytometry and OD
measurements still sit ~49–66× above `N_inoc`, and the provenance of `N_inoc`
is still unrecovered. Anyone reading D2 should treat that finding as live.

> Recomputing D2's headline numbers under the new constant and default route is
> out of scope for a documentation change and has **not** been done here. It is
> the obvious next check.

## Re-rendering these reports

**Both reports rebuild from a clean clone.** Verified 2026-08-07 by cloning the
branch into a scratch directory and rendering each one there:

```bash
quarto render reports/D1_baseline/D1_baseline.qmd  --to pdf   # 18 pages
quarto render reports/D2_n0_routes/D2_n0_routes.qmd --to pdf   # 16 pages
```

### What is tracked under `runs/`, and what is not

The 48 small text tables the two `.qmd` files actually read — **133 KB**, in
`runs/D1_baseline/comparisons/` and `runs/D2_analysis/`. Determined by parsing
every read call in both sources, not by taking whole directories.

Everything else under `runs/` is excluded: the figures and `.rds` model objects
(~47 MB), `runs/D1_baseline/tables/`, `runs/D1_originals/`, and the three D2
route trees. `runs/MANIFEST.md` records the breakdown per tree.

Consequence: the **reports** rebuild, but the **comparison scripts** in each
report's `scripts/` folder do not — they consume the excluded run outputs.
Re-deriving the tables from scratch means re-running the pipeline first
(`RUNBOOK.md` §5–§6). Rendering the reports needs none of that.

`env/versions.json` feeds the provenance footer only; it arrives with the
packaging change. Both reports read it defensively, so a clone without it still
renders and that footer reads "not recorded".

### Known: the committed D2 PDF is one sentence behind its source

A fresh render of `D2_n0_routes.qmd` contains a sentence the committed PDF does
not — *"Only CUE can move, and it moves by 0.07 °C (backward) and 0.13 °C (lower
bound)."* Both the `.qmd` and the `.pdf` were committed together in `f914b88`,
so the PDF was rendered just before the final edit to the source and never
re-rendered.

**No number is contradicted.** All 752 substantive numbers in the report match
between the committed and freshly rendered PDFs, and both 0.07 °C and 0.13 °C
already appear elsewhere in the committed PDF; the new sentence restates them.
The committed PDF has deliberately **not** been overwritten — it is the record
of what was published at the time. Re-render it whenever D2 is next revisited.
