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
| D2 | [N₀ reconciliation and dual route](D2_n0_routes/D2_n0_routes.pdf) | 2026-07-26 | The published *Pseudomonas* temperature result does **not** need correcting: growth and respiration T<sub>opt</sub> are structurally independent of N₀ and the CUE optimum moves 0.07 °C. But the two independent measurements of cell density agree with each other and both sit 49–66× above `N_inoculation_cells_per_L`, so absolute per-cell respiration is high by a median factor of 57, and the provenance of `N_inoc` could not be recovered from the deposit. |
