# Prompts — OxygenModel

Each file here is a prompt to paste into **Claude Code** (run from the project root
`.../work/OxygenModel`). Most are autonomous and modify code, run analyses, and commit;
launch Claude Code in an auto-approving mode (accept edits + allow commands, e.g.
`--dangerously-skip-permissions`) to run them unattended. They build on each other, so
run them in the order below.

Convention follows the `etcGEMs` and `Candidas` repos: the **`D` series** is the
oxygen-dynamics estimator — reproducing the published method, then upgrading it.

Status legend: ✅ run · ▶ in progress · ⏳ pending

## What this repo is

The reference implementation of the method published as **Cakin et al. 2026,
*ISME Communications* 6:ycag024** — simultaneous growth and respiration from a single
dissolved-oxygen time series. It is the **engine**; applications (currently `Candidas`)
should depend on a pinned version of it rather than re-implementing the estimator.

Three decisions frame the D-series:

- **Multi-dataset from the start.** The engine should ingest both the bacterial series
  here and *Candida*-style temperature series, so the maintenance/yield split has more
  than one taxon to identify from.
- **Curation provenance is fixed in the engine**, not per-application.
- **The Bayesian upgrade prototypes in `brms` and ships as hand-written Stan.**

## Reporting convention

Every prompt that performs a **substantive analysis** writes a self-contained Quarto
report rendered to PDF, in its own folder under `reports/`:

```
reports/
├── _shared/           header/theme, references.bib, _template.qmd
├── D1_baseline/       D1_baseline.qmd  +  D1_baseline.pdf  +  figures/
├── D2_n0_routes/      ...
└── README.md          one line per report: ID, title, date, conclusion
```

The `.qmd` is committed beside the PDF so every report is regenerable, and each report
carries a provenance footer — commit, environment, run date, and a table mapping every
number to the file it came from.

Prompts that produce only reference material — runbooks, scoping notes, environment
records — write plain markdown instead. A rendered PDF is for findings, not documentation.

D1 sets up `reports/_shared/` and the template; later prompts copy it.

## Run order

**Baseline**

1. `D1_baseline_reproduction_prompt.md` — ✅ pin the environment, make the pipeline run
   unattended (add the temperature-CUE stage to the runner; guard the two Shiny apps),
   regenerate everything into `runs/D1_baseline/`, and run three comparisons: modular
   pipeline vs committed outputs, modular vs `scripts/original_scripts/` (the code that
   produced the paper — the restructure at `ed5e1a6` has never been checked against it),
   and the paper's headline numbers (R² vs OD600 and flow cytometry, Bland–Altman bias,
   the *Pseudomonas* thermal parameters, where the CUE optimum sits). Also audits the
   data, reports whether the QC gate or the manual curation is doing the work, and
   scopes what stands in the way of multi-dataset ingestion. Changes no analysis
   decision. **Run first.**

   > **Outcome** — the pipeline reproduces from the deposited data in 63 s, unattended;
   > 19 of 22 regenerable committed tables IDENTICAL, 10 at exactly 0.0.
   > But `results/tables/` came from the **modular** pipeline, not `original_scripts/`,
   > and **Figure 6 is produced by no script in the repository**.
   > → [`reports/D1_baseline/D1_baseline.pdf`](../reports/D1_baseline/D1_baseline.pdf)

**The N₀ question**

2. `D2_n0_reconciliation_and_dual_route_prompt.md` — ✅ reconcile the biomass measurements,
   then test both back-calculations. Three independent statements of cell density currently
   disagree: the paper's stated OD₆₀₀ = 0.0005 at aliquoting, the per-replicate flow-cytometry
   counts, and `N_inoculation_cells_per_L` — whose provenance D1 could not recover (median
   ratio 6.6 against FC, and not a simple OD conversion). Once they're on one scale:

   ```
   forward   N₀ = N_inoc · exp(r·δ)              assumes growth at rate r through stabilisation
   backward  N₀ = N_end  · exp(−r·(t_end − t_0))  assumes only what the O₂ model already assumes
   ```

   `FC_Final` is a genuine per-replicate endpoint count, so their ratio **measures** growth
   during stabilisation instead of bounding it. Then re-run the *Pseudomonas* CUE analysis
   under both routes and quantify how far the published result moves — the formal assessment
   that decides whether the journal needs contacting. **Run after D1.**

   > **Outcome** — flow cytometry and OD agree with each other to 1.26× but sit 66× and 49×
   > above `N_inoc`, whose provenance could not be recovered; backward/forward N₀ differs by
   > a median factor of **57**.
   > The *Pseudomonas* temperature **conclusion holds** under all three routes (growth and
   > respiration optima identical to 8 s.f., CUE optimum moves 0.07 °C) — what moves is the
   > **level** (peak CUE 0.209 → 0.946). No journal correction warranted; the absolute
   > respiration scale is unresolved.
   > → [`reports/D2_n0_routes/D2_n0_routes.pdf`](../reports/D2_n0_routes/D2_n0_routes.pdf)

**Structure and provenance**

> **No prompt file exists yet for D3 or D4.** They are specified below but not
> written; there is nothing in this folder to run for either. D5 is the next
> prompt with a file.

3. `D3` — ⏳ *(not yet written)* the data contract: generalise `config.R` and the numbered
   stages so the engine ingests both the bacterial series and a *Candida*-style
   temperature series, per the scoping D1 produces.
4. `D4` — ⏳ *(not yet written)* curation provenance: instrument `03_trim_selector.R`
   (~1,500 lines) with a decision log — operator, timestamp, action, reason code from a
   fixed vocabulary — make the displayed window identical to the fitted one, and measure
   inter-operator agreement on a re-trimmed subset. Fix it here so applications inherit it.

**The Bayesian upgrade**

5. `D5_estimator_uncertainty_prompt.md` — ⏳ **next to run.** v0: faithful port of the
   current estimator into a testable form; regression baseline; the validation suite
   assembled from the cutoff-sensitivity, Monte-Carlo-N₀ and simulation-recovery stages,
   extended to cover the *r*–*K* joint behaviour and misspecification (lag, drift,
   saturation) rather than only recovery under the model it fits.

   > Note: `12_joint_rK_estimator.R` already landed on `main` independently and
   > propagates the within-curve *r*–*K* covariance through a hierarchical
   > measurement-error model in Stan. D5 should build on it rather than duplicate it.
6. `D6` — ⏳ v1: one-stage hierarchical fit to the traces in `brms` — correlated random
   effects, AR(1) residuals. Propagates the *r*–*K* covariance; ends the measurement-error
   understatement; non-respiring wells become `r ≈ 0` rather than deletions.
7. `D7` — ⏳ v2: explicit lag/equilibration parameter in hand-written Stan, so N₀ is
   estimated rather than assumed.
8. `D8` — ⏳ v3: reparameterise for the ridge (rotated basis, non-centred hierarchy).
9. `D9` — ⏳ v4: maintenance vs growth-associated respiration, `q = m + Y⁻¹r` — the split
   the paper's introduction describes but the model does not implement.

**Housekeeping** — not analyses; they move files and open pull requests. They change no
number and write no report.

- `handover_tidy_and_push_prompt.md` — ✅ tidy the working tree, write the run manifest
  and gitignore policy, push the review branch.
- `packaging_pr_prompt.md` — reapply the environment and execution changes (pinned
  `renv.lock`, headless guards, complete `run_all.sh`) onto the upstream `main` as a
  science-free PR.
- `prompts_and_reports_pr_prompt.md` — this folder and `reports/` onto `main`, so the
  specification record stops being branch-dependent.

**Blocked**

- Temperature compensation. `T_internal [°C]` is recorded by the PreSens exports and
  stripped at conversion; in the *Candida* plates the sample runs up to 4.8 K below the
  set point the optode compensates at, over exactly the pre-window interval. Needs the
  original `.xlsx` exports (requested).
- ~~Packaging. Deliberately last — no point refactoring a moving target.~~ Done: see
  `packaging_pr_prompt.md` above. The pinned environment, the headless guards and the
  complete runner are now a separate science-free PR.

## The through-line

The method infers growth and respiration simultaneously from one oxygen trace, and its
growth estimates are externally validated against flow cytometry and OD. What is not yet
validated is the respiration side, which depends on a cell density that cannot be measured
in a sealed vial and is therefore back-calculated. The D-series makes the estimator
reproducible, settles which back-calculation the data support, and replaces the two-stage
point-estimate pipeline with one hierarchical model that propagates its own uncertainty.
