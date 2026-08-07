# Claude Code prompt — Bring prompts/ and reports/ onto main (autonomous)

Run from the project root (`.../work/OxygenModel`). Documentation only. It runs no analysis, refits
nothing, renders nothing new and changes no result. It puts the specification record and the two
rendered reports onto `main` so they stop being branch-dependent.

NOTE TO USER: launch in an auto-approving mode. Run AFTER PR #1 (packaging) has merged, so this
branches off a main that already has the pinned environment. Small and fast — file moves and a PR.

CONTEXT. `prompts/` and `reports/` are currently tracked only on `gyd/d1-d2-review`. The packaging
branch deliberately excluded them, so on any other checkout they vanish — which is a poor property
for the documents that record how every analysis in this repository was commissioned and what it
concluded. `reports/_shared/` in particular is the report template and style that D5 onward are
specified to follow, so it needs to be on `main` before the next analysis prompt runs.

What is on `gyd/d1-d2-review`:

    prompts/D1_baseline_reproduction_prompt.md
    prompts/D2_n0_reconciliation_and_dual_route_prompt.md
    prompts/README.md
    prompts/handover_tidy_and_push_prompt.md
    reports/D1_baseline/            (D1_baseline.qmd, .pdf, figures/, scripts/)
    reports/D2_n0_routes/           (D2_n0_routes.qmd, .pdf, figures/, scripts/)
    reports/_shared/                (header.tex, references.bib, _template.qmd)
    reports/README.md
    reports/RUNBOOK.md
    reports/MULTIDATASET_SCOPING.md

Also present on disk but UNTRACKED on the current branch, and to be included:

    prompts/D5_estimator_uncertainty_prompt.md
    prompts/packaging_pr_prompt.md
    prompts/prompts_and_reports_pr_prompt.md   (this file)

DELIBERATELY EXCLUDED: `runs/` — 58 MB of non-canonical output trees. `runs/MANIFEST.md` records
what each tree was and what claim it underpins, and that manifest travels with `reports/`; the trees
themselves do not belong in the repository. Do not bring them across in any form.

THIS IS A SHARED REPOSITORY. Push a branch and open a PR against `main`. Do not push to `main`.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first: on `gyd/d1-d2-review` —
`prompts/README.md`, `reports/README.md`, `reports/_shared/`, and both rendered reports; on
`origin/main` — the current `README.md` and repository layout. Do NOT re-run any analysis, re-render
any report, or alter any number.

PART A - branch
- `git fetch origin`. Confirm PR #1 has merged by checking that `SETUP.md`, `env/versions.json`,
  `renv.lock`, `scripts/00_install.R` and `scripts/run_all.sh` are present on `origin/main`. If they
  are NOT, STOP and report it — this prompt is meant to run after the packaging lands.
- `git switch -c gyd/docs origin/main`. Print the base commit.

PART B - bring the files across
- Take `prompts/` and `reports/` from `gyd/d1-d2-review` (`git checkout gyd/d1-d2-review -- prompts
  reports`, or equivalent), plus the three untracked prompt files currently on disk.
- Verify `runs/` is absent and stays absent. Verify no `results/`, `data/` or `scripts/` file is
  touched by this branch — `git diff --stat origin/main...gyd/docs` should show only `prompts/`,
  `reports/`, `.gitignore` and `README.md`.
- Remove any LaTeX build artefacts that came along (`*.aux`, `*.log`, `*.toc`, `*.tex` under
  `reports/`) and confirm `.gitignore` covers them. Do NOT add a global `*.log` rule — it would
  swallow the pipeline logs `run_all.sh` writes.
- Confirm both PDFs are present and open, and that each `.qmd` is committed beside its PDF.

PART C - make the documents true on main
The prompts and reports were written before Ilgaz's fifteen commits and before the packaging merge.
Correct anything now stale, WITHOUT rewriting history or re-running analysis:
- `prompts/README.md`: mark D1 and D2 ✅ with their outcomes; add D5 as the next pending item; add the
  two housekeeping prompts; make sure the reporting convention section is accurate and points at
  `reports/_shared/`. Note explicitly that D3 (data contract) and D4 (curation provenance) are
  described but not yet written.
- `reports/README.md`: one line per report — ID, title, date, conclusion.
- `reports/RUNBOOK.md`: update for the current script numbering (01-12, with 08 window sensitivity,
  09 montecarlo, 10 simulation recovery, 11 temperature CUE, 12 joint estimator) and for
  `run_all.sh` now existing on main.
- Add a short section to the root `README.md` pointing at `prompts/`, `reports/` and `SETUP.md`, and
  describing the convention in two lines. Do not rewrite the existing content.
- ADD A DATED NOTE, do not silently edit, wherever a report's conclusion has since been superseded.
  Specifically: `reports/D2_n0_routes/` concluded on a `Ninoc.csv` that has since been regenerated
  and on `FC_TO_CELLS_PER_L = 909916`, which was a placeholder and is now derived as 10,255,100.
  Add a clearly marked "Subsequent developments" note at the top of the report README entry — and, if
  it can be done without re-rendering, as a preamble in the .qmd — recording that the 24x residual
  was traced to a stale file and the constant has since been corrected. The report stands as the
  record of what was found at the time; the note stops a later reader taking its numbers as current.

PART D - PR
- Push `gyd/docs` and open a PR against `origin/main` with `gh`. Title: "Docs: prompt series and
  rendered reports". Body: what is included, what is deliberately excluded (`runs/`), that it changes
  no code and no number, and that `reports/_shared/` is the template the next analysis prompt uses.
- Do not push to `main`. Do not force-push. Verify `origin/main` is unchanged at the end.
- After the PR is open, delete the local `gyd/d1-d2-review` branch ONLY if its entire content is now
  either merged (packaging), included here (prompts, reports), or deliberately excluded (`runs/`).
  Confirm that inventory explicitly before deleting, and do not delete the remote copy if one exists.

VERIFY (report all)
1. Confirmation the packaging PR had merged before branching; the base commit used.
2. `git diff --stat origin/main...gyd/docs` — confirm only prompts/, reports/, .gitignore, README.md.
3. `runs/` absent; no results/, data/ or scripts/ file touched.
4. Both PDFs present with their .qmd beside them; LaTeX artefacts removed and gitignored; no global
   `*.log` rule.
5. The stale-content corrections made in PART C, listed individually.
6. The "Subsequent developments" note recording the Ninoc regeneration and the corrected constant.
7. PR opened; `origin/main` unchanged.
8. The `gyd/d1-d2-review` inventory, and whether it was deleted.

CONSTRAINTS
- No analysis, no refits, no re-rendering, no number changes. If bringing the reports across reveals
  a numerical inconsistency, record it in the "Subsequent developments" note — do not resolve it by
  recomputing.
- `runs/` does not come across, in any form.
- Never modify `results/`, `data/`, `scripts/` or `scripts/original_scripts/` on this branch.
- Reports are a record of what was found at the time. Correct stale *documentation*; annotate, never
  rewrite, superseded *findings*.
- Keep the diff reviewable in ten minutes. Its value is that it is obviously safe.
- Autonomous; commit in parts: "docs: bring prompts/ and reports/ onto main",
  "docs: refresh prompt and report indices for the current state",
  "docs: note subsequent developments on the D2 report", "docs: runbook for current numbering".
```
