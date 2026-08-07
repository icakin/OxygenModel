# Claude Code prompt — Handover: tidy the repository, write the handover note, and push for review (autonomous)

Run from the project root (`.../work/OxygenModel`). Housekeeping and handover only. It runs no
analysis, refits nothing, and changes no result. The purpose is to leave the repository in a state
where its author can pick the work up cold, and to get the D1/D2 work onto GitHub for review.

NOTE TO USER: launch in an auto-approving mode. Fast — file moves, a written note, and git. The push
step needs network and credentials for `git@github.com:icakin/OxygenModel`.

CONTEXT. D1 and D2 are complete, committed in eleven parts on local `main`, which is now 11 commits
ahead of `origin/main`. `prompts/` is untracked. The working tree has also accumulated `runs/` (58 MB
across `D1_baseline/` and the three D2 route trees), `logs/`, `Rplots.pdf` and `.DS_Store`. The two
rendered reports (`reports/D1_baseline/D1_baseline.pdf`, `reports/D2_n0_routes/D2_n0_routes.pdf`) are
committed and are the substantive record.

THIS IS SOMEONE ELSE'S REPOSITORY. The remote belongs to the first author of the published paper.
Push a BRANCH and open a pull request; do not push to `main`. The point is to hand work over for
review, not to land it.

The findings being handed over, for the handover note:
- D1: the pipeline reproduces from the deposited data in 63 s. But `results/tables/` was produced by
  the modular pipeline, not `scripts/original_scripts/`, and Figure 6 is produced by no script in the
  repository. The QC gate passes 75/75 because `manual_fit_windows.csv` sets both edges on all 75
  curves; window choice alone moves `K` and `R` by 23% and `C_tot` by 69%.
- D2: `N_inoculation_cells_per_L` cannot be reconciled with either independent measurement — flow
  cytometry and OD agree with each other to 1.26x but sit 66x and 49x above it, and in several taxa
  its largest value is exactly 3.000000x its smallest to 7-8 significant figures. Backward/forward N0
  differs by a median factor of 57, against a published Monte Carlo that spanned CV = 0.10-0.40. The
  *Pseudomonas* temperature CONCLUSION holds under all three routes (growth and respiration optima
  identical to 8 s.f., CUE optimum moves 0.07 C) because r is independent of N0 by construction; what
  moves is the LEVEL (peak CUE 0.209 -> 0.946).
- Unresolved: CUE 0.946 is not physiologically credible either, so the reconciliation is not closed —
  the carbon side (SYBR counting non-viable particles, 100 fg C um^-3 described in the paper as "the
  upper end", RQ = 1) needs the same treatment.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first: prompts/README.md,
reports/README.md, reports/{RUNBOOK.md,MULTIDATASET_SCOPING.md}, both rendered reports and their
.qmd sources, SETUP.md, .gitignore, and `git log --oneline origin/main..HEAD`. Do NOT re-run any
analysis, refit anything, or alter any number.

PART A - tidy the working tree
- Delete `Rplots.pdf` and any `.DS_Store`. Add both to `.gitignore`.
- `runs/` is 58 MB of non-canonical output trees that the two reports cite as evidence. Apply the
  same policy the Candidas repo uses: keep TABLES tracked (small CSVs, they are the evidence), and
  gitignore the bulk (figures, .rds). Write `runs/MANIFEST.md` recording, for each tree: which prompt
  produced it, which commit, when, the environment from `env/versions.json`, and what claim it
  underpins. State the policy in one line at the top of that .gitignore block.
- Check `renv/`: commit `renv.lock`, gitignore `renv/library/`, `renv/staging/`, `renv/sandbox/` and
  the other machine-specific subdirectories per renv's own recommendation. Verify a fresh
  `renv::restore()` would still work from what remains tracked.
- Gitignore `logs/` but keep the directory (a `.gitkeep`), since `run_all.sh` writes there.
- `git add prompts/` — it is currently untracked and it is the record of how this work was specified.
- NEVER modify `results/`, `data/`, `scripts/original_scripts/`, or either rendered report. Verify by
  checksum afterwards that all three are byte-identical to their pre-run state.

PART B - the handover note
Write `HANDOVER.md` at the project root. It is the single document someone picks this up from, so
write it for a reader returning cold, not as a changelog. It must contain:
- **Status**: what has run (D1, D2), what it concluded in three sentences each, and where the full
  detail lives (the two PDFs — give paths).
- **What is solid**: the pipeline reproduces in 63 s from deposited data; the growth-rate validation
  against OD600 and flow cytometry is untouched by any of this; the temperature conclusion holds
  under all three N0 routes. Say this FIRST and plainly.
- **What is open**, each with the specific next action:
  * `N_inoculation_cells_per_L` provenance — unrecovered; blocks the absolute respiration/CUE scale.
  * The parallel-vial 45-min reference density described in the methods — not in the deposit.
  * The carbon side — CUE 0.946 under the backward route is not credible, so the reconciliation is
    not closed; FC viability, the carbon-density choice and RQ all need checking.
  * `T_internal` — recorded by the PreSens exports, stripped at conversion; the original .xlsx
    exports are needed. Note that the published methods already describe this drift qualitatively.
  * Figure 6 — produced by no script in the repository.
- **Decisions already taken**, so they are not relitigated: multi-dataset engine; curation provenance
  fixed here rather than per-application; Bayesian prototypes in brms and ships as Stan; the
  back-calculation is a design choice forced by the sealed vial, so the question is which one.
- **How to resume**: the exact commands (`bash scripts/00_install.R` / `run_all.sh`), what to read
  first, and which prompt is next (`prompts/README.md` has the ordered series; D3 onward are
  described but not yet written).
- **The reporting convention**, in two lines, so the next report matches.
- Keep it under two pages. Link out rather than restate.

PART C - update the indices
- `prompts/README.md`: mark D1 and D2 ✅ with a two-line outcome each, exactly as the Candidas
  prompts README does. Leave D3-D9 ⏳.
- `reports/README.md`: confirm both reports have their one-line entry with ID, title, date and
  conclusion.
- `README.md` (project root): add a short section pointing at `HANDOVER.md`, `SETUP.md` and
  `reports/`, and describing the `runs/` convention. Do not rewrite the existing content.

PART D - push for review
- Create a branch off current `main` named `gyd/d1-d2-review` and push it to `origin`.
- Do NOT push to `origin/main`. Do NOT force-push anything.
- Open a pull request against `origin/main` if `gh` is available and authenticated; if it is not,
  print the exact URL and the PR title and body to paste. The PR body should be a condensed version
  of HANDOVER.md's status and open-questions sections, and must state plainly that it is for review
  rather than to be merged unexamined.
- Verify the push landed: report the branch name, the commit range, and the file count and size
  pushed. Confirm `origin/main` is unchanged.

VERIFY (report all)
1. Stray files removed and gitignored; `runs/` policy applied with sizes before and after; what is
   now tracked versus ignored.
2. `runs/MANIFEST.md` covers every tree with producer, commit, date, environment and the claim it
   underpins.
3. `renv.lock` tracked, machine-specific renv subdirectories ignored, and a restore would still work.
4. `prompts/` tracked.
5. `results/`, `data/`, `scripts/original_scripts/` and both rendered PDFs byte-identical to their
   pre-run state, verified by checksum.
6. `HANDOVER.md` exists and covers all six required sections; state its length.
7. `prompts/README.md`, `reports/README.md` and the root `README.md` updated.
8. Branch pushed, commit range, PR opened or the URL and body printed; `origin/main` untouched.

CONSTRAINTS
- No analysis, no refits, no number changes. If tidying reveals a numerical inconsistency, record it
  in HANDOVER.md under open questions — do not resolve it by recomputing.
- Move, don't copy-then-delete. List every deletion before it happens.
- Branch and PR only. `origin/main` must be untouched at the end, verified.
- HANDOVER.md must lead with what is solid before what is open. The reader is the person who built
  this and published the method; the note should be useful to them, not a list of faults.
- Autonomous; commit in parts: "handover: tidy working tree, runs/ manifest and gitignore policy",
  "handover: track prompts/", "handover: HANDOVER.md + index updates",
  "handover: push gyd/d1-d2-review for review".
```
