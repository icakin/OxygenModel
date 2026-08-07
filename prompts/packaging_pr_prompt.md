# Claude Code prompt — Packaging PR: reapply the environment and execution changes onto Ilgaz's main (autonomous)

Run from the project root (`.../work/OxygenModel`). Packaging and execution only. It changes NO
analysis decision, NO constant, NO script logic and NO number. The purpose is a small, science-free
pull request the repository owner has already agreed to merge.

NOTE TO USER: launch in an auto-approving mode. Needs network and credentials for
`git@github.com:icakin/OxygenModel.git`. Fast — this is a re-application of existing work, not new work.

CONTEXT. `gyd/d1-d2-review` carries the D1/D2 work, of which the packaging half has been reviewed and
accepted: *"the renv pin, the headless guards, adding stage 10 to the runner, and the run_all.sh
logging are improvements I'm keeping… send the science-free PR and I'll merge it straight away."*

DO NOT CHERRY-PICK. The two branches diverged at `ed5e1a6` and he has since renumbered every script.
Our changes were written against a layout where temperature-CUE was `10_temperature_cue.R`; on his
main it is `11_temperature_cue.R`, `08` is now `08_window_sensitivity.R`, and `12_joint_rK_estimator.R`
is new. A cherry-pick will either conflict or silently drop a stage. Reapply by hand against his head.

His current head is `893e6d3` (fetch first; it may have moved). His script layout:

    01_longdata.R  02_trimming.R  03_trim_selector.R  04_experiment_inputs.R
    05_oxygen_fits.R  06_main_figures.R  07_cutoff_sensitivity.R
    08_window_sensitivity.R  09_montecarlo_N0.R  10_simulation_recovery.R
    11_temperature_cue.R  12_joint_rK_estimator.R  config.R  run_all.R

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first: on `gyd/d1-d2-review` —
`SETUP.md`, `env/versions.json`, `renv.lock`, `scripts/00_install.R`, `scripts/run_all.sh`, the
headless guards in `scripts/03_trim_selector.R` and `scripts/04_experiment_inputs.R`, and `.gitignore`.
Then on `origin/main` — `scripts/run_all.R`, `scripts/config.R`, and the two app scripts as they now
stand. Understand both sides before writing anything.

PART A - branch
- `git fetch origin`, then branch off his CURRENT head: `git switch -c gyd/packaging origin/main`.
- Confirm you branched off his head and not off our review branch. Print the base commit.

PART B - reapply, adapted to his layout
Bring across ONLY these, rewritten as needed for his numbering:
- `renv.lock` and the renv activation (`.Rprofile`, `renv/activate.R`) — the pinned library. Verify
  `renv::status()` is sane and that the lockfile covers every package his 01-12 actually load,
  INCLUDING anything new since we branched (`rstan` for 12, `lme4`/`lmerTest`/`multcomp` for 06,
  whatever 08 needs). Add what is missing; do not remove anything.
- `env/versions.json` — regenerate on this machine rather than copying ours; it must describe the
  environment as it is now.
- `scripts/00_install.R` — idempotent installer. Update it for the packages above.
- `scripts/run_all.sh` — `set -euo pipefail`, per-stage timing, tee'd log to `logs/`. Must call HIS
  script names and numbers, and must include `11_temperature_cue.R`, which is the stage that was
  missing from the runner.
- Headless guards on `03_trim_selector.R` and `04_experiment_inputs.R` so they are sourceable without
  launching a browser. Reapply the guard to HIS versions of those files; do not overwrite them with
  ours. Their committed outputs remain inputs and must not be regenerated.
- `SETUP.md` — prerequisites, the install command, expected runtime per stage, disk. Rewrite the stage
  list for his numbering.
- `.gitignore` — `logs/`, `.DS_Store`, `Rplots.pdf`, LaTeX build artefacts (`*.aux`, `*.log`, `*.toc`,
  `*.tex` under `reports/` only — do NOT ignore `.log` globally, it will swallow real logs), and the
  renv machine-specific subdirectories.

EXPLICITLY EXCLUDED — none of these may appear in the diff:
  reports/ · prompts/ · runs/ · HANDOVER.md · any results/ file · any data/ file ·
  any change to config.R constants · any change to model, fitting or figure logic ·
  scripts/original_scripts/ (he has restored it with a provenance README; leave it alone)

PART C - verify it actually runs
- `bash scripts/run_all.sh` end to end on this branch, unattended. Report total wall time and
  per-stage timings. All twelve stages must run or be documented as intentionally skipped.
- CRITICAL: confirm the run reproduces HIS committed `results/` — this PR must change no number.
  Compare every regenerated table against the committed copy and report max absolute and max relative
  difference per file. If ANY table moves beyond floating-point noise, STOP and report it rather than
  committing; it means a packaging change altered behaviour, which is exactly what this PR must not do.
- Note for the report, do not fix: `N0_METHOD = "depletion"` uses `FC_TO_CELLS_PER_L = 10255100` while
  `data/Ninoc.csv` still encodes the old placeholder (implied constant 913,120), so the two N0 routes
  currently differ by 11.23x from the constant alone. That is his to regenerate and is out of scope here.

PART D - the pull request
- Push `gyd/packaging` to `origin`. Open a PR against `origin/main` with `gh` if available; otherwise
  print the URL, title and body to paste.
- Title: "Packaging: pinned environment, headless guards, complete runner"
- Body must state, briefly: what is included; that it changes no analysis decision, constant or
  number; the verification that his committed results regenerate unchanged; the total runtime; and
  that it is deliberately separate from the N0 work.
- Do NOT push to `origin/main`. Do NOT force-push. Verify `origin/main` is untouched at the end.

VERIFY (report all)
1. Base commit branched from, and confirmation it is his head not ours.
2. Every file added or modified, with a one-line reason each; confirmation the excluded list is absent
   from `git diff origin/main...gyd/packaging --stat`.
3. `bash scripts/run_all.sh` completes unattended; total and per-stage wall time; all twelve stages
   accounted for.
4. His committed `results/` regenerate unchanged — per-file max absolute and relative difference.
5. `renv.lock` covers every package 01-12 load, including rstan/lme4/lmerTest/multcomp.
6. PR opened (or URL/title/body printed); `origin/main` unchanged, verified.

CONSTRAINTS
- No analysis. No constant changed. No number changed. If making the pipeline run requires changing a
  number, STOP and report it — that is a finding, not a fix to make here.
- Reapply, do not cherry-pick. Adapt to his numbering rather than reverting it.
- Do not touch `results/`, `data/`, `scripts/original_scripts/`, or any of the excluded paths.
- Keep the diff small enough to review in ten minutes. This PR's value is that it is obviously safe.
- Autonomous; commit in parts: "packaging: pinned renv environment + 00_install + SETUP",
  "packaging: headless guards on the two Shiny scripts",
  "packaging: run_all.sh with per-stage timing and logging; include 11_temperature_cue",
  "packaging: gitignore hygiene".
```
