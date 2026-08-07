# Claude Code prompt — Amend PR #2 so the reports actually re-render (autonomous)

Run from the project root (`.../work/OxygenModel`). A small amendment to an OPEN pull request. It
runs no analysis, re-renders nothing by default, and changes no result. It adds the input tables the
two committed reports need in order to rebuild from a clean clone.

NOTE TO USER: launch in an auto-approving mode. Small — file additions and a commit onto an existing
branch. PR #2 (`gyd/docs`) is already open; this pushes to it rather than creating anything new.

CONTEXT. PR #2 brought `prompts/` and `reports/` onto a branch off `main`, and deliberately excluded
`runs/` — 58 MB of non-canonical output trees. That exclusion was right for the figures and the
`.rds` model objects. It was wrong for the tables: both `reports/D1_baseline/D1_baseline.qmd` and
`reports/D2_n0_routes/D2_n0_routes.qmd` read from `runs/`, so as things stand the `.qmd` files are
committed but cannot rebuild, which defeats the point of committing them.

The convention says every report is regenerable. Right now that is aspirational rather than true.
This prompt makes it true, at the smallest cost that achieves it.

---

```
Work AUTONOMOUSLY end to end; commit onto the EXISTING branch `gyd/docs`; print a summary. Read
first: both report `.qmd` sources and every path they read; `runs/MANIFEST.md`; the current
`.gitignore`; and `reports/README.md`. Do NOT re-run any analysis, refit anything, or alter any number.

PART A - establish exactly what the reports need
- Parse both `.qmd` files for every file path they read. List them, with size, and mark each as
  present or absent in the repository as PR #2 currently stands.
- Do NOT guess from directory names. If a report reads a whole directory, enumerate what it actually
  consumes at render time.
- Report the total size of the required set.

PART B - bring across the tables, and only the tables
- Add the required TABLES from `runs/` (CSV, TSV, TXT — the small text outputs the reports read).
- Do NOT add: `.rds` model objects, `.tiff`/`.png`/`.pdf` figure files, or any table a report does
  not actually read. If a required file is large (say >5 MB), flag it and ask rather than committing
  it silently.
- Preserve the `runs/<tree>/tables/...` paths the reports already reference, so no `.qmd` edit is
  needed. If a path must change, change the `.qmd` rather than the layout, and say so.
- Update `.gitignore` to express the policy precisely: `runs/` excluded EXCEPT the table
  subdirectories, with a one-line comment stating why. Keep the existing
  `!reports/_shared/header.tex` negation and do not introduce a global `*.log` rule.
- `runs/MANIFEST.md` must remain and must now also record, per tree, which of its files are tracked
  and which are excluded, so the omission is documented rather than implicit.

PART C - prove it
- From a CLEAN CLONE of the branch in a scratch directory (not the working tree), render both reports
  and confirm they build. This is the whole point of the amendment — a render that only works because
  of untracked files in your working directory does not count.
- Compare the freshly rendered PDFs against the committed ones. They should be materially identical;
  small differences from timestamps or font embedding are fine and should be reported as such. If any
  NUMBER differs, STOP and report it — that would mean the committed PDF and its `.qmd` disagree.
- Do not replace the committed PDFs with the fresh renders unless a number differs, in which case
  report rather than overwrite.

PART D - update the record
- `reports/README.md`: state plainly that both reports regenerate from a clean clone, and what is
  and is not tracked under `runs/`.
- Add the same two lines to the PR #2 description via `gh pr edit`, so a reviewer sees why the diff
  grew.

PART E - push
- Commit onto `gyd/docs` and push. Do NOT open a new PR — PR #2 updates automatically.
- Do NOT push to `main`. Verify `origin/main` is unchanged.

VERIFY (report all)
1. The full list of paths both reports read, with sizes, and which were missing.
2. What was added, with total size added to the repository; confirmation no `.rds` and no figure
   binaries came across.
3. The `.gitignore` rule as written, and confirmation the `_shared/header.tex` negation survives and
   no global `*.log` rule exists.
4. Clean-clone render: both reports build; page counts; any difference from the committed PDFs, and
   whether any NUMBER differs.
5. `runs/MANIFEST.md` updated with the tracked/excluded breakdown per tree.
6. PR #2 description updated; `origin/main` unchanged.

CONSTRAINTS
- No analysis, no refits, no number changes. If the clean-clone render produces a different number
  from the committed PDF, that is a finding to report, not to fix by overwriting.
- Tables only. No `.rds`, no figure binaries, nothing a report does not read.
- Amend the existing branch and PR. Do not open a new one, do not force-push, do not touch `main`.
- Keep the size increase small and state it. If it exceeds ~10 MB, stop and report rather than commit.
- Autonomous; commit in parts: "docs: add the runs/ tables the reports read",
  "docs: scope the runs/ gitignore to exclude bulk but keep tables",
  "docs: verify both reports render from a clean clone".
```
