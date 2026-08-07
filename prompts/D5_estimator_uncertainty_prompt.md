# Claude Code prompt — D5: measure what the estimator's uncertainty actually is (autonomous)

Run from the project root (`.../work/OxygenModel`). Diagnostic only. It changes no pipeline number,
no constant and no analysis decision. Its purpose is to replace an assertion — that the two-stage
pipeline understates its own uncertainty — with a measured number, obtained in this repository's own
validation harness on the deposited data. It is the evidence base for the one-stage rebuild (D6).

NOTE TO USER: launch in an auto-approving mode. Work on a branch and open a PR; do not push to `main`.
The pipeline runs in 63 s, so the cost here is the bootstraps and the simulation grid, not the fits.

CONTEXT.

The pipeline estimates each curve with `nlsLM` on the normalised oxygen model
`O2*(t) = O2*_0 + (k/r)(1 - exp(r t))`, then carries point estimates (and, in
`12_joint_rK_estimator.R`, the `vcov()` covariance) into a taxon-level model. Three observations
motivate this prompt:

1. `12_joint_rK_estimator.R` returns `sd_logr_joint` = 0.0577 for ALL fifteen taxa (range
   0.05768-0.05774) and `sd_logK_joint` = 0.0936-0.0937 for all fifteen. Under
   `y_i ~ MVN(mu_t, Tau + S_i)`, if `Tau >> S_i` then `sd(mu_t) ~ sqrt(Tau/n)`, and every taxon has
   n = 5 — so the same number every time. Back-solving gives between-replicate
   `sd(log r) ~ 0.129`. The per-curve covariance being propagated is contributing nothing.

2. The likely cause is that `vcov(nls)` is computed on ~200 residuals from a smooth curve sampled at
   1-minute intervals. If those residuals are strongly autocorrelated, the independence assumption
   behind the naive standard error is badly violated and the SE is too small by roughly
   `sqrt((1+rho)/(1-rho))`. At rho = 0.98 that is a factor of ~10 — the order of the gap observed
   between `se_r/r` and replicate-to-replicate reproducibility.

3. The script's own header records that a direct one-stage fit will not converge: *"the smooth curves
   let r and K trade off... The posterior is a long ridge and chains never mix."* That diagnosis is
   correct. The remedy is reparameterisation plus honest residual structure, not two stages — but the
   case for that has to be measured, not argued.

D5 measures all three. It does NOT build the replacement model; that is D6.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first: scripts/config.R (the
model function, FIT_LOWER/FIT_UPPER, QC thresholds), scripts/05_oxygen_fits.R (the per-curve fit and
what it retains), scripts/09_montecarlo_N0.R, scripts/10_simulation_recovery.R,
scripts/12_joint_rK_estimator.R (both stages), scripts/08_window_sensitivity.R, and
results/tables/{oxygen_fit_results.csv, oxygen_fit_curves.csv, oxygen_results_with_R.csv,
joint_vs_twostage_uncertainty.csv}.

BRANCH: off `gyd/packaging`, not `main` — PR #1 is open but not yet merged, and this prompt needs the
pinned environment (renv, 00_install.R, run_all.sh) that only exists there.
`git switch -c gyd/d5-uncertainty gyd/packaging`. It will rebase cleanly once PR #1 and PR #2 land.

PART 0 - complete the runner
`run_all.R` on main covers 01, 02 and 05-10 only. The packaging branch adds 11. Script 12
(`12_joint_rK_estimator.R`) is in NEITHER, and this prompt re-runs it twice, so add it — as the last
stage, guarded so a missing `rstan` degrades to a clear skip message rather than aborting the run.
Note in the commit message that 03 and 04 remain deliberately outside the runner (they are the Shiny
apps, whose committed outputs are inputs). Verify `bash scripts/run_all.sh` still completes and still
reproduces the committed `results/` unchanged before going further.

PART A - residual autocorrelation, and what it does to vcov(nls)
- Refit every curve exactly as `05` does, and retain the residual series in fit-window time order.
- Per curve report: lag-1 residual autocorrelation, the Durbin-Watson statistic, an AR(1) rho fitted
  to the residuals, and the implied effective sample size `n_eff = n (1-rho)/(1+rho)`.
- Report the variance-inflation factor `sqrt((1+rho)/(1-rho))` per curve and its distribution overall
  and by taxon.
- Then compute, for each curve, BOTH the naive `vcov(nls)` standard errors and an
  autocorrelation-corrected version (either a sandwich/Newey-West estimator with a stated bandwidth,
  or refitting under an explicit AR(1) error model — state which and why). Report the ratio.
- State plainly: by what factor does the naive SE understate the corrected one, median and range.

PART B - the naive SE against actual reproducibility
- Per taxon, compare the median `se_r/r` from `summary(nls)` with `sd(log r)` across that taxon's
  five replicates; likewise for K. Report both, their ratio, and the overall median ratio.
- Do the same for the PART A corrected SEs. If the correction is the right explanation, the corrected
  SEs should land much closer to the replicate spread; report how much of the gap it closes and how
  much remains (the remainder is genuine biological/handling variation between replicates, and
  quantifying that split is a result in itself).

PART C - the r-K ridge, measured three ways
Report `corr(log r, log K)` under each, with the distribution across curves:
  (i)   the `vcov(nls)` correlation per curve, as the pipeline currently sees it;
  (ii)  a parametric bootstrap: simulate from each fitted curve at its own residual scale, ONCE with
        i.i.d. noise and ONCE with the AR(1) structure estimated in PART A, refit each replicate
        dataset, and take the sampling correlation. The difference between the two is the point;
  (iii) the observed within-(taxon) correlation across the five real replicates.
- Also report the direction of the ridge in the rotated basis: the sampling variance of
  `log r + log K` versus `log K - log r`. State which combination is well determined and which is
  not, and note that CUE and per-cell R depend on the poorly determined one.
- Keep the bootstrap replicate count modest (a few hundred) and state it; this is a diagnostic, not a
  precision exercise.

PART D - extend the existing validation suite (ADDITIVELY - do not change current behaviour or outputs)
- `10_simulation_recovery.R`: it currently simulates from the fitted model and reports bias and RMSE
  on r and K separately. Add (a) the JOINT recovery — the sampling covariance and correlation of the
  recovered `(log r, log K)`, not just marginals; and (b) MISSPECIFICATION arms, where the truth
  differs from the fitted model: an initial lag/equilibration phase, a slow linear optode drift, and
  saturating rather than exponential growth late in the window. Fit the standard model to each and
  report the bias in r, K and per-cell R. This tests robustness, which recovery-under-truth cannot.
- `09_montecarlo_N0.R`: it currently varies N0 MAGNITUDE over a lognormal CV grid. Add a FUNCTIONAL
  FORM arm — compare N0 under the depletion anchor against the initial-count route against a
  no-back-calculation lower bound, at a matched constant, and report the spread in per-cell R. The
  distinction matters: magnitude rescales R, form changes its dependence on r.
- Both must keep their existing outputs byte-identical; new results go to new files.

PART E - the consequence for the joint estimator
- Re-run `12_joint_rK_estimator.R` unchanged to reproduce its current output.
- Then re-run it with `S_i` replaced by the PART A autocorrelation-corrected covariance, writing to a
  NEW output file. Do not modify the committed script's default behaviour; add an option or a copy.
- Report whether the posterior SDs become taxon-specific, and by how much the taxon-level intervals
  widen. Report max R-hat for both runs.
- If they remain identical across taxa even after correction, that is the more important finding and
  must be reported as such — it would mean the between-replicate variance genuinely dominates and the
  covariance is irrelevant at this design, which changes what D6 should aim at.

PART F - what it means downstream
- Take the taxon-level contrasts the pipeline currently reports and recompute their intervals under
  the corrected uncertainty. State how many contrasts that are currently "resolved" remain so.
- One paragraph, written to be quoted: is the two-stage pipeline over-confident, by how much, and
  does the one-stage rebuild follow from the numbers or not.

PART G - report
- `reports/D5_estimator_uncertainty/D5_estimator_uncertainty.qmd` rendered to PDF, following the
  convention established in `reports/_shared/`: self-contained, readable by someone who has not
  followed this work, figures in its own `figures/`, and a provenance footer mapping every number to
  the file it came from.
- Include: the autocorrelation distribution; naive vs corrected SE; SE vs replicate reproducibility;
  the three ridge estimates and the rotated-basis variances; the misspecification biases; the joint
  model before and after; and the downstream contrast count.
- Be even-handed. If the correction turns out to be small, or the ridge less severe than expected,
  say so as plainly as the alternative. A null result here means the two-stage pipeline is fine and
  D6 should be rescoped, which is worth knowing.
- Add its one-line entry to `reports/README.md`.

PART H - PR
- Push `gyd/d5-uncertainty` and open a PR against `origin/main`. Body: what was measured, the headline
  numbers, and an explicit statement that no pipeline number changes. Do not push to `main`.

VERIFY (report all)
0. Script 12 added to `run_all.R`; `bash scripts/run_all.sh` completes and still reproduces the
   committed `results/` unchanged (per-file, no table beyond floating-point noise).
1. Residual lag-1 autocorrelation, AR(1) rho, n_eff and variance-inflation factor: distribution
   overall and by taxon.
2. Naive vs autocorrelation-corrected SE: median and range of the ratio, for r and for K.
3. Per taxon: median se_r/r, sd(log r) across replicates, and both ratios (naive and corrected);
   how much of the gap the correction closes.
4. corr(log r, log K) under all three methods; sampling variance of `log r + log K` vs `log K - log r`.
5. Joint recovery and the three misspecification biases from the extended `10`; the N0 functional-form
   spread from the extended `09`; confirmation both scripts' existing outputs are byte-identical.
6. Joint estimator before and after the corrected S_i: posterior SDs per taxon, interval widths, R-hat.
7. Contrast count before and after.
8. Report renders; page count; no missing figure or unresolved reference; .qmd committed beside it.
9. `results/`, `data/` and all existing outputs byte-identical, verified by checksum.

CONSTRAINTS
- Diagnostic only. No pipeline number changes, no constant changes, no analysis decision changes.
  `N0_METHOD`, `FC_TO_CELLS_PER_L`, `DEPLETION_FRAC`, the QC thresholds, the fit windows and the
  exclusions are all untouched.
- Additive only when extending `09` and `10`: their current outputs must be byte-identical afterwards.
  Do not modify `12_joint_rK_estimator.R`'s default behaviour.
- Write only to new paths and new output files, plus additive sections in `09`/`10`. Never overwrite
  `results/`, `data/` or `scripts/original_scripts/`.
- Write the new code so it would work on a *Candida*-style temperature series too (taxon x
  temperature x replicate), not only this design — the engine is intended to be multi-dataset. Do not
  implement that ingestion here; just do not hard-code assumptions that would block it.
- Report ACTUAL numbers everywhere. State bootstrap counts, bandwidths and any estimator choice, with
  the reason.
- Autonomous; commit in parts: "D5: residual autocorrelation and corrected per-curve SEs",
  "D5: SE vs replicate reproducibility",
  "D5: r-K ridge under vcov, parametric bootstrap and observed replicates",
  "D5: extend 09/10 with N0 form and joint/misspecified recovery",
  "D5: joint estimator under corrected measurement covariance",
  "D5: report (Quarto -> PDF) + PR".
```
