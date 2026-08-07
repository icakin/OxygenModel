# Claude Code prompt — D6: the two misspecifications that actually bite — lag and saturation (autonomous)

Run from the project root (`.../work/OxygenModel`). Diagnostic first, then a model extension. It does
not change the default estimator or any committed number; it establishes whether the model's two
measured vulnerabilities are present in the real curves, and if so what a corrected model does to the
estimates.

NOTE TO USER: launch in an auto-approving mode. Work on a branch off `gyd/d5-uncertainty` and open a
PR. Moderate cost — refits over a grid, plus model comparison. No large MCMC unless PART D escalates.

CONTEXT — THIS PROMPT REPLACES THE ORIGINAL D6, WHICH WAS WRONG.

D6 was originally specified as a one-stage hierarchical fit whose purpose was to propagate the
within-curve r-K covariance, on the premise that the two-stage pipeline understated its uncertainty.
D5 tested that premise and refuted it:

- Residuals are NOT autocorrelated (median AR(1) rho = -0.034, DW 2.05, 0/75 curves failing
  Ljung-Box). `05` fits the RAW `Oxygen` column; the spline in `02_trimming.R` only selects the
  window. The AR(1)-corrected SE is 0.97x the naive one.
- The within-curve SE is 37x smaller than the between-replicate spread, but that gap is genuine
  replicate variation, not understated measurement error - the correction closes -0.2% of it, and
  inflating S_i a hundredfold moves the posterior by ~10%.
- The r-K ridge is real (corr = -0.985) but PROTECTIVE: the `(e^{rT}-1)/r` factors cancel, giving
  `log R = log K + r*dt`, and the anti-correlation makes that combination LESS variable than if r and
  K were independent (sd 0.0035 vs 0.0092).
- Contrast resolution is unchanged under every covariance treatment tried: 104/105 either way.

So uncertainty propagation is not the problem. The measured problem is BIAS from misspecification:

    none              +0.09%      (in per-cell R)
    optode drift      -1.1%
    initial lag      +18.6%
    saturating growth -78.1%

Drift is a non-issue. Lag matters. Saturation is catastrophic and had never been tested. D6 asks
whether either is present in the real data, and what to do about it.

A HYPOTHESIS TO TEST, NOT ASSUME. `08_window_sensitivity.R` found that a +/-6 min shift of the window
edges moves K and per-cell R by 15-34% while moving r by ~0.3%. If growth or oxygen uptake saturates
near the end of some windows, then moving the end changes how much saturation is included - which
would produce exactly that asymmetry between r and K. If that link holds, the window sensitivity and
the saturation bias are one cause, and the fix is a model term rather than tighter curation.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first: scripts/config.R (the
model, bounds, DEPLETION_FRAC, the N0 block), scripts/05_oxygen_fits.R, scripts/02_trimming.R (what
O2_fit is and how the window is chosen), scripts/08_window_sensitivity.R,
scripts/10_simulation_recovery.R (including the misspecification arms D5 added - reuse their
generators, do not rewrite them), and reports/D5_estimator_uncertainty/. Branch:
`git switch -c gyd/d6-model-structure gyd/d5-uncertainty`.

PART A - is saturation present in the real curves?
- For every curve, fit the current model over its committed window and examine the residuals as a
  function of (i) time within the window and (ii) ABSOLUTE oxygen concentration. Systematic curvature
  in the residuals late in the window, or below some O2 level, is the signature.
- Report per curve: the residual trend over the last third of the window, and the minimum O2 reached.
- Distinguish the two candidate mechanisms rather than assuming one:
  * OXYGEN LIMITATION - respiration becomes O2-dependent as O2 falls, i.e.
    dO2/dt = -q N(t) O2/(Km + O2). Predicts the deviation scales with absolute O2 and appears in every
    curve that gets low enough, regardless of taxon.
  * GROWTH SATURATION - the culture approaches carrying capacity, so N(t) stops being exponential.
    Predicts the deviation scales with elapsed growth (r*t or N/N_max), not with O2 level.
  These make different predictions; report which the residuals support, or that they cannot be
  separated in this data.
- Report how many of the 75 curves show a detectable deviation, and how large.

PART B - test the window-sensitivity hypothesis
- Take the per-curve window sensitivity already computed in `08` (dK_pct, dR_pct) and correlate it
  against the PART A saturation measures - residual curvature, minimum O2 reached, elapsed growth at
  the window end.
- If the correlation is strong, state plainly that the window sensitivity is a symptom of
  misspecification rather than a curation problem, and quantify how much of it saturation explains.
- If it is weak, say so - that would mean the two findings are independent and both need addressing.
- This is the pivot of the whole prompt. Report it before going further.

PART C - extend the model
Implement as ALTERNATIVE fits alongside the current one. Do NOT change the default estimator.
- M0: current model, as reference.
- M1: + explicit lag/equilibration. The onset interval is instrumental, not biological (the paper's
  own methods describe a ~45 min stabilisation during which the reported values drift upward as the
  optode settles), so a lag term here is a nuisance parameter absorbing that transient, not a
  biological lag. Say which you have implemented and why.
- M2: + the mechanism PART A supports (O2 limitation via a Km term, or saturating growth). If PART A
  cannot separate them, implement both and compare.
- M3: M1 + M2, if both are supported.
- Fit with `nlsLM` first. Only if identification fails - non-convergence, parameters at bounds, or a
  ridge that makes the extra term unidentifiable - escalate to a Bayesian fit with weak priors, and
  say explicitly that the escalation is for IDENTIFICATION, not for uncertainty. D5 established that
  uncertainty propagation is not the issue here; do not reintroduce it as a justification.

PART D - compare
- Per curve and per taxon: AIC/BIC, residual structure, and the fitted r, K and per-cell R under each
  model. Report the systematic shift in per-cell R from M0 to each alternative.
- THE KEY TEST: re-run the `08` window sweep under the best-supported model. If the extension is
  right, window sensitivity should FALL substantially. Report dK_pct and dR_pct under M0 and under
  the alternative, side by side. This is the strongest available evidence that the extension is
  capturing something real rather than just adding parameters.
- Report what happens to the between-taxon ordering of per-cell R and to the Fig 6 growth-respiration
  slope under each model. State whether any published conclusion moves.
- Guard against overfitting: an extra parameter will always improve fit. Use the window-sensitivity
  reduction and the residual structure as the real tests, not the likelihood alone.

PART E - what it means
- One paragraph, written to be quoted: is the current model adequate for the follow-on work; if not,
  which extension is needed; and what it would change about the published estimates.
- Be explicit about the scope of any recommendation. Changing the default estimator is a decision for
  the repository owner, not this prompt.

PART F - report
- `reports/D6_model_structure/D6_model_structure.qmd` rendered to PDF, following
  `reports/_shared/`: self-contained, provenance footer, figures in its own `figures/`, and every
  input it reads tracked so it rebuilds from a clean clone (see the runs/ policy in .gitignore).
- Include: the residual diagnostics, the mechanism discrimination, the window-sensitivity
  correlation, the model comparison table, the window sweep under each model, and the effect on
  ordering and on Fig 6.
- Be even-handed. If saturation turns out to be rare or small in the real curves, say so as plainly
  as the alternative - that would mean the -78.1% simulation arm was pessimistic and the current
  model is fine, which is a good result and closes the question.
- Add its entry to `reports/README.md`.

PART G - PR
- Push and open a PR against `main`. Body: what was tested, whether the vulnerabilities are present,
  and an explicit statement of what does and does not change in the default pipeline. Do not push to
  `main`, `gyd/packaging`, `gyd/docs` or `gyd/d5-uncertainty`.

VERIFY (report all)
1. Residual structure per curve against time-in-window and against absolute O2; how many curves show
   a detectable deviation and how large.
2. Which mechanism the residuals support - O2 limitation, growth saturation, or indistinguishable -
   with the evidence for the discrimination.
3. The correlation between per-curve window sensitivity and the saturation measures; how much of the
   15-34% K/R movement saturation explains.
4. Model comparison M0-M3: AIC/BIC, convergence, any parameters at bounds; whether Bayesian
   escalation was needed and for which model.
5. Window sweep under M0 vs the best alternative: dK_pct and dR_pct side by side.
6. Systematic shift in per-cell R; between-taxon ordering; Fig 6 slope under each model.
7. Report renders from a clean clone; page count; .qmd and its inputs committed.
8. `results/`, `data/`, `config.R` and `scripts/original_scripts/` byte-identical, verified by checksum.

CONSTRAINTS
- The default estimator does not change. Alternatives are additional fits writing to new files.
- No constant changes: N0_METHOD, FC_TO_CELLS_PER_L, DEPLETION_FRAC, the QC thresholds, the fit
  windows and the exclusions are untouched.
- Reuse D5's misspecification generators in `10_simulation_recovery.R` rather than writing new ones,
  so the simulated and observed diagnostics are on the same footing.
- Do NOT justify a Bayesian fit by uncertainty propagation. D5 refuted that. If you escalate, it is
  for identification, and say so.
- Test the hypothesis in PART B honestly. A weak correlation is a real answer and must be reported as
  one; do not go looking for a stronger framing.
- Report ACTUAL numbers everywhere.
- Autonomous; commit in parts: "D6: residual diagnostics - is saturation present",
  "D6: discriminate O2 limitation from growth saturation",
  "D6: window sensitivity against saturation measures",
  "D6: lag and saturation model variants", "D6: window sweep under each model + comparison",
  "D6: report (Quarto -> PDF) + PR".
```
