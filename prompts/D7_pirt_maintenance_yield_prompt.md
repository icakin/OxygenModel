# Claude Code prompt — D7: can maintenance and yield be separated from the data already in hand? (autonomous)

Run from the project root (`.../work/OxygenModel`). Analysis on existing outputs. It fits no new
oxygen curves, changes no estimator, and alters no committed number. It asks one question: does the
existing data support a Pirt decomposition of respiration into a maintenance term and a
growth-associated term — and if so, is that a new result or a reparameterisation of Figure 6?

NOTE TO USER: launch in an auto-approving mode. Cheap — regressions on tables that already exist,
plus a simulation to check identifiability. Branch off `gyd/d6-model-structure` and open a PR.

CONTEXT.

The paper's introduction frames respiration as having two components: *"a maintenance component that
is independent of the instantaneous growth rate but varies with factors such as temperature, and a
growth-associated component that is proportional to biomass synthesis."* The implemented model has a
single lumped per-cell respiration R. Separating the two is the Pirt relation:

    R = m + mu / Y

where `m` is maintenance respiration (the intercept, respiration at zero growth) and `Y` the yield
(carbon into biomass per carbon respired). In this repository's units that is per-cell respiration
`R_C_fg_cell_h` against per-cell growth `G_C_fg_cell_h`, on LINEAR axes.

Figure 6 already regresses these two quantities — but as `log_G ~ log_Rc`, log-log, with growth as
the response. That is close to a Pirt analysis and is not one: a log-log slope of 0.283 does not
yield `m` and `Y`. The question this prompt asks is whether simply reparameterising onto linear axes,
with respiration as the response, recovers the two parameters from data already collected.

WHY THIS IS THE RIGHT NEXT STEP. D5 and D6 both returned negative results: the estimator's
uncertainty is honest (D5) and its structure is adequate (D6). The Bayesian rebuild originally
proposed as a FIX is therefore not justified, and should not be revived here. What remains is whether
the method can deliver a NEW capability — the maintenance/yield split — and the cheapest possible
test of that is a reparameterisation, not a new estimator.

TWO IDENTIFIABILITY WORRIES TO TAKE SERIOUSLY, NOT WAVE AT.
- Pirt needs variation in `mu` at otherwise fixed conditions. The 15-taxon experiment varies taxon at
  a single temperature (30 C); the *Pseudomonas* experiment varies temperature within one taxon.
  These give `m` and `Y` quite different meanings, and pooling them would be wrong.
- `G = r * q` where `q` is a per-taxon constant, and `R = K/N0` where `N0` contains `exp(r*delta)`.
  So both axes carry `r`. Whether the regression is measuring physiology or that shared dependence
  has to be established, not assumed. D6 found an analogous case where an apparent predictor turned
  out to be `r` in disguise.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first: the paper's Materials and
methods on respiration and CUE (Eq. 12-14) if a copy is in the repo; scripts/config.R (the carbon
constants, N0 block); scripts/05_oxygen_fits.R (how G and R are formed);
scripts/06_main_figures.R (the Fig 6 RIS model, from ~line 441); scripts/11_temperature_cue.R;
results/tables/{oxygen_results_with_R.csv, data_Fig6_RIS_COLLAPSED_centeredX.csv,
mixed_model_Fig6_RIS_centeredX.txt}; and reports/D5_estimator_uncertainty/ and
reports/D6_model_structure/. Branch: `git switch -c gyd/d7-pirt gyd/d6-model-structure`.

PART A - is the shared r dependence fatal before we start?
- Write out algebraically what `G` and `R` each contain, under the current default
  (`N0_METHOD = "depletion"`). Establish whether a regression of R on G is measuring a physiological
  relationship or a construction.
- Test it: recompute both axes with `N0` held FIXED per taxon (no r dependence in N0 at all) and
  refit. If the Pirt parameters survive that, the relationship is not an artefact of shared `r`. If
  they collapse, report that and stop — there is no decomposition to be had here and the rest of the
  prompt is moot.
- Report this before proceeding. It is the gate.

PART B - the two datasets, kept separate
Fit `R = m + G/Y` (linear, respiration as response) in both settings, and DO NOT pool them:
- (B1) 15 taxa at 30 C, n = 75. Here `m` and `Y` are BETWEEN-TAXON quantities: the intercept is the
  respiration of a hypothetical non-growing taxon, which is a strong extrapolation. Say so.
- (B2) *Pseudomonas* across 20-40 C, from `11_temperature_cue.R`. Here `m` and `Y` are
  WITHIN-taxon and vary with temperature, which is what the introduction actually describes. This is
  the scientifically meaningful one.
- For each: intercept, slope, their intervals, R2, and the implied `m` and `Y` in interpretable units.
  Report the leverage of the extrapolation to G = 0 — how far outside the observed range it sits.
- For (B2), also fit `m` and `Y` as functions of temperature if the data support it, and report
  whether maintenance rises with temperature as the introduction asserts.

PART C - reconcile with Figure 6
- Show explicitly how the published `log_G ~ log_Rc` slope of 0.283 relates to the Pirt
  parameterisation. A log-log slope of ~1 would mean growth-proportional respiration with negligible
  maintenance; a slope well below 1 implies a substantial maintenance term. State what 0.283 implies,
  and whether that is consistent with the `m` estimated directly in PART B.
- If they disagree, that is a finding: the two parameterisations of the same data should not tell
  different stories, and the discrepancy would locate the problem.
- Note that Fig 6's model is `(1 | Taxon)` after D5/D6 found the random-slope form singular; keep the
  same random-effects structure when comparing so the difference is the parameterisation and not the
  hierarchy.

PART D - identifiability
- Simulate from a known Pirt relationship at the observed design (n, spread of G, observed residual
  scale) and check that `m` and `Y` are recovered. Report bias and the sampling correlation between
  the intercept and slope estimates — if they are strongly anti-correlated the split is weakly
  identified regardless of fit quality, which is the standard failure mode for this analysis.
- Report the minimum spread in `G` that would be needed to separate `m` from `Y` at useful precision,
  and whether each dataset meets it.
- Be explicit if the honest answer is "the data constrain the product but not the split".

PART E - what it means
- One paragraph, written to be quoted: can maintenance and yield be separated from data already in
  hand; if so from which experiment and with what precision; and if not, what measurement would be
  needed. Concretely for the last part — a wider range of growth rates at fixed temperature (e.g. a
  substrate-limited series or a dilution series) is the classical route; say whether that is what is
  required.
- State plainly whether this is a NEW result or a reparameterisation of an existing figure. Both are
  legitimate outcomes; conflating them is not.

PART F - report
- `reports/D7_pirt/D7_pirt.qmd` rendered to PDF following `reports/_shared/`, with every input it
  reads tracked so it rebuilds from a clean clone. Add its entry to `reports/README.md`.
- Include: the algebra from PART A and the fixed-N0 test; both regressions with their extrapolation
  leverage shown on the figures; the Fig 6 reconciliation; the identifiability simulation.
- Be even-handed. "The split is not identifiable from this design" is a perfectly good result and
  closes the question cleanly — report it as plainly as the alternative. Do not stretch a weak
  intercept into a maintenance estimate.

PART G - PR
- Push and open a PR against `main`. Do not push to `main` or to any upstream branch in the stack.

VERIFY (report all)
1. The algebra of what G and R contain; the fixed-N0 test result; whether the gate passed.
2. B1 and B2 separately: intercept, slope, intervals, R2, implied m and Y, extrapolation leverage.
3. Whether maintenance rises with temperature in B2, with the numbers.
4. The relationship between the 0.283 log-log slope and the Pirt parameters; whether they agree.
5. Identifiability: recovery bias, intercept-slope sampling correlation, the spread in G needed, and
   whether each dataset meets it.
6. The one-paragraph verdict, and whether this is new or a reparameterisation.
7. Report renders from a clean clone; page count; inputs tracked.
8. `results/`, `data/`, `config.R` and `scripts/original_scripts/` byte-identical, verified by checksum.

CONSTRAINTS
- No new curve fits, no estimator change, no constant change, no committed number altered.
- Do NOT revive the Bayesian rebuild as a fix. D5 and D6 refuted the premises for it. If a Bayesian
  fit helps here it is for IDENTIFICATION of a weakly-determined intercept, and must be justified on
  that basis alone and reported as such.
- Do not pool the 15-taxon and *Pseudomonas* datasets. They estimate different things.
- The gate in PART A is real. If the relationship collapses when N0 is held fixed, stop and report;
  do not proceed to fit parameters to an artefact.
- Report ACTUAL numbers everywhere, including the extrapolation distance to G = 0.
- Autonomous; commit in parts: "D7: algebra and the fixed-N0 gate",
  "D7: Pirt fits, 15-taxon and Pseudomonas, kept separate",
  "D7: reconcile with the Fig 6 log-log parameterisation",
  "D7: identifiability simulation", "D7: report (Quarto -> PDF) + PR".
```
