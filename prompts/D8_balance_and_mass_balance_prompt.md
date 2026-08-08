# Claude Code prompt — D8: the balance question, answered without the absolute scale — plus an independent mass-balance check (autonomous)

Run from the project root (`.../work/OxygenModel`). Analysis on existing outputs and the deposited
data. It fits no new oxygen curves by default, changes no estimator, and alters no committed number.
It does two things: establishes the quantities the assay was actually built to measure, using only
parameters that carry none of the unresolved scale assumptions; and computes a completely independent
estimate of carbon-use efficiency from mass balance, as a check on the kinetic one.

NOTE TO USER: launch in an auto-approving mode. Cheap — arithmetic and regressions on tables that
already exist. Branch off `gyd/d7-pirt` and open a PR.

CONTEXT — WHY THIS IS THE RIGHT ANALYSIS.

The purpose of this assay is to measure growth and respiration simultaneously so that the BALANCE
between them can be tracked against drivers such as temperature. D5, D6 and D7 all returned negative
results on the estimator itself — its uncertainty is honest, its structure is adequate, and the
maintenance/yield split is not identifiable here. Meanwhile the one open problem, from D2, is the
absolute respiration scale: the N_inoc discrepancy, the FC-to-cells conversion, the carbon quota and
RQ.

The key observation this prompt tests is that THE BALANCE QUESTION DOES NOT NEED ANY OF THAT.

With `G = r * q` and `R = c * K` (q the per-cell carbon quota, c collecting the O2->C conversion and
1/N0):

    CUE = G / (G + R) = 1 / (1 + (c/q) * K/r)

Maximising CUE means MINIMISING `K/r`, and the constant `c/q` drops out of the argmin entirely. So
`T_opt(CUE)` is identifiable WITHOUT the FC constant, the carbon quota, RQ, or the absolute N0 scale.
Only the LEVEL of CUE depends on those. Likewise the difference in thermal sensitivities,
`E_K - E_r`, comes straight from the fitted parameters.

THE ONE CONDITION, which must be tested rather than assumed: this holds only if N0 is
TEMPERATURE-INDEPENDENT. It is not — cells genuinely grow during the interval before the fit window,
and faster when warm. With `delta_t` around 124 min and r spanning roughly 0.01-0.03 min^-1, the true
N0 may differ severalfold between cold and warm vials. The current back-projection corrects for this
using the final-temperature r, which over-corrects, worst at the extremes — precisely the axis of
interest. PART B quantifies that, rather than waving at it.

SECOND STRAND. `data/OD_r_FC_r.csv` carries `FC_Initial` and `FC_Final` per replicate, and the raw
traces give total oxygen drawdown. Biomass produced against carbon respired over the same run is an
INTEGRATED CUE obtained by mass balance — independent of the kinetic model, of N0, and of the
back-calculation. Comparing it against the kinetic CUE is a genuine external check on the method, and
any systematic gap is informative, because the integrated version includes the low-growth tail where
maintenance dominates.

---

```
Work AUTONOMOUSLY end to end; commit in parts; print a summary. Read first: scripts/config.R (the
carbon constants, the N0 block, DEPLETION_FRAC), scripts/05_oxygen_fits.R (how G, R and CUE are
formed), scripts/11_temperature_cue.R (INCLUDING line ~40, which sets its own INOC_DELAY_MIN <- 45,
overriding config — this is a known inconsistency and must be handled explicitly, not silently),
results/tables/{oxygen_fit_results.csv, oxygen_results_with_R.csv, Oxygen_All_Long.csv},
data/{OD_r_FC_r.csv, taxon_cell_sizes.csv, Ninoc.csv}, and the D5/D6/D7 reports. Branch:
`git switch -c gyd/d8-balance gyd/d7-pirt`.

PART A - the scale-free balance quantities
Work from the FITTED parameters r and K only. Do not convert to per-cell carbon anywhere in this part.
- For the *Pseudomonas* temperature series: fit Boltzmann-Arrhenius (and a unimodal alternative where
  the data support it) to r(T) and to K(T) separately. Report `E_r`, `E_K`, and `E_K - E_r` with
  intervals. State plainly which is more temperature-sensitive and by how much.
- Compute `K/r` against temperature and locate its minimum — this is `T_opt(CUE)` and it is
  independent of `c/q`. Report it with an interval (bootstrap over curves), and compare it against
  the `T_opt(CUE)` the current pipeline reports via `11_temperature_cue.R`. They should agree; if
  they do not, that locates a problem.
- Report the fold-change in `K/r` across the measured temperature range — the shift in the balance —
  which is also scale-free.
- Do the same for the 15-taxon set at 30 C as a between-taxon comparison, being explicit that a
  single temperature gives no thermal sensitivity there, only a ranking.
- State, in one paragraph, exactly which of these quantities carry NO dependence on
  FC_TO_CELLS_PER_L, the carbon quota, RQ, or N0. That paragraph is the point of the whole prompt.

PART B - how much does the N0 temperature dependence contaminate them?
- Quantify the true problem: report `delta_t` and `r` per curve, the implied `exp(r*delta_t)` factor,
  and how much it varies across the temperature range. This is the differential in true N0 between
  cold and warm vials that any N0 treatment must handle.
- Recompute `E_K` and `T_opt(CUE)` under three N0 treatments — the current default, a fixed
  temperature-independent N0, and the ramp-aware variant if it can be reconstructed cheaply — and
  report how far each quantity moves. `E_r` and the argmin of `K/r` should be UNAFFECTED by
  construction; verify that and report it as a check, not an assumption.
- Handle the `11_temperature_cue.R` INOC_DELAY_MIN override explicitly: report results both as the
  script currently computes them and with the config value, and state which the published figure used.
- Then answer the design question with numbers: if `delta_t` were reduced to ~10 min by
  pre-equilibrating the medium at temperature before inoculation, what would `exp(r*delta_t)` become,
  and how much of the contamination would remain? Give the figure for the coldest and warmest points.

PART C - integrated CUE from mass balance
- Per replicate: biomass carbon produced = `(FC_Final - FC_Initial) * FC_TO_CELLS_PER_L * cell_carbon`;
  carbon respired = total O2 drawdown over the run * the O2->C conversion * vial volume. Integrated
  CUE = biomass C / (biomass C + respired C).
- Derive it explicitly in the report, stating every assumption and where each number comes from. Be
  careful about the interval: the kinetic CUE is instantaneous over the fit window, the mass-balance
  CUE integrates the whole run including the post-window decline. They are NOT the same estimand and
  the comparison must say so.
- Report both per replicate and per taxon, with the ratio. Test whether the gap correlates with how
  much of the run lies outside the fit window, and with the fraction of oxygen consumed post-window —
  if maintenance dominates the tail, the integrated value should be systematically lower.
- Note which inputs each version needs. The mass-balance version still needs FC_TO_CELLS_PER_L and
  the carbon quota for the ABSOLUTE value, but the RATIO of the two CUEs cancels the quota — report
  that ratio as the scale-free comparison.
- Flag any replicate where the mass balance is physically impossible (CUE outside [0,1], or biomass
  carbon exceeding total carbon). Those are diagnostic of the conversion constants and must be
  reported, not dropped.

PART D - what it means
- One paragraph, written to be quoted: which conclusions about the growth-respiration balance are
  secure without resolving the absolute scale; which are not; and what the mass-balance check says
  about the kinetic CUE.
- A short recommendation on reporting: which quantities the papers should lead with, and which should
  carry a stated uncertainty rather than a point value.
- If PART C shows the two CUE estimates disagree badly, say so plainly — that would be the most
  important finding in the prompt and it must not be softened.

PART E - report
- `reports/D8_balance/D8_balance.qmd` rendered to PDF following `reports/_shared/`, every input
  tracked so it rebuilds from a clean clone (see the runs/ policy in .gitignore). Add its entry to
  `reports/README.md`.
- Figures: `E_r` and `E_K` on one panel; `K/r` against temperature with its minimum marked;
  `T_opt(CUE)` under each N0 treatment; kinetic against mass-balance CUE per replicate.
- Be even-handed. If the scale-free quantities turn out to be materially contaminated by the N0
  temperature dependence, that is the finding and it must be reported as plainly as the alternative.

PART F - PR
- Push and open a PR against `main`. Do not push to `main` or to any upstream branch in the stack.

VERIFY (report all)
1. `E_r`, `E_K`, `E_K - E_r` with intervals; which is more temperature-sensitive.
2. `argmin(K/r)` with an interval, and its agreement with the pipeline's `T_opt(CUE)`.
3. Fold-change in `K/r` across the temperature range.
4. The explicit list of quantities carrying no dependence on FC_TO_CELLS_PER_L, the quota, RQ or N0.
5. `exp(r*delta_t)` per curve and its range across temperature; the movement in `E_K` and
   `T_opt(CUE)` under each N0 treatment; confirmation `E_r` and `argmin(K/r)` are unaffected.
6. The INOC_DELAY_MIN override: results both ways, and which the published figure used.
7. The pre-equilibration projection at delta_t ~10 min, coldest and warmest.
8. Integrated vs kinetic CUE per replicate and per taxon, their ratio, and whether the gap tracks the
   post-window oxygen fraction; any physically impossible replicate named.
9. Report renders from a clean clone; page count; inputs tracked.
10. `results/`, `data/`, `config.R` and `scripts/original_scripts/` byte-identical, by checksum.

CONSTRAINTS
- No estimator change, no constant change, no committed number altered. New outputs to new files.
- PART A must not convert to per-cell carbon anywhere. If a quantity needs the quota or N0, it does
  not belong in PART A — put it in PART B or C and say so.
- Do not silently reconcile the `11_temperature_cue.R` delay override. Report both.
- The two CUE estimates are different estimands. Do not present their difference as an error in
  either without establishing which interval each refers to.
- Report ACTUAL numbers everywhere, with intervals where a bootstrap is cheap.
- Autonomous; commit in parts: "D8: scale-free balance quantities from r and K",
  "D8: N0 temperature dependence and the pre-equilibration projection",
  "D8: integrated CUE from mass balance", "D8: report (Quarto -> PDF) + PR".
```
