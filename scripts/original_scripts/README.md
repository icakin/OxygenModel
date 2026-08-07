# original_scripts/ — provenance record (NOT pipeline code)

These are the original monolithic scripts that produced the numbers in the
published paper (Cakin *et al.* 2026, *ISME Communications* 6(1):ycag024).
They are kept here as a historical / provenance record only. They are **not**
part of the current analysis pipeline and are **not** run by `run_all.R`.

The reproducible pipeline is the numbered scripts in the parent `scripts/`
directory (`01_longdata.R` … `13_depletion_frac_sensitivity.R`). Use those to
regenerate every figure and table. The files here are frozen and should not be
edited.

## Contents

- `OxygenModel.R` — original end-to-end oxygen-dynamics model + curve fitting
- `CUE.R` — original carbon-use-efficiency computation
- `MC.R` — original Monte-Carlo N0 sensitivity analysis
- `Simulation.R` — original synthetic parameter-recovery simulation
- `05CUTOFF.R` — original O2-cutoff sensitivity check

## Note on scale

The current pipeline's **absolute** respiration and CUE scale differs from these
originals: the flow-cytometry-to-cells constant (`FC_TO_CELLS_PER_L` in
`config.R`) is now **derived from the dilution chain** rather than tuned to a
target. Slopes, thermal shapes and among-taxon orderings are unchanged. See
`config.R` and the review correspondence for details.
