# data/ — data dictionary

Input files for the oxygen pipeline. For every column: what it measures, when in
the protocol it was taken, its units, and which script consumes it.

**Timebase.** t = 0 of every oxygen series is **inoculation/sealing on the
reader**: the recording includes the ~45-min stabilization period, visible in
the deposited traces as an initial upward drift that peaks at ~43-52 min (the
warm-up of the reader optics). `02_trimming.R` starts each fit window at that
peak, so `delta_Ninoc_to_N0_min` (= the peak time) spans inoculation to fit
start, consistent with `INOC_DELAY_MIN = 0` in `config.R`. All times below are
minutes on this clock.

---

## Oxygen_Data_Long.csv
Raw dissolved-oxygen time series, one row per reading. Recorded continuously at
~1-min intervals over ~200 min, starting after the ~45-min stabilization.

| column | meaning | units | consumed by |
|---|---|---|---|
| Time | minutes since recording start | min | 01 |
| Taxon | bacterial taxon (15 levels) | - | 01 |
| Replicate | vial (R1-R5) | - | 01 |
| Oxygen | dissolved O2 concentration | mg O2 L^-1 | 01 |

`01_longdata.R` tidies this into `results/tables/Oxygen_All_Long.csv`; the full
trace is also read by `config.R::load_depletion_table()` to find each vial's O2
depletion time (the N0 anchor of the depletion route in 05).

## OD_r_FC_r.csv
Per-replicate biomass reference measurements. **The primary file** for the
OD/flow-cytometry growth validation (Figs 3-5) and for the depletion-route N0.

| column | meaning | when taken | units | consumed by |
|---|---|---|---|---|
| Taxon, Replicate | identifiers | - | - | 05, 06 |
| OD_Initial | optical density (OD600) | ~45 min post-inoculation, **parallel vials** | dimensionless | 06 (r_OD600) |
| OD_Final | optical density (OD600) | end of run, **SensorDish vials** | dimensionless | 06 (r_OD600) |
| FC_Initial | flow-cytometry count of the measured (diluted, SYBR-stained) sample | ~45 min post-inoculation, **parallel vials** | events uL^-1 | 06 (r_FC) |
| FC_Final | flow-cytometry count, same protocol | end of run, **SensorDish vials** | events uL^-1 | 05 (depletion N0), 06 (r_FC) |
| Oxygen_r | oxygen-model growth rate for this curve (legacy copy; 06 recomputes it from results/) | - | min^-1 | none (legacy) |
| Time | elapsed time between the initial and final biomass measurements (renamed `Duration` inside 06) | - | min | 06 |

**FC_Initial is NOT an at-inoculation count.** No flow-cytometry measurement was
taken from the SensorDish vials at inoculation. FC_Initial/OD_Initial come from
the parallel vials at ~45 min (matching the stabilization period); FC_Final and
OD_Final come from the SensorDish vials themselves at the end of the run.
Conversion of FC events uL^-1 to cells L^-1 of culture uses `FC_TO_CELLS_PER_L`
in `config.R` (derived there from the dilution chain).

## Ninoc.csv
Per-replicate starting cell densities for the "initial" N0 route.
**Derived, not measured** - regenerated in code by `scripts/regenerate_ninoc.R`
from the depletion-route N0: `N_inoculation = N0_depletion * exp(-r * delta)`,
so the forward route `N0 = N_inoc * exp(r * delta)` reproduces the depletion N0
exactly.

| column | meaning | units | consumed by |
|---|---|---|---|
| Taxon, Replicate | identifiers | - | 05 |
| N_inoculation_cells_per_L | back-projected cell density at t = 0 of the recording clock | cells L^-1 | 05 (only when `N0_METHOD = "initial"`) |
| delta_Ninoc_to_N0_min | delay from t = 0 to the fit-window start for that curve (from 02's trimming metadata) | min | 05 |

Timepoint note: this column is anchored at t = 0 = inoculation. The measured
FC_Initial/OD_Initial reference sits at ~45 min, which is close to each curve's
fit-window start (delta ~ 43-52 min), so the natural matched comparison is
FC_Initial against N0 at fit start (equivalently N_inoc * exp(r * delta)), not
against N_inoculation directly.

## Ninoc_preRegen_backup.csv
Frozen copy of the previous Ninoc.csv, which was derived offline with the old
placeholder conversion constant (909,916) and carried a 1e6 clamp. Kept as a
provenance record only; consumed by nothing.

## Cell_Counts.csv
Legacy summary of OD_r_FC_r.csv: each taxon's FC_Initial/OD values replaced by
the rounded per-taxon arithmetic mean (identical value repeated across that
taxon's replicate rows). One measurement summarised - not a second sampling.
Consumed by nothing in the pipeline; retained for provenance.

## taxon_cell_sizes.csv
Per-taxon cell dimensions and carbon content (literature-based).

| column | meaning | units | consumed by |
|---|---|---|---|
| Taxon | identifier | - | config |
| cell_width_um, cell_length_um | cell dimensions | um | config |
| cell_volume_um3 | cell volume from dimensions | um^3 | config |
| carbon_density_fg_per_um3 | assumed carbon density (100 fg C um^-3; see note) | fg um^-3 | config |
| cell_carbon_fg | per-cell carbon = volume x density | fg | 05, 06 (via `cell_carbon_of()`) |

Note on the 100 fg C um^-3 conversion: this is an ASSUMED conversion factor, not a
measured one, and it sets the absolute scale of growth carbon flux (and hence the
CUE level, though not its thermal shape or optimum). It should not be described as
an upper bound for cultured, exponentially growing cells: live-cell dry-mass
density measurements of order 300 fg um^-3, combined with measured bacterial
carbon fractions of 45-50% of dry mass, put a rich-medium heterotroph nearer
135-150 fg C um^-3. Treat 100 as a conservative value within the published range
and report 100 / 150 / 200 as a sensitivity where the absolute CUE level matters.

## Oxygen_Data_Filtered_CUE.csv
Temperature-gradient oxygen series (Pseudomonas sp., 20-40 C).

| column | meaning | units | consumed by |
|---|---|---|---|
| Taxon | "Pseudomonas" | - | 11 |
| Temperature | incubation setpoint | C | 11 |
| Replicate | vial | - | 11 |
| Time | minutes since recording start | min | 11 |
| Oxygen | dissolved O2 | mg O2 L^-1 | 11 |

Note: `11_temperature_cue.R` currently sets its own `INOC_DELAY_MIN <- 45`
(overriding config's 0) for this experiment's N0; see the script header for the
reconciliation status of that override.
