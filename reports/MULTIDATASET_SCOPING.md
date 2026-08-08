# Multi-dataset scoping — what stands in the way of one engine, two datasets

**Status: scoping only. Nothing here is implemented.** Written during D1
(commit `ed5e1a6` + the D1 changes) against `scripts/config.R` and `01`–`10`.

The goal is one engine that ingests both the bacterial series in this repo
(15 taxa × 5 replicates, one temperature) and a *Candida*-style temperature
series (four *C. auris* clades + *C. parapsilosis*, 12 temperatures).

---

## 1. Where the code assumes THIS dataset

### 1.1 The single-temperature split runs straight down the middle

This is the biggest structural problem, and it is not a config item.

| | `01`–`08` (bacterial) | `10` (temperature) |
|---|---|---|
| Series key | `(Taxon, Replicate)` | `(Taxon, Temperature, Replicate)` |
| Input | `data/Oxygen_Data_Long.csv` | `data/Oxygen_Data_Filtered_CUE.csv` |
| Trimming | `02_trimming.R`, spline peak → hybrid end, writes per-curve metadata | inline `which(sm < -1e-7)[1] + 15`, no metadata, no diagnostics |
| Manual windows / exclusions | honoured (`MANUAL_FIT_WINDOWS`, `PLOT_EXCLUDE_POINTS`) | **ignored entirely** |
| Outlier step | none | one round of ±2 SD trimming then refit |
| QC gate | recorded as `fit_ok`, all rows kept | `if (!fit_ok) next` — failures silently dropped |
| N₀ | `N_inoc` per curve from `data/Ninoc.csv`, δ per curve | a single scalar `INOC_CELLS_PER_uL <- 550`, a single `INOC_DELAY_MIN <- 45` |
| Cell carbon | per taxon, from `data/taxon_cell_sizes.csv` | re-derived inline from local `cell_width`/`cell_length`/`C_density_fg_per_um3` |
| Respiration → carbon | `O2_to_C_mass` from molar masses (0.37537) | `mgL_to_mol_per_mL()` then `× 12e15`, a different route |

`10_temperature_cue.R` is effectively a **second, parallel implementation** of
stages 02+05 that shares only `resp_model()` and the QC thresholds from
`config.R`. A *Candida* temperature dataset would land in `10`'s world and
inherit none of `02`/`03`/`05`'s machinery. Generalisation means collapsing
these two paths, not adding a third.

### 1.2 Hard-coded taxon names and counts

| Location | What | Effect on a second dataset |
|---|---|---|
| `config.R:320-323` `FIG2_SELECTED_COMBOS` | literal `Bacillus/R4`, `Burkholderia/R2`, `Arthrobacter/R2`, `Yersinia/R1` | Fig 2 silently prints "none matched … skipped" |
| `config.R:301` `TAXON_ORDER <- NULL` | no ordering | *Candida* clades would sort alphabetically, not by clade |
| `config.R:275-289` `TAXON_PALETTE` | 16 colours, recycled | fine for 15 taxa; recycles silently past 16 |
| `04_experiment_inputs.R:62` | fallback `c("Bacillus", "Yersinia")` | wrong defaults when no long table exists yet |
| `03_trim_selector.R` header/comments | describes "OTU = Candida isolate code 1..18" | already stale for *this* dataset — the app was ported from the Candida pipeline |
| `06_main_figures.R` | `scale_colour_brewer(palette = "Dark2")` on `Replicate` | breaks past 8 replicates |

None of these are counts of taxa — the pipeline genuinely does not care how
many taxa there are. The binding constraint is *identity*, not *cardinality*.

### 1.3 Replicate naming

`02_trimming.R:277,356` do `as.integer(str_remove(Replicate, "^R"))` to sort
curves. Replicates that are not `R<integer>` sort as `NA` — no error, just
arbitrary ordering, and therefore arbitrary `curve_code` assignment (`C001`…),
which is the key the 03 app and the manual-window CSVs are anchored to. Several
places `toupper()` the replicate (`config.R:292-293`, `05:74`, `03:493`), so
replicate labels are effectively case-insensitive — a dataset with `a`/`A` as
distinct replicates would silently merge.

### 1.4 Plate / well layout

There is **none**. No script reads a plate map, a well ID, or a channel number.
The long table is already `(Taxon, Replicate, Time, Oxygen)`. That is a feature:
plate geometry is somebody else's problem, upstream of `data/`. A *Candida*
dataset must arrive already melted the same way.

### 1.5 Table schemas that would need a per-dataset identity

| File | Schema | Assumption to break |
|---|---|---|
| `data/Ninoc.csv` | `Taxon, Replicate, N_inoculation_cells_per_L, delta_Ninoc_to_N0_min` | no `Temperature` column, so one inoculum per taxon×replicate only |
| `data/taxon_cell_sizes.csv` | `Taxon, cell_width_um, cell_length_um, cell_volume_um3, carbon_density_fg_per_um3, cell_carbon_fg` | rod geometry. *Candida* is an ellipsoid/budding yeast — the rod formula in `config.R:206-211` and `04_experiment_inputs.R::rod_volume()` is wrong for it, and cell carbon is ~2–3 orders of magnitude larger |
| `data/OD_r_FC_r.csv` | `Taxon, Replicate, OD_Initial, OD_Final, FC_Initial, FC_Final, Oxygen_r, Time` | growth rates are formed as `(log(Final) − log(Initial)) / Time` in `06`; a dataset with a full growth curve rather than two endpoints has nowhere to go |
| `data/Cell_Counts.csv` | as above minus `Oxygen_r`, `Time` | **read by no automated stage at all** (see §3) |

### 1.6 Units

All implicit, none declared anywhere in `data/`:

- `Time` — **minutes**. `MIN_TO_H <- 60` and every `r_per_hour` depend on it.
- `Oxygen` — **mg O₂ L⁻¹**. Only the *ratio* enters the fit (`Oxygen/O0`), but
  `C_tot_mg_per_L`, `R`, and `10`'s `mgL_to_mol_per_mL()` all assume mg L⁻¹.
- `Temperature` — **°C**. `10` adds 273.15 in six places.
- `N_inoculation_cells_per_L` — **cells L⁻¹**, but `10` works in cells µL⁻¹ and
  `04`'s help text says counts are cells µL⁻¹ converted by `1e6`. Three unit
  conventions for one quantity (see the D1 report, §E2).
- Cell carbon — **fg C cell⁻¹**.

### 1.7 What in `config.R` would have to become per-dataset

| Currently a global | Should be per-dataset |
|---|---|
| `OXYGEN_LONG_CSV`, `CUE_CSV`, `NINOC_CSV`, `TAXON_CELL_SIZES_CSV` | one input block per dataset |
| `N_inoculation_cells_per_L`, `INOC_DELAY_MIN` | per dataset (and, for a temperature series, potentially per temperature) |
| `CELL_WIDTH_UM`, `CELL_LENGTH_UM`, `CARBON_DENSITY_FG_PER_UM3`, and the rod-volume formula | per dataset — the *shape model* changes, not just its parameters |
| `RESPIRATORY_QUOTIENT` | per dataset |
| `FIT_LOWER` / `FIT_UPPER` (`r` capped at 0.1 min⁻¹) | per dataset — yeast at low temperature will sit near the lower bound |
| `FIG2_SELECTED_COMBOS`, `TAXON_ORDER`, `TAXON_PALETTE` | per dataset |
| `MANUAL_FIT_WINDOWS_CSV`, `PLOT_EXCLUDE_CSV` | per dataset (curation is not transferable) |
| `N0_BACKPROJECT`, the 8 QC thresholds, `USE_DRAWDOWN_WINDOW` | **keep global** — these are method decisions, and making them per-dataset would let two datasets be analysed under different rules |

That last row matters. The line to draw is: *dataset facts* become per-dataset,
*method decisions* stay global. Otherwise "multi-dataset" quietly becomes
"two different methods".

---

## 2. A minimal data contract

Both datasets would have to supply exactly these four things. Everything else is
derived.

### C1 — the O₂ series (required)

One long CSV, one row per reading:

| column | type | units | notes |
|---|---|---|---|
| `dataset_id` | chr | — | e.g. `bacteria_2026`, `candida_auris` |
| `group` | chr | — | the biological grouping variable (today's `Taxon`; a *Candida* clade) |
| `replicate` | chr | — | free text; **not** required to be `R<n>` |
| `temperature_c` | dbl | °C | required; a single-temperature dataset repeats one value |
| `time_min` | dbl | min | 0 = first reading of that series |
| `oxygen` | dbl | declared by C4 | |

Series key is always `(dataset_id, group, temperature_c, replicate)`. The
single-temperature case stops being a special case.

### C2 — inoculum (required)

| column | type | units |
|---|---|---|
| `dataset_id`, `group`, `replicate`, `temperature_c` | as C1 | |
| `n_inoc_cells_per_L` | dbl | cells L⁻¹ |
| `delta_inoc_to_first_reading_min` | dbl | min |

`temperature_c` may be `NA` to mean "applies at every temperature". This
replaces both `data/Ninoc.csv` **and** `10`'s scalar `INOC_CELLS_PER_uL` /
`INOC_DELAY_MIN`.

### C3 — cell carbon (required)

| column | type | units |
|---|---|---|
| `dataset_id`, `group` | as C1 | |
| `cell_carbon_fg` | dbl | fg C cell⁻¹ |
| `shape_model` | chr | `rod`, `sphere`, `ellipsoid`, or `direct` |
| `dim_a_um`, `dim_b_um`, `carbon_density_fg_per_um3` | dbl | optional; only if `shape_model != "direct"` |

Supplying `cell_carbon_fg` directly is the escape hatch that lets a yeast
dataset in without teaching the engine yeast geometry.

### C4 — the dataset manifest (required)

One row per dataset, declaring what is currently implicit:

| column | example |
|---|---|
| `dataset_id` | `candida_auris` |
| `oxygen_units` | `mg_O2_per_L` |
| `time_units` | `min` |
| `respiratory_quotient` | `1` |
| `r_lower_per_min`, `r_upper_per_min` | `1e-5`, `0.05` |
| `group_order` | `Clade_I|Clade_III|Clade_IV|Clade_V|parapsilosis` |
| `manual_windows_csv`, `exclusions_csv` | paths, may be empty |
| `fig2_examples` | `Clade_I/R2|parapsilosis/R1` |

### Optional C5 — an independent growth-rate reference

Only needed to reproduce Figs 3/4/5:

`dataset_id, group, replicate, temperature_c, method, r_per_min`

with `method ∈ {OD600, flow_cytometry, ...}`. Long, not wide — this removes
`06`'s hard-wired `r_OD600`/`r_FC` columns and the two-endpoint
`(log(Final) − log(Initial)) / Time` assumption, which does not survive contact
with a dataset that has real growth curves.

### What each dataset would have to supply

| | bacterial (this repo) | *Candida* temperature series |
|---|---|---|
| C1 | rename `Oxygen_Data_Long.csv` columns; add `temperature_c` constant; add `dataset_id` | melt to long; already has a temperature axis |
| C2 | `Ninoc.csv` renamed; already per-replicate | currently a scalar — must be measured or declared per clade |
| C3 | `taxon_cell_sizes.csv` renamed, `shape_model = rod` | new; `shape_model = ellipsoid` or `direct` |
| C4 | one row | one row; wider `r` bounds |
| C5 | `OD_r_FC_r.csv` melted to long | optional |

---

## 3. Things worth fixing while doing this

1. **`data/Cell_Counts.csv` is consumed by no automated stage.** Only the `04`
   Shiny app touches it, and only as a fallback list of taxon names. It carries
   `FC_Initial` that disagrees with `data/OD_r_FC_r.csv` in 72 of 75 rows.
2. **Three unit conventions for the inoculum** (cells L⁻¹ in `Ninoc.csv`,
   cells µL⁻¹ in `10`, "counts × 1e6" in `04`'s help text) is exactly the sort of
   thing C2 + C4 exist to stop.
3. **`10` drops QC failures with `next`** while `01`–`08` record `fit_ok` and
   keep the row. Under one engine these must agree, or the two datasets are
   filtered differently without anybody noticing.
4. **Curve codes (`C001`…) are positional**, assigned by sort order in `02`.
   Adding a dataset, or a replicate, renumbers every curve after it — and the
   manual-window CSVs are keyed on `(Taxon, Replicate)` but the 03 app displays
   and reasons in curve codes. A stable, content-derived key would be safer.

---

## 4. Rough shape of the work (not a commitment)

1. Land the contract as four CSV schemas + a validator that fails loudly.
2. Collapse `10`'s inline trim/fit into the `02`+`05` path, keyed on
   `(group, temperature_c, replicate)`. This is the only genuinely hard step and
   must be done with the D1 baseline as the regression test: the bacterial
   numbers must not move.
3. Make `config.R` read the manifest instead of holding dataset constants, while
   keeping the method decisions (`N0_BACKPROJECT`, the QC thresholds, the model)
   global.
4. Only then add the *Candida* dataset.

Step 2 is where reproducibility is at risk. Every number in
`reports/D1_baseline/D1_baseline.pdf` is the baseline it has to reproduce.
