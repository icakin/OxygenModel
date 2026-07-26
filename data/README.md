# data/

Drop the raw / prepared input CSVs here. Expected files:

- `Oxygen_Data_Long.csv`          columns: Taxon, Replicate, Time, Oxygen
- `Oxygen_Data_Filtered_CUE.csv`  columns: Taxon, Temperature, Replicate, Time, Oxygen
                                  (temperature experiment, used by 10 only)

N0 is reconstructed from config's N_inoculation_cells_per_L plus the per-curve
inoculation->O2-start delay recorded by 02_trimming, so no separate Ninoc CSV
is required. Add these files before running the pipeline.
