# =============================================================================
# run_all.R - Master pipeline runner (bacterial O2 respiration, figures only)
# =============================================================================
# Numbered running order; each non-app script sources config.R for shared paths.
#
#   01 longdata            raw Oxygen_Data_Long.csv -> tidy long table
#   02 trimming            spline hybrid start/end trimming + diagnostics + meta
#   03 trim_selector  (App) review/override each curve's fit window + exclusions   [optional]
#   04 experiment_inputs (App) per-taxon Ninoc + cell size                           [optional]
#   05 oxygen_fits         nlsLM fits (honours 03 windows/exclusions), N0/R recon
#   06 main_figures        Fig 2, Supp Fig 3, Supp S1, Fig 6 (r vs R RIS)
#   07 fig6_growth_resp    Figure 6: per-cell respiration vs growth (r vs R)
#   08 cutoff_sensitivity  O2_norm >= 0.5 sensitivity fits + dynamics PDFs
#   09 window_sensitivity  fit-window robustness of K, R, r (reviewer test)
#   10 montecarlo_N0       Monte-Carlo of R vs N0 uncertainty + figures
#   11 simulation_recovery synthetic parameter recovery (Fig S2)          [no inputs]
#   12 temperature_cue     temperature-response growth/resp/CUE (Fig 7)
#
# 03 and 04 are interactive Shiny apps ("Run App"): run them BETWEEN 02 and 05 to
# hand-tune fit windows / cell sizes. They are optional - the pipeline runs fully
# automatically without them. 12 uses its own input data/Oxygen_Data_Filtered_CUE.csv
# and is independent of 01-11.
#
# This pipeline produces FIGURES only; the OD600 / flow-cytometry growth-rate
# comparison from the original OxygenModel.R is intentionally not included.
# =============================================================================

script_dir <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  d <- if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE)) else
         tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (is.null(d) || is.na(d) || !nzchar(d)) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
        nzchar(rstudioapi::getActiveDocumentContext()$path))
      d <- dirname(rstudioapi::getActiveDocumentContext()$path) else d <- getwd()
  }
  d
})

run_script <- function(name) {
  message("\n", strrep("=", 70), "\nRunning: ", name, "\n", strrep("=", 70), "\n")
  source(file.path(script_dir, name), local = TRUE)
}

# Non-interactive scripts only (skip the two Shiny apps 03 / 04).
for (s in c("01_longdata.R", "02_trimming.R", "05_oxygen_fits.R",
            "06_main_figures.R", "07_fig6_growth_resp.R",
            "08_cutoff_sensitivity.R", "09_window_sensitivity.R",
            "10_montecarlo_N0.R", "11_simulation_recovery.R")) {
  run_script(s)
}

message("\n", strrep("=", 70))
message("Done: trimming, fits, figures, cutoff + window sensitivity, MC (Fig 2, 6, S1, S2, Supp 3).")
message("Optional interactive tuning: 03_trim_selector.R and 04_experiment_inputs.R (Run App),")
message("then re-run 05 onward. Run 12_temperature_cue.R separately (temperature data).")
message(strrep("=", 70))
