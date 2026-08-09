# =============================================================================
# regenerate_ninoc.R - Regenerate data/Ninoc.csv in code from the depletion N0.
# =============================================================================
# Ninoc.csv (the input for the "initial" N0 route) was derived OFFLINE with the
# old placeholder FC_TO_CELLS_PER_L and was never updated after the constant was
# corrected to the dilution-chain value (10,255,100). So the initial and
# depletion routes currently differ by the constant ratio (~11.2x) from the file
# alone, and the offline derivation additionally clamped some per-curve implied
# constants at 1e6.
#
# This regenerates Ninoc.csv IN CODE from the depletion N0 that 05 now computes
# with the corrected constant, so the forward route reproduces the depletion N0
# by construction (no typed constant, no clamp):
#
#     N_inoculation = N0_depletion * exp(-r * delta_Ninoc_to_N0_min)
#
# The per-curve delta (inoculation -> fit start) is a geometric quantity that
# does not depend on the constant, so it is preserved from the existing file;
# only the density is recomputed.
#
# RUN ORDER: after 05 has run with N0_METHOD = "depletion" (so
#            oxygen_results_with_R.csv reflects the corrected constant), and
#            after the packaging PR is merged (so the "nothing moved" baseline
#            stays clean, per Gab). Then re-run the pipeline to refresh outputs.
#
# INPUT  results/tables/oxygen_results_with_R.csv  (N0_cells_per_L, r_per_minute)
#        data/Ninoc.csv                            (delta_Ninoc_to_N0_min)
# OUTPUT data/Ninoc.csv                            (regenerated in place)
#        data/Ninoc_preRegen_backup.csv            (copy of the stale file)
# RUN    Rscript scripts/regenerate_ninoc.R
# =============================================================================

.this_dir <- local({
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
source(file.path(.this_dir, "config.R"))
suppressPackageStartupMessages({ library(dplyr) })

stopifnot(file.exists(RESULTS_FINAL_CSV), file.exists(NINOC_CSV))
if (!identical(N0_METHOD, "depletion"))
  warning("N0_METHOD is not 'depletion'; oxygen_results_with_R.csv N0 may not reflect the ",
          "corrected constant. Run 05 with N0_METHOD = 'depletion' before regenerating.")

res <- readr::read_csv(RESULTS_FINAL_CSV, show_col_types = FALSE) %>%
  dplyr::transmute(Taxon = as.character(Taxon), Replicate = as.character(Replicate),
                   r_per_minute = as.numeric(r_per_minute), N0_cells_per_L = as.numeric(N0_cells_per_L))
old <- readr::read_csv(NINOC_CSV, show_col_types = FALSE) %>%
  dplyr::transmute(Taxon = as.character(Taxon), Replicate = as.character(Replicate),
                   old_N = as.numeric(N_inoculation_cells_per_L),
                   delta_Ninoc_to_N0_min = as.numeric(delta_Ninoc_to_N0_min))

new <- old %>%
  dplyr::inner_join(res, by = c("Taxon", "Replicate")) %>%
  dplyr::mutate(N_inoculation_cells_per_L = N0_cells_per_L * exp(-r_per_minute * delta_Ninoc_to_N0_min)) %>%
  dplyr::filter(is.finite(N_inoculation_cells_per_L), N_inoculation_cells_per_L > 0)

# report the shift vs the stale file (should be ~ the constant ratio, ~11.2x)
cat(sprintf("\nRegenerated %d rows.  new/old N_inoculation ratio: median %.2fx (range %.2f-%.2f)\n",
            nrow(new), median(new$N_inoculation_cells_per_L / new$old_N),
            min(new$N_inoculation_cells_per_L / new$old_N),
            max(new$N_inoculation_cells_per_L / new$old_N)))

out <- new %>% dplyr::arrange(Taxon, Replicate) %>%
  dplyr::select(Taxon, Replicate, N_inoculation_cells_per_L, delta_Ninoc_to_N0_min)

file.copy(NINOC_CSV, file.path(data_dir, "Ninoc_preRegen_backup.csv"), overwrite = TRUE)
readr::write_csv(out, NINOC_CSV)
cat("Backed up stale file to data/Ninoc_preRegen_backup.csv; wrote regenerated data/Ninoc.csv.\n")
cat("Forward route (N_inoc * exp(r * delta)) now reproduces the depletion N0 exactly; no 1e6 clamp.\n")
