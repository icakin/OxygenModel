###############################################################################
# MONTE CARLO SENSITIVITY (NEW MODEL):
#   O2_norm(t) = O2_0 + (resp_tot / r) * (1 - exp(r * t))
#
# Goal:
#   For each Taxon × Replicate, quantify how uncertainty in N0
#   (initial cell density at O2 start) propagates into the inferred
#   per-cell respiration rate:
#
#       R = resp_tot * O2_ref / N0
#
#   where:
#     - resp_tot is the dimensionless total respiration scaling parameter
#       from the new oxygen model (eq. 5 in the Methods).
#     - O2_ref is the initial dissolved oxygen concentration (mg O2 L^-1)
#       used for normalising the oxygen series (eq. 4).
#     - N0 is the initial cell density at the start of the oxygen time series
#       (cells L^-1), estimated from inoculation + r (eq. 7).
#
# Inputs required (from main pipeline):
#   - Tables/oxygen_model_results.csv
#       * columns: Taxon, Replicate, resp_tot (or resp_rate -> renamed),
#                  r_per_minute, fit_ok, etc.
#
#   - Tables/O2start_density_estimates_all.csv
#       * columns: Taxon, Replicate, N0_O2start_cells_per_L, ...
#
#   - Tables/Oxygen_Data_Filtered.csv
#       * columns: Taxon, Replicate, Time, Oxygen
#       * already trimmed to peak–second-inflection, but NOT derivative-trimmed.
#
# Output:
#   - Tables/N0_MC_newmodel_raw_R_draws.csv
#   - Tables/N0_MC_newmodel_summary_per_series.csv
#   - Tables/N0_MC_newmodel_taxon_summary.csv
#   - Diagnostic plots of relative SD of R due to N0 uncertainty
###############################################################################

library(tidyverse)
library(zoo)

# ───────────────────────── 0. Project paths ──────────────────────────────── #

data_dir   <- "data"
plots_dir  <- "plots"
tables_dir <- "Tables"

dir.create(data_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# ───────────────────────── 1. Load model results (new model) ─────────────── #

# Should contain: Taxon, Replicate, resp_tot, r_per_minute, fit_ok, ...
model_results <- read_csv(
  file.path(tables_dir, "oxygen_model_results.csv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate)
  )

# Safety: if the column is still called resp_rate, rename to resp_tot
if ("resp_rate" %in% names(model_results) && !"resp_tot" %in% names(model_results)) {
  model_results <- model_results %>%
    dplyr::rename(resp_tot = resp_rate)
}

# Optional: keep only good fits if you want
model_results <- model_results %>%
  dplyr::filter(is.finite(resp_tot))

print(head(model_results))

# ───────────────── 2. Load N0 estimates at O2 start (cells/L) ────────────── #

N0_raw <- read_csv(
  file.path(tables_dir, "O2start_density_estimates_all.csv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate),
    N0_cellsL = as.numeric(N0_O2start_cells_per_L)
  )

print(head(N0_raw))

# Taxon-level mean & SD of N0 (used for Monte Carlo draws)
N0_stats_taxon <- N0_raw %>%
  dplyr::filter(is.finite(N0_cellsL), N0_cellsL > 0) %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(
    N0_mean = mean(N0_cellsL, na.rm = TRUE),
    N0_sd   = sd(N0_cellsL,   na.rm = TRUE),
    n_N0    = dplyr::n(),
    .groups = "drop"
  )

write_csv(
  N0_stats_taxon,
  file.path(tables_dir, "N0_estimated_stats_by_taxon_MC_newmodel.csv")
)
print(N0_stats_taxon)

# ───────────────── 3. Load trimmed oxygen data & compute O2_ref ─────────── #

# This file is from the main trimming script; columns: Taxon, Replicate, Time, Oxygen
oxygen_filtered <- read_csv(
  file.path(tables_dir, "Oxygen_Data_Filtered.csv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate)
  )

print(head(oxygen_filtered))

# We now reconstruct O2_ref for each Taxon × Replicate in a way
# that mirrors the main fitting pipeline:
#
#   - derivative-based trimming to detect onset of decline
#   - then O2_ref = mean of first 3 Oxygen values in the trimmed df

O2_ref_series <- oxygen_filtered %>%
  dplyr::group_by(Taxon, Replicate) %>%
  dplyr::group_modify(~{
    df <- .x %>% dplyr::arrange(Time)
    if (nrow(df) < 5) {
      return(tibble(O2_ref = NA_real_))
    }
    
    df <- df %>%
      dplyr::mutate(
        dO2   = c(NA, diff(Oxygen)),
        dt    = c(NA, diff(Time)),
        dO2dt = dO2 / dt,
        sm    = zoo::rollmean(dO2dt, 3, fill = NA, align = "right")
      )
    
    idx <- which(df$sm < -1e-7)[1]
    idx <- if (is.na(idx) || idx > nrow(df) - 2) 1 else idx + 5
    df_trim <- df[idx:nrow(df), ]
    
    if (nrow(df_trim) < 3) {
      return(tibble(O2_ref = NA_real_))
    }
    
    O0 <- mean(head(df_trim$Oxygen, 3), na.rm = TRUE)
    tibble(O2_ref = O0)
  }) %>%
  dplyr::ungroup()

write_csv(
  O2_ref_series,
  file.path(tables_dir, "O2_ref_series_newmodel_MC.csv")
)

print(head(O2_ref_series))

# ───────────────── 4. Build series-level info for MC ────────────────────── #

# series_info: one row per Taxon × Replicate with:
#   - resp_tot (from model)
#   - O2_ref (from oxygen data)
#   - N0_mean, N0_sd (from taxon-level N0 stats)

series_info <- model_results %>%
  dplyr::select(Taxon, Replicate, resp_tot) %>%
  dplyr::left_join(
    O2_ref_series,
    by = c("Taxon", "Replicate")
  ) %>%
  dplyr::left_join(
    N0_stats_taxon,
    by = "Taxon"
  ) %>%
  dplyr::filter(
    is.finite(resp_tot),
    is.finite(O2_ref),
    is.finite(N0_mean),
    N0_mean > 0
  )

write_csv(
  series_info,
  file.path(tables_dir, "N0_MC_newmodel_series_info.csv")
)

print(head(series_info))

# ───────────────── 5. Monte Carlo settings ──────────────────────────────── #

N_MC <- 300  # number of N0 draws per Taxon × Replicate

set.seed(123)  # reproducible

# ───────────────── 6. Monte Carlo: N0 uncertainty → R uncertainty ───────── #

# We will treat each row of series_info as a parameter combination:
#   (Taxon, Replicate, resp_tot, O2_ref, N0_mean, N0_sd).
#
# For each row:
#   - Draw N0_draw ~ Normal(N0_mean, N0_sd), truncated > 0
#   - Compute R_draw = resp_tot * O2_ref / N0_draw
#
# This directly encodes how plausible variability in N0
# affects the per-cell respiration rate R, while
# holding resp_tot and O2_ref fixed to their fitted/observed values.

mc_results <- purrr::pmap_dfr(
  series_info,
  function(Taxon, Replicate, resp_tot, O2_ref, N0_mean, N0_sd, n_N0) {
    
    # if N0_sd is missing or zero, give small relative SD so MC can run
    if (!is.finite(N0_sd) || N0_sd <= 0) {
      N0_sd <- 0.05 * N0_mean
    }
    
    purrr::map_dfr(seq_len(N_MC), function(mc_id) {
      # Draw N0 from Normal(mean, sd), truncated at > 0
      N0_draw <- NA_real_
      tries <- 0
      while (is.na(N0_draw) || N0_draw <= 0) {
        N0_draw <- rnorm(1, mean = N0_mean, sd = N0_sd)
        tries <- tries + 1
        if (tries > 50) {
          N0_draw <- N0_mean
          break
        }
      }
      
      R_draw <- resp_tot * O2_ref / N0_draw  # mg O2 cell^-1 min^-1
      
      tibble::tibble(
        Taxon     = Taxon,
        Replicate = Replicate,
        mc_id     = mc_id,
        resp_tot  = resp_tot,
        O2_ref    = O2_ref,
        N0_draw   = N0_draw,
        R_draw    = R_draw
      )
    })
  }
)

write_csv(
  mc_results,
  file.path(tables_dir, "N0_MC_newmodel_raw_R_draws.csv")
)

print(head(mc_results))

# ───────────────── 7. Summarise uncertainty in R per series ─────────────── #

mc_summary <- mc_results %>%
  dplyr::filter(is.finite(R_draw), R_draw > 0) %>%
  dplyr::group_by(Taxon, Replicate) %>%
  dplyr::summarise(
    n_MC        = dplyr::n(),
    N0_mean_MC  = mean(N0_draw, na.rm = TRUE),
    N0_sd_MC    = sd(N0_draw,   na.rm = TRUE),
    R_mean      = mean(R_draw,  na.rm = TRUE),
    R_sd        = sd(R_draw,    na.rm = TRUE),
    R_q025      = stats::quantile(R_draw, 0.025, na.rm = TRUE),
    R_q975      = stats::quantile(R_draw, 0.975, na.rm = TRUE),
    R_rel_sd    = R_sd / R_mean,  # relative SD (dimensionless)
    .groups     = "drop"
  )

write_csv(
  mc_summary,
  file.path(tables_dir, "N0_MC_newmodel_summary_per_series.csv")
)

print(mc_summary)

# ───────────────── 8. Taxon-level summary ───────────────────────────────── #

mc_summary_taxon <- mc_summary %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(
    n_reps        = dplyr::n(),
    R_rel_sd_med  = median(R_rel_sd, na.rm = TRUE),
    R_rel_sd_min  = min(R_rel_sd, na.rm = TRUE),
    R_rel_sd_max  = max(R_rel_sd, na.rm = TRUE),
    .groups       = "drop"
  )

write_csv(
  mc_summary_taxon,
  file.path(tables_dir, "Table_S1_N0_MC_newmodel_taxon_summary.csv")
)

print(mc_summary_taxon)

# ───────────────── 9. Global summary for Methods text ───────────────────── #

overall_stats <- mc_summary %>%
  dplyr::summarise(
    n_series    = dplyr::n(),
    R_rel_sd_med = median(R_rel_sd, na.rm = TRUE),
    R_rel_sd_min = min(R_rel_sd, na.rm = TRUE),
    R_rel_sd_max = max(R_rel_sd, na.rm = TRUE)
  )

print(overall_stats)

# ───────────────── 10. Diagnostic plots: relative SD of R ───────────────── #

mc_long <- mc_summary %>%
  dplyr::select(Taxon, Replicate, R_rel_sd)

# Density of relative SD across all series
p_rel_sd <- ggplot(mc_long, aes(x = R_rel_sd)) +
  ggplot2::geom_density(fill = "steelblue", alpha = 0.4) +
  ggplot2::scale_x_continuous(labels = scales::scientific) +
  ggplot2::labs(
    x = "Relative standard deviation of R (SD / mean)",
    y = "Density",
    title = "MC sensitivity of per-cell respiration rate R to N0 (new model)"
  ) +
  ggplot2::theme_classic(base_size = 12)

ggplot2::ggsave(
  file.path(plots_dir, "N0_MC_newmodel_R_rel_sd_density.png"),
  p_rel_sd, width = 6, height = 4, dpi = 300
)
ggplot2::ggsave(
  file.path(plots_dir, "N0_MC_newmodel_R_rel_sd_density.pdf"),
  p_rel_sd, width = 6, height = 4
)

# Boxplot of R_rel_sd by taxon
p_box <- ggplot(mc_long, aes(x = Taxon, y = R_rel_sd)) +
  ggplot2::geom_boxplot(outlier.size = 0.8, fill = "grey70") +
  ggplot2::scale_y_continuous(labels = scales::scientific) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = "Taxon",
    y = "Relative standard deviation of R (SD / mean)",
    title = "MC relative uncertainty in R by taxon (N0 sensitivity, new model)"
  ) +
  ggplot2::theme_classic(base_size = 11)

ggplot2::ggsave(
  file.path(plots_dir, "N0_MC_newmodel_R_rel_sd_by_taxon.png"),
  p_box, width = 7, height = 6, dpi = 300
)
ggplot2::ggsave(
  file.path(plots_dir, "N0_MC_newmodel_R_rel_sd_by_taxon.pdf"),
  p_box, width = 7, height = 6
)

###############################################################################
# End of Monte Carlo sensitivity script (new model, R vs N0)
###############################################################################
