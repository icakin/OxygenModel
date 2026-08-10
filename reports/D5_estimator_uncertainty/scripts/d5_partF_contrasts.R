# =============================================================================
# D5 PART F - what the corrected uncertainty does to the taxon-level contrasts
# =============================================================================
# The pipeline reports taxon-level growth and respiration. The question a reader
# cares about is which taxa are actually distinguishable. All 105 pairwise
# contrasts among the 15 taxa are recomputed under three uncertainty models:
#
#   within-curve only   the two-stage pipeline's own se_r/r, propagated to the
#                       taxon mean as se/sqrt(n) - i.e. pretending replicate
#                       spread does not exist
#   AR(1)-corrected     the same with PART A's corrected SEs
#   joint posterior     the taxon-level posterior SDs from PART E (naive S_i),
#                       which DO include between-replicate variance
#
# A contrast counts as resolved if its 95% interval excludes zero, with
# Holm correction across the 105 comparisons.
#
# OUTPUT  runs/D5_analysis/F1_*.csv .. F3_*.csv
# =============================================================================

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d5_common.R"))

message("\n== D5 PART F: taxon contrasts under each uncertainty model ==")

A1 <- readr::read_csv(d5f("A1_per_curve_residual_diagnostics.csv"), show_col_types = FALSE)
E1 <- readr::read_csv(d5f("E1_joint_posterior_by_variant.csv"),     show_col_types = FALSE)
ID  <- id_cols(A1); GRP <- setdiff(ID, "Replicate")

# Taxon means and the three standard errors of that mean.
tx <- A1 |>
  dplyr::group_by(dplyr::across(dplyr::all_of(GRP))) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_logr = mean(log(r)),
    # (1) within-curve only: average per-curve se(log r), then /sqrt(n)
    se_within_naive = mean(se_r_naive / r) / sqrt(dplyr::n()),
    se_within_corr  = mean(se_r_gls   / r) / sqrt(dplyr::n()),
    .groups = "drop")

jt <- E1 |> dplyr::filter(variant == "naive") |>
  dplyr::select(Taxon, mean_logr_joint = m_lr, se_joint = sd_lr)
tx <- tx |> dplyr::left_join(jt, by = "Taxon")
d5_write(tx, "F1_taxon_level_uncertainty.csv")

# All pairwise contrasts.
pairs <- utils::combn(nrow(tx), 2)
mk <- function(mean_col, se_col, label) {
  i <- pairs[1, ]; j <- pairs[2, ]
  d  <- tx[[mean_col]][i] - tx[[mean_col]][j]
  se <- sqrt(tx[[se_col]][i]^2 + tx[[se_col]][j]^2)
  z  <- d / se
  p  <- 2 * stats::pnorm(-abs(z))
  tibble::tibble(
    model = label,
    taxon_a = tx$Taxon[i], taxon_b = tx$Taxon[j],
    diff_logr = d, se = se, z = z, p = p,
    p_holm = stats::p.adjust(p, method = "holm"),
    resolved_raw  = p < 0.05,
    resolved_holm = stats::p.adjust(p, method = "holm") < 0.05)
}

F2 <- dplyr::bind_rows(
  mk("mean_logr", "se_within_naive", "within-curve only (naive)"),
  mk("mean_logr", "se_within_corr",  "within-curve only (AR(1)-corrected)"),
  mk("mean_logr_joint", "se_joint",  "joint posterior (includes replicate spread)")
)
d5_write(F2, "F2_pairwise_contrasts.csv")

F3 <- F2 |> dplyr::group_by(model) |>
  dplyr::summarise(
    n_contrasts = dplyr::n(),
    resolved_uncorrected = sum(resolved_raw),
    resolved_holm = sum(resolved_holm),
    pct_resolved_holm = 100 * sum(resolved_holm) / dplyr::n(),
    median_se = stats::median(se),
    .groups = "drop")
d5_write(F3, "F3_contrast_counts.csv")

print(as.data.frame(F3), row.names = FALSE, digits = 4)
message("D5 PART F done.")
