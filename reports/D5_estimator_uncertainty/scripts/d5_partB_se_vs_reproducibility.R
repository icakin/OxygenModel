# =============================================================================
# D5 PART B - the naive SE against actual replicate reproducibility
# =============================================================================
# Per taxon: the median within-curve relative SE from summary(nls) against the
# between-replicate sd(log r) and sd(log K) actually observed across that taxon's
# five replicates. Then the same with PART A's autocorrelation-corrected SEs, to
# see how much of the gap the correction closes.
#
# A within-curve relative SE se_r/r IS the sd of log r to first order, so the two
# are directly comparable on the log scale.
#
# INPUT   runs/D5_analysis/A1_per_curve_residual_diagnostics.csv
# OUTPUT  runs/D5_analysis/B1_*.csv .. B4_*.csv
# =============================================================================

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d5_common.R"))

message("\n== D5 PART B: within-curve SE vs between-replicate spread ==")

A1 <- readr::read_csv(d5f("A1_per_curve_residual_diagnostics.csv"), show_col_types = FALSE)

# Grouping for "replicates of the same thing": everything in the key except
# Replicate. On a temperature dataset that is Taxon x Temperature automatically.
ID    <- id_cols(A1)
GRP   <- setdiff(ID, "Replicate")
message("  replicate groups defined by: ", paste(GRP, collapse = " x "))

B1 <- A1 |>
  dplyr::group_by(dplyr::across(dplyr::all_of(GRP))) |>
  dplyr::summarise(
    n_rep = dplyr::n(),
    # within-curve: median relative SE across the group's curves
    within_se_logr_naive = stats::median(se_r_naive / r, na.rm = TRUE),
    within_se_logK_naive = stats::median(se_K_naive / K, na.rm = TRUE),
    within_se_logr_corr  = stats::median(se_r_gls / r, na.rm = TRUE),
    within_se_logK_corr  = stats::median(se_K_gls / K, na.rm = TRUE),
    # between-replicate: the sd actually observed
    between_sd_logr = stats::sd(log(r)),
    between_sd_logK = stats::sd(log(K)),
    median_r = stats::median(r), median_K = stats::median(K),
    .groups = "drop") |>
  dplyr::mutate(
    ratio_logr_naive = between_sd_logr / within_se_logr_naive,
    ratio_logK_naive = between_sd_logK / within_se_logK_naive,
    ratio_logr_corr  = between_sd_logr / within_se_logr_corr,
    ratio_logK_corr  = between_sd_logK / within_se_logK_corr) |>
  dplyr::arrange(dplyr::desc(ratio_logr_naive))

d5_write(B1, "B1_within_se_vs_between_replicate_by_taxon.csv")

med <- function(x) stats::median(x, na.rm = TRUE)
B2 <- tibble::tibble(
  quantity = c("median within-curve se(log r) naive",
               "median within-curve se(log r) AR(1)-corrected",
               "median between-replicate sd(log r)",
               "RATIO between/within naive (log r)",
               "RATIO between/within corrected (log r)",
               "median within-curve se(log K) naive",
               "median within-curve se(log K) AR(1)-corrected",
               "median between-replicate sd(log K)",
               "RATIO between/within naive (log K)",
               "RATIO between/within corrected (log K)"),
  value = c(med(B1$within_se_logr_naive), med(B1$within_se_logr_corr),
            med(B1$between_sd_logr), med(B1$ratio_logr_naive), med(B1$ratio_logr_corr),
            med(B1$within_se_logK_naive), med(B1$within_se_logK_corr),
            med(B1$between_sd_logK), med(B1$ratio_logK_naive), med(B1$ratio_logK_corr)))
d5_write(B2, "B2_se_vs_reproducibility_overall.csv")

# ---- how much of the gap does the correction close? -------------------------
# Gap measured on the log scale: log(between/within). The correction closes
# log(within_corr/within_naive) of it.
gap_close <- function(within_naive, within_corr, between) {
  gap0 <- log(between / within_naive)
  gap1 <- log(between / within_corr)
  ifelse(gap0 == 0, NA_real_, (gap0 - gap1) / gap0)
}
B3 <- tibble::tibble(
  parameter = c("log r", "log K"),
  median_gap_naive_x = c(med(B1$ratio_logr_naive), med(B1$ratio_logK_naive)),
  median_gap_corrected_x = c(med(B1$ratio_logr_corr), med(B1$ratio_logK_corr)),
  median_fraction_of_gap_closed =
    c(med(gap_close(B1$within_se_logr_naive, B1$within_se_logr_corr, B1$between_sd_logr)),
      med(gap_close(B1$within_se_logK_naive, B1$within_se_logK_corr, B1$between_sd_logK))),
  interpretation = "fraction of the log-scale gap removed by the AR(1) correction")
d5_write(B3, "B3_gap_closed_by_correction.csv")

# ---- the variance split ------------------------------------------------------
# Total replicate-to-replicate variance = measurement (within-curve) + genuine
# biological/handling variation. Since between_sd is the TOTAL observed spread,
# the residual component is sqrt(max(between^2 - within^2, 0)).
B4 <- B1 |>
  dplyr::transmute(
    dplyr::across(dplyr::all_of(GRP)),
    total_sd_logr = between_sd_logr,
    measurement_sd_logr = within_se_logr_naive,
    biological_sd_logr = sqrt(pmax(between_sd_logr^2 - within_se_logr_naive^2, 0)),
    pct_variance_measurement =
      100 * within_se_logr_naive^2 / between_sd_logr^2,
    total_sd_logK = between_sd_logK,
    measurement_sd_logK = within_se_logK_naive,
    biological_sd_logK = sqrt(pmax(between_sd_logK^2 - within_se_logK_naive^2, 0)),
    pct_variance_measurement_K =
      100 * within_se_logK_naive^2 / between_sd_logK^2) |>
  dplyr::arrange(dplyr::desc(pct_variance_measurement))
d5_write(B4, "B4_variance_split_by_taxon.csv")

message(sprintf(
  "\n  log r: within-curve se median %.4f | between-replicate sd median %.4f | ratio %.1fx",
  med(B1$within_se_logr_naive), med(B1$between_sd_logr), med(B1$ratio_logr_naive)))
message(sprintf(
  "  log K: within-curve se median %.4f | between-replicate sd median %.4f | ratio %.1fx",
  med(B1$within_se_logK_naive), med(B1$between_sd_logK), med(B1$ratio_logK_naive)))
message(sprintf(
  "  measurement noise accounts for a median of %.2f%% of replicate variance in log r",
  med(B4$pct_variance_measurement)))
message(sprintf(
  "  AR(1) correction closes %.1f%% of the log r gap",
  100 * med(gap_close(B1$within_se_logr_naive, B1$within_se_logr_corr, B1$between_sd_logr))))
message("D5 PART B done.")
