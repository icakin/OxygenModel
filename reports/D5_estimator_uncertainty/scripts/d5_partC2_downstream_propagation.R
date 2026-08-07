# =============================================================================
# D5 PART C (continued) - which rotated combination do R and CUE actually use?
# =============================================================================
# The brief asserted that per-cell R and CUE depend on the POORLY determined
# combination log K - log r. Deriving it from the pipeline instead:
#
#   C_tot            = (K/r)(exp(r T_end) - 1) * O2_ref              [05, l.170]
#   biomass_integral = N0 (exp(r T_end) - 1) / r                     [05, l.246]
#   R = C_tot / biomass_integral
#     = [(K/r)(e^{rT}-1) O2ref] / [N0 (e^{rT}-1)/r]
#     = K * O2_ref / N0
#
# The (e^{rT}-1)/r factors CANCEL EXACTLY. r does not enter R through the
# integral at all. It re-enters only through the depletion anchor:
#
#   N0 = FC_Final * FC_TO_CELLS_PER_L * exp(-r (t_depletion - fit_start))
#   =>  log R = log K + r * dt + const,     dt = t_depletion - fit_start
#
# So the combination that matters is log K + dt*r, NOT log K - log r. Whether
# that combination is well or badly determined is an empirical question, and the
# per-curve vcov answers it directly by the delta method:
#
#   var(log R) = var(log K) + dt^2 var(r) + 2 dt cov(log K, r)
#
# The third term is where the r-K ridge either helps or hurts. Comparing against
# the same expression with cov set to zero isolates its contribution.
#
# Growth carbon G = r * 60 * cell_carbon, so log G = log r + const and
# var(log G) = var(r)/r^2 exactly.
#
# OUTPUT  runs/D5_analysis/C6_*.csv, C7_*.csv
# =============================================================================

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d5_common.R"))

message("\n== D5 PART C (cont.): propagation into R, G and CUE ==")

curves <- readr::read_csv(FITCURVES_CSV,     show_col_types = FALSE)
res    <- readr::read_csv(RESULTS_FINAL_CSV, show_col_types = FALSE)
ID     <- id_cols(res)

rows <- vector("list", nrow(res))
for (i in seq_len(nrow(res))) {
  key <- res[i, ID, drop = FALSE]
  g <- curves |> dplyr::semi_join(key, by = ID) |> dplyr::arrange(Time0_min)
  if (nrow(g) < 10) next
  r0 <- res$r_per_minute[i]; K0 <- res$K[i]; O0 <- res$O2_0[i]
  dt <- res$t_depletion_min[i] - res$fit_start_min[i]
  f <- fit_ols(g$Time0_min, g$Oxygen_norm, list(r = r0, K = K0, O2_0 = O0))
  if (inherits(f, "try-error")) next
  V <- tryCatch(stats::vcov(f), error = function(e) NULL)
  if (is.null(V)) next

  var_logK  <- V["K","K"] / K0^2
  var_r     <- V["r","r"]
  cov_logK_r<- V["r","K"] / K0
  var_logG  <- V["r","r"] / r0^2

  var_logR      <- var_logK + dt^2 * var_r + 2 * dt * cov_logK_r
  var_logR_indep<- var_logK + dt^2 * var_r          # covariance term removed

  rows[[i]] <- tibble::tibble(
    key, dt_min = dt, r = r0, K = K0,
    R_C_fg_cell_h = res$R_C_fg_cell_h[i], G_C_fg_cell_h = res$G_C_fg_cell_h[i],
    sd_logK = sqrt(max(var_logK, 0)), sd_logG = sqrt(max(var_logG, 0)),
    term_logK = var_logK, term_dt2_var_r = dt^2 * var_r,
    term_cross = 2 * dt * cov_logK_r,
    sd_logR = sqrt(max(var_logR, 0)),
    sd_logR_if_independent = sqrt(max(var_logR_indep, 0)))
}

C6 <- dplyr::bind_rows(rows) |>
  dplyr::mutate(
    ratio_actual_over_independent = sd_logR / sd_logR_if_independent,
    CUE = G_C_fg_cell_h / (G_C_fg_cell_h + R_C_fg_cell_h))
d5_write(C6, "C6_R_propagation_per_curve.csv")

GRP <- setdiff(ID, "Replicate")
obs <- res |>
  dplyr::filter(is.finite(R_C_fg_cell_h), R_C_fg_cell_h > 0) |>
  dplyr::group_by(dplyr::across(dplyr::all_of(GRP))) |>
  dplyr::summarise(observed_sd_logR_replicates = stats::sd(log(R_C_fg_cell_h)),
                   .groups = "drop")

C7 <- C6 |>
  dplyr::group_by(dplyr::across(dplyr::all_of(GRP))) |>
  dplyr::summarise(
    median_dt_min = stats::median(dt_min),
    within_sd_logR = stats::median(sd_logR),
    within_sd_logR_if_indep = stats::median(sd_logR_if_independent),
    within_sd_logG = stats::median(sd_logG),
    .groups = "drop") |>
  dplyr::left_join(obs, by = GRP) |>
  dplyr::mutate(ratio_observed_over_within = observed_sd_logR_replicates / within_sd_logR)
d5_write(C7, "C7_R_propagation_by_taxon.csv")

md <- function(x) stats::median(x, na.rm = TRUE)
C8 <- tibble::tibble(
  quantity = c("median dt = t_depletion - fit_start (min)",
               "median within-curve sd(log R)",
               "median within-curve sd(log R) if r and K were independent",
               "ratio actual / independent",
               "median within-curve sd(log G) = sd(log r)",
               "median observed sd(log R) across replicates",
               "ratio observed / within-curve",
               "median share of var(log R) from var(log K)",
               "median share of var(log R) from dt^2 var(r)",
               "median share of var(log R) from the 2 dt cov term"),
  value = c(md(C6$dt_min), md(C6$sd_logR), md(C6$sd_logR_if_independent),
            md(C6$ratio_actual_over_independent), md(C6$sd_logG),
            md(C7$observed_sd_logR_replicates), md(C7$ratio_observed_over_within),
            md(C6$term_logK / (C6$sd_logR^2)),
            md(C6$term_dt2_var_r / (C6$sd_logR^2)),
            md(C6$term_cross / (C6$sd_logR^2))))
d5_write(C8, "C8_R_propagation_headline.csv")

message(sprintf("\n  dt median %.0f min | within-curve sd(log R) %.4f | if independent %.4f (ratio %.2f)",
  md(C6$dt_min), md(C6$sd_logR), md(C6$sd_logR_if_independent),
  md(C6$ratio_actual_over_independent)))
message(sprintf("  observed sd(log R) across replicates %.4f -> %.0fx the within-curve value",
  md(C7$observed_sd_logR_replicates), md(C7$ratio_observed_over_within)))
message("D5 PART C (cont.) done.")
