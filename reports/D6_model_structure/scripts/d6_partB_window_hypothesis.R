# =============================================================================
# D6 PART B - is the window sensitivity a symptom of saturation?
# =============================================================================
# THE PIVOT OF THIS STUDY, and it is a hypothesis to test, not to confirm.
#
# 08_window_sensitivity.R found that a +/-6 min shift of the window edges moves
# K and per-cell R by 15-34% while moving r by ~0.3%. The hypothesis is that if
# growth or oxygen uptake saturates near the end of some windows, moving the end
# changes how much saturation is included - producing exactly that asymmetry.
#
# If true: per-curve window sensitivity should correlate with the PART A
# saturation measures, and the fix is a model term rather than tighter curation.
# If false: the two findings are independent and both need addressing separately.
#
# A WEAK CORRELATION IS A REAL ANSWER and is reported as one.
#
# INPUT   results/tables/window_sensitivity_percurve.csv (from 08, unchanged)
#         runs/D6_analysis/A1_residual_structure_per_curve.csv
# OUTPUT  runs/D6_analysis/B1_*.csv .. B4_*.csv
# =============================================================================

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d6_common.R"))

message("\n== D6 PART B: window sensitivity vs saturation measures ==")

ws <- readr::read_csv(file.path(tables_dir, "window_sensitivity_percurve.csv"),
                      show_col_types = FALSE)
A1 <- readr::read_csv(d6f("A1_residual_structure_per_curve.csv"), show_col_types = FALSE)
ID <- id_cols(A1)

# Per-curve sensitivity: the largest absolute move across the shift grid.
B1 <- ws |>
  dplyr::filter(!is_base, ok) |>
  dplyr::group_by(dplyr::across(dplyr::all_of(ID))) |>
  dplyr::summarise(
    n_shifts = dplyr::n(),
    max_abs_dK_pct = max(abs(dK_pct), na.rm = TRUE),
    max_abs_dR_pct = max(abs(dR_pct), na.rm = TRUE),
    max_abs_dr_pct = max(abs(dr_pct), na.rm = TRUE),
    .groups = "drop") |>
  dplyr::inner_join(A1, by = ID)
d6_write(B1, "B1_window_sensitivity_vs_saturation.csv")

# Correlations, Pearson and Spearman (the latter in case of a monotone but
# nonlinear relation), with a p-value so "weak" can be stated with a number.
ct <- function(x, y, method) {
  o <- suppressWarnings(stats::cor.test(x, y, method = method))
  c(est = unname(o$estimate), p = unname(o$p.value))
}
mk <- function(pred_name, pred, mechanism) {
  a <- ct(pred, B1$max_abs_dK_pct, "pearson")
  b <- ct(pred, B1$max_abs_dK_pct, "spearman")
  c_ <- ct(pred, B1$max_abs_dR_pct, "pearson")
  d_ <- ct(pred, B1$max_abs_dR_pct, "spearman")
  tibble::tibble(
    predictor = pred_name, mechanism = mechanism,
    r_pearson_dK = a[1], p_pearson_dK = a[2],
    rho_spearman_dK = b[1], p_spearman_dK = b[2],
    r_pearson_dR = c_[1], p_pearson_dR = c_[2],
    rho_spearman_dR = d_[1], p_spearman_dR = d_[2])
}

B2 <- dplyr::bind_rows(
  mk("late-window mean residual (sd)",   B1$late_mean_sd,       "saturation signature"),
  mk("residual curvature",               B1$curvature,          "saturation signature"),
  mk("RSS gain from growth-sat term",    B1$rss_gain_growth,    "growth saturation"),
  mk("RSS gain from O2-limitation term", B1$rss_gain_o2,        "O2 limitation"),
  mk("min absolute O2 reached",          B1$min_O2_abs,         "O2 limitation"),
  mk("fraction O2 remaining at end",     B1$O2_frac_remaining,  "O2 limitation"),
  mk("elapsed growth r*t at end",        B1$growth_at_end,      "growth saturation"),
  mk("N-fold increase at end",           B1$Nfold_at_end,       "growth saturation"),
  # controls: things that are NOT saturation, to show what does predict it
  mk("window length (min)",              B1$t_end,              "control"),
  mk("residual sd (noise level)",        B1$sigma,              "control"),
  mk("fitted r",                         B1$r,                  "control"),
  mk("fitted K",                         B1$K,                  "control"))
d6_write(B2, "B2_correlations.csv")

# How much of the K/R movement does saturation explain? R^2 of the best
# saturation predictor, against the best control.
r2 <- function(y, x) {
  m <- stats::lm(y ~ x); summary(m)$r.squared
}
sat_preds <- list(`late-window mean residual` = B1$late_mean_sd,
                  `residual curvature` = B1$curvature,
                  `RSS gain growth-sat` = B1$rss_gain_growth,
                  `elapsed growth at end` = B1$growth_at_end,
                  `fraction O2 remaining` = B1$O2_frac_remaining)
ctl_preds <- list(`window length` = B1$t_end,
                  `residual sd` = B1$sigma,
                  `fitted r` = B1$r)
B3 <- dplyr::bind_rows(
  tibble::tibble(predictor = names(sat_preds), family = "saturation",
                 R2_dK = vapply(sat_preds, function(x) r2(B1$max_abs_dK_pct, x), 0),
                 R2_dR = vapply(sat_preds, function(x) r2(B1$max_abs_dR_pct, x), 0)),
  tibble::tibble(predictor = names(ctl_preds), family = "control",
                 R2_dK = vapply(ctl_preds, function(x) r2(B1$max_abs_dK_pct, x), 0),
                 R2_dR = vapply(ctl_preds, function(x) r2(B1$max_abs_dR_pct, x), 0))
) |> dplyr::arrange(dplyr::desc(R2_dR))
d6_write(B3, "B3_variance_explained.csv")

# A multiple regression of the sensitivity on ALL saturation measures at once,
# so the verdict does not rest on any single predictor.
mfull <- stats::lm(max_abs_dR_pct ~ late_mean_sd + curvature + rss_gain_growth +
                     rss_gain_o2 + O2_frac_remaining + growth_at_end, data = B1)
B4 <- tibble::tibble(
  quantity = c("R^2, all saturation measures -> max|dR_pct|",
               "adjusted R^2",
               "model p-value",
               "R^2, window length alone -> max|dR_pct|",
               "median max|dK_pct| across curves",
               "median max|dR_pct| across curves",
               "median max|dr_pct| across curves"),
  value = c(summary(mfull)$r.squared,
            summary(mfull)$adj.r.squared,
            stats::pf(summary(mfull)$fstatistic[1], summary(mfull)$fstatistic[2],
                      summary(mfull)$fstatistic[3], lower.tail = FALSE),
            r2(B1$max_abs_dR_pct, B1$t_end),
            stats::median(B1$max_abs_dK_pct), stats::median(B1$max_abs_dR_pct),
            stats::median(B1$max_abs_dr_pct)))
d6_write(B4, "B4_headline.csv")

# ---- the alternative explanation, and it is arithmetic ----------------------
# Shifting a window edge by Delta on a curve growing at rate r changes the
# oxygen consumed over the window by roughly exp(r*Delta), so the fitted K must
# move by about r*Delta to compensate. That is a property of the exponential,
# not a symptom of misspecification. If true, max|dK_pct| ~ 100*r*Delta with a
# coefficient near 1, and r alone should explain nearly all the variance.
DELTA <- 6
arith <- stats::lm(max_abs_dK_pct ~ 0 + I(100 * r * DELTA), data = B1)
B5 <- tibble::tibble(
  quantity = c("corr(max|dK_pct|, fitted r)",
               "R^2 of fitted r alone on max|dK_pct|",
               "median fitted r (per min)",
               "arithmetic prediction 100*r*6 (%)",
               "observed median max|dK_pct| (%)",
               "regression coefficient vs the arithmetic prediction (1.0 = exact)",
               "R^2 of ALL saturation measures on max|dR_pct|"),
  value = c(stats::cor(B1$max_abs_dK_pct, B1$r),
            stats::cor(B1$max_abs_dK_pct, B1$r)^2,
            stats::median(B1$r),
            100 * stats::median(B1$r) * DELTA,
            stats::median(B1$max_abs_dK_pct),
            unname(stats::coef(arith)[1]),
            summary(mfull)$r.squared))
d6_write(B5, "B5_arithmetic_explanation.csv")

message(sprintf("\n  ARITHMETIC: corr(max|dK|, r) = %.4f (R^2 %.3f); 100*r*6 = %.1f%% vs observed %.1f%%",
  stats::cor(B1$max_abs_dK_pct, B1$r), stats::cor(B1$max_abs_dK_pct, B1$r)^2,
  100 * stats::median(B1$r) * DELTA, stats::median(B1$max_abs_dK_pct)))

message(sprintf("\n  median per-curve sensitivity: |dK| %.1f%%, |dR| %.1f%%, |dr| %.2f%%",
  stats::median(B1$max_abs_dK_pct), stats::median(B1$max_abs_dR_pct),
  stats::median(B1$max_abs_dr_pct)))
message(sprintf("  R^2 of ALL saturation measures on |dR_pct|: %.3f (adj %.3f, p = %.3g)",
  summary(mfull)$r.squared, summary(mfull)$adj.r.squared,
  stats::pf(summary(mfull)$fstatistic[1], summary(mfull)$fstatistic[2],
            summary(mfull)$fstatistic[3], lower.tail = FALSE)))
message("  strongest single predictors of |dR_pct|:")
print(as.data.frame(head(B3, 5)), row.names = FALSE, digits = 3)
message("D6 PART B done.")
