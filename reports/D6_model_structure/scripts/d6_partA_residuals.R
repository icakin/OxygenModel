# =============================================================================
# D6 PART A - is saturation present in the real curves?
# =============================================================================
# Fits the current model over each committed window and examines the residuals
# as a function of (i) time within the window and (ii) ABSOLUTE oxygen.
#
# The two candidate mechanisms make DIFFERENT predictions, which is what makes
# them separable:
#
#   OXYGEN LIMITATION   deviation scales with absolute O2, appears in any curve
#                       that gets low enough, regardless of taxon
#   GROWTH SATURATION   deviation scales with elapsed growth (r*t), not with O2
#
# So the discriminating test is: does the late-window residual trend track the
# minimum O2 reached, or the elapsed growth at the window end? Both are computed
# per curve and correlated against the residual curvature.
#
# INPUT   results/tables/oxygen_fit_curves.csv, oxygen_fit_results.csv,
#         oxygen_results_with_R.csv
# OUTPUT  runs/D6_analysis/A1_*.csv .. A5_*.csv
# =============================================================================

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d6_common.R"))

message("\n== D6 PART A: residual structure in the real curves ==")

curves <- readr::read_csv(FITCURVES_CSV,     show_col_types = FALSE)
fits   <- readr::read_csv(RESULTS_FIT_CSV,   show_col_types = FALSE)
res    <- readr::read_csv(RESULTS_FINAL_CSV, show_col_types = FALSE)
ID     <- id_cols(fits)

rows <- list(); prof <- list()

for (i in seq_len(nrow(fits))) {
  key <- fits[i, ID, drop = FALSE]
  g <- curves |> dplyr::semi_join(key, by = ID) |> dplyr::arrange(Time0_min)
  if (nrow(g) < 15) next
  t <- g$Time0_min; y <- g$Oxygen_norm; n <- length(t)
  start <- list(r = fits$r_per_minute[i], K = fits$K[i], O2_0 = fits$O2_0[i])
  f <- fit_m0(t, y, start)
  if (inherits(f, "try-error")) next
  cf <- stats::coef(f); e <- y - stats::predict(f)
  O2ref <- fits$O2_ref[i]

  # absolute oxygen (mg/L) on the normalised trace
  o2_abs <- y * O2ref
  # elapsed growth at each point, and at the window end
  growth <- cf[["r"]] * t

  # last third of the window: the signature region
  cut <- t >= (min(t) + (2 / 3) * (max(t) - min(t)))
  lt  <- if (sum(cut) >= 5) stats::coef(stats::lm(e[cut] ~ t[cut]))[2] else NA_real_
  # residual mean in the last third, in units of the residual sd (effect size)
  sig <- stats::sd(e)
  late_mean_sd <- if (sum(cut) >= 5) mean(e[cut]) / sig else NA_real_

  # curvature: quadratic term of residuals vs time (a systematic bend)
  qq <- try(stats::coef(stats::lm(e ~ poly(t, 2, raw = TRUE))), silent = TRUE)
  curv <- if (inherits(qq, "try-error")) NA_real_ else unname(qq[3])

  # DEVIATION TEST. A RESET-style regression of the residuals on the squared
  # fitted value has no power here: for this 3-parameter model fitted^2 is very
  # nearly in the span of the Jacobian, so the residuals are orthogonal to it by
  # construction and every p-value comes back ~1. It was tried and discarded.
  #
  # The valid test is a NESTED F-test against each saturating alternative, since
  # M0 is M2g at Nmax_rel -> Inf and M2o at Km -> 0.
  rss0 <- sum(e^2)
  d <- data.frame(t = t, y = y)
  ftest <- function(fit_alt, k_extra) {
    if (inherits(fit_alt, "try-error") || is.null(fit_alt)) return(c(NA_real_, NA_real_))
    rss1 <- sum(stats::residuals(fit_alt)^2)
    if (!is.finite(rss1) || rss1 <= 0 || rss1 >= rss0) return(c(1, (rss0 - rss1) / rss0))
    df2 <- n - (3 + k_extra)
    Fst <- ((rss0 - rss1) / k_extra) / (rss1 / df2)
    c(stats::pf(Fst, k_extra, df2, lower.tail = FALSE), (rss0 - rss1) / rss0)
  }
  f2g <- fit_try(y ~ m2g_model(r, K, t, O2_0, Nmax_rel), d,
                 start = c(as.list(cf), Nmax_rel = 50),
                 lower = c(FIT_LOWER, Nmax_rel = 1.01),
                 upper = c(FIT_UPPER, Nmax_rel = 1e6))
  f2o <- fit_try(y ~ m2o_predict(r, K, t, O2_0, Km), d,
                 start = c(as.list(cf), Km = 0.05),
                 lower = c(FIT_LOWER, Km = 1e-6),
                 upper = c(FIT_UPPER, Km = 10))
  g_res <- ftest(f2g, 1); o_res <- ftest(f2o, 1)
  p_sat_growth <- g_res[1]; rss_gain_growth <- g_res[2]
  p_sat_o2     <- o_res[1]; rss_gain_o2     <- o_res[2]

  rows[[length(rows) + 1]] <- tibble::tibble(
    key, n = n, sigma = sig,
    r = cf[["r"]], K = cf[["K"]], O2_0 = cf[["O2_0"]], O2_ref = O2ref,
    t_end = max(t), min_O2_norm = min(y), min_O2_abs = min(o2_abs),
    O2_frac_remaining = min(y) / cf[["O2_0"]],
    growth_at_end = max(growth), Nfold_at_end = exp(max(growth)),
    late_slope = lt, late_mean_sd = late_mean_sd,
    curvature = curv,
    p_sat_growth = p_sat_growth, rss_gain_growth = rss_gain_growth,
    p_sat_o2 = p_sat_o2, rss_gain_o2 = rss_gain_o2)

  prof[[length(prof) + 1]] <- tibble::tibble(
    key, t = t, resid = e, resid_sd = e / sig, o2_abs = o2_abs, growth = growth,
    frac_window = (t - min(t)) / (max(t) - min(t)))
}

A1 <- dplyr::bind_rows(rows)
d6_write(A1, "A1_residual_structure_per_curve.csv")

PR <- dplyr::bind_rows(prof)
# Binned residual profile against the two candidate axes.
A2 <- PR |>
  dplyr::mutate(bin = cut(frac_window, breaks = seq(0, 1, by = 0.1),
                          include.lowest = TRUE, labels = FALSE)) |>
  dplyr::group_by(bin) |>
  dplyr::summarise(frac_window_mid = (bin[1] - 0.5) / 10,
                   n = dplyr::n(), mean_resid_sd = mean(resid_sd),
                   median_resid_sd = stats::median(resid_sd),
                   .groups = "drop")
d6_write(A2, "A2_residual_profile_by_window_position.csv")

A3 <- PR |>
  dplyr::mutate(bin = cut(o2_abs, breaks = stats::quantile(o2_abs, seq(0, 1, 0.1)),
                          include.lowest = TRUE, labels = FALSE)) |>
  dplyr::group_by(bin) |>
  dplyr::summarise(o2_abs_mid = stats::median(o2_abs), n = dplyr::n(),
                   mean_resid_sd = mean(resid_sd), .groups = "drop")
d6_write(A3, "A3_residual_profile_by_absolute_O2.csv")

# ---- how many curves show a detectable deviation? ---------------------------
A4 <- tibble::tibble(
  test = c("nested F-test vs growth saturation, p < 0.05",
           "nested F-test vs growth saturation, p < 0.01",
           "nested F-test vs O2 limitation, p < 0.05",
           "nested F-test vs O2 limitation, p < 0.01",
           "|late-window mean residual| > 0.5 sd",
           "|late-window mean residual| > 1.0 sd",
           "window ends below 20% of starting O2"),
  n_curves = c(sum(A1$p_sat_growth < 0.05, na.rm = TRUE),
               sum(A1$p_sat_growth < 0.01, na.rm = TRUE),
               sum(A1$p_sat_o2 < 0.05, na.rm = TRUE),
               sum(A1$p_sat_o2 < 0.01, na.rm = TRUE),
               sum(abs(A1$late_mean_sd) > 0.5, na.rm = TRUE),
               sum(abs(A1$late_mean_sd) > 1.0, na.rm = TRUE),
               sum(A1$O2_frac_remaining < 0.20, na.rm = TRUE)),
  of_total = nrow(A1))
d6_write(A4, "A4_deviation_counts.csv")

# ---- mechanism discrimination -----------------------------------------------
# Which better explains the late-window residual signature: how LOW the oxygen
# got, or how MUCH GROWTH elapsed?
cc <- function(x, y) suppressWarnings(stats::cor(x, y, use = "complete.obs"))
A5 <- tibble::tibble(
  predictor = c("min absolute O2 reached (mg/L)",
                "fraction of O2 remaining at window end",
                "elapsed growth r*t at window end",
                "N-fold increase at window end"),
  mechanism = c("O2 limitation", "O2 limitation",
                "growth saturation", "growth saturation"),
  corr_with_late_mean_resid = c(
    cc(A1$min_O2_abs, A1$late_mean_sd),
    cc(A1$O2_frac_remaining, A1$late_mean_sd),
    cc(A1$growth_at_end, A1$late_mean_sd),
    cc(A1$Nfold_at_end, A1$late_mean_sd)),
  corr_with_curvature = c(
    cc(A1$min_O2_abs, A1$curvature),
    cc(A1$O2_frac_remaining, A1$curvature),
    cc(A1$growth_at_end, A1$curvature),
    cc(A1$Nfold_at_end, A1$curvature)))
d6_write(A5, "A5_mechanism_discrimination.csv")

message(sprintf("\n  nested F-test vs growth saturation: p<0.05 in %d/%d curves; vs O2 limitation: %d/%d",
  sum(A1$p_sat_growth < 0.05, na.rm = TRUE), nrow(A1),
  sum(A1$p_sat_o2 < 0.05, na.rm = TRUE), nrow(A1)))
message(sprintf("  median RSS reduction from the extra term: growth %.3f%%, O2 %.3f%%",
  100*stats::median(A1$rss_gain_growth, na.rm = TRUE),
  100*stats::median(A1$rss_gain_o2, na.rm = TRUE)))
message(sprintf("  late-window mean residual (in sd): median %.3f, range %.3f to %.3f",
  stats::median(A1$late_mean_sd, na.rm = TRUE), min(A1$late_mean_sd, na.rm = TRUE),
  max(A1$late_mean_sd, na.rm = TRUE)))
message(sprintf("  min O2 remaining fraction: median %.3f, range %.3f-%.3f",
  stats::median(A1$O2_frac_remaining), min(A1$O2_frac_remaining), max(A1$O2_frac_remaining)))
message(sprintf("  N-fold growth at window end: median %.2f, range %.2f-%.2f",
  stats::median(A1$Nfold_at_end), min(A1$Nfold_at_end), max(A1$Nfold_at_end)))
print(as.data.frame(A5), row.names = FALSE, digits = 3)
message("D6 PART A done.")
