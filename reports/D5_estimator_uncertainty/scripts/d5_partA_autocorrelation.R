# =============================================================================
# D5 PART A - residual autocorrelation, and what it does to vcov(nls)
# =============================================================================
# Refits every curve exactly as 05_oxygen_fits.R does, keeps the residual series
# in fit-window time order, and measures:
#
#   * lag-1 residual autocorrelation, Durbin-Watson, an AR(1) rho
#   * the implied effective sample size n_eff = n (1-rho)/(1+rho)
#   * the variance-inflation factor sqrt((1+rho)/(1-rho))
#   * the naive vcov(nls) standard errors, and two autocorrelation-corrected
#     alternatives: an AR(1) GLS refit (headline) and Newey-West (cross-check)
#
# INPUT   results/tables/oxygen_fit_curves.csv, oxygen_fit_results.csv
# OUTPUT  runs/D5_analysis/A1_*.csv .. A4_*.csv
# RUN     Rscript reports/D5_estimator_uncertainty/scripts/d5_partA_autocorrelation.R
#
# DIAGNOSTIC ONLY - writes nothing outside runs/D5_analysis/.
# =============================================================================

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d5_common.R"))

message("\n== D5 PART A: residual autocorrelation ==")

curves <- readr::read_csv(FITCURVES_CSV,   show_col_types = FALSE)
fits   <- readr::read_csv(RESULTS_FIT_CSV, show_col_types = FALSE)
ID     <- id_cols(fits)
message("  curve key: ", paste(ID, collapse = " x "))

# Newey-West bandwidth. Two are reported: the textbook rule, and one matched to
# the AR(1) correlation length (the lag at which the AR(1) ACF falls below 0.05,
# i.e. log(0.05)/log(rho)), capped at n/4 so the estimator stays sane.
nw_bandwidth_matched <- function(n, rho) {
  if (!is.finite(rho) || rho <= 0) return(1L)
  L <- ceiling(log(0.05) / log(rho))
  as.integer(max(1L, min(L, floor(n / 4))))
}

rows <- vector("list", nrow(fits))
acf_rows <- vector("list", nrow(fits))   # ACF at lags 1..MAXLAG, per curve
MAXLAG <- 20L

for (i in seq_len(nrow(fits))) {
  key <- fits[i, ID, drop = FALSE]
  g <- curves |> dplyr::semi_join(key, by = ID) |> dplyr::arrange(Time0_min)
  if (nrow(g) < 10) next

  t <- g$Time0_min; y <- g$Oxygen_norm; n <- length(y)
  start <- list(r = fits$r_per_minute[i], K = fits$K[i], O2_0 = fits$O2_0[i])

  f <- fit_ols(t, y, start)
  if (inherits(f, "try-error")) next
  cf <- stats::coef(f)
  e  <- y - stats::predict(f)

  # ---- naive, exactly what the pipeline sees -------------------------------
  se_naive <- tryCatch(summary(f)$coefficients[, "Std. Error"],
                       error = function(err) rep(NA_real_, 3))
  V_naive  <- tryCatch(stats::vcov(f), error = function(err) NULL)

  # ---- autocorrelation of the residuals ------------------------------------
  rho  <- ar1_rho(e)
  acf1 <- suppressWarnings(stats::acf(e, lag.max = 1, plot = FALSE)$acf[2])
  dw   <- durbin_watson(e)
  # Formal test: Ljung-Box on the first 10 lags. The null is "no residual
  # autocorrelation"; a small p-value is evidence the residuals are correlated.
  lb   <- suppressWarnings(stats::Box.test(e, lag = 10, type = "Ljung-Box"))
  # Full ACF profile, so the report can show there is no structure at ANY lag,
  # not merely at lag 1.
  ac   <- suppressWarnings(stats::acf(e, lag.max = MAXLAG, plot = FALSE)$acf[-1])
  acf_rows[[i]] <- tibble::tibble(key, lag = seq_along(ac), acf = as.numeric(ac))
  neff <- n_eff_ar1(n, rho)
  vif  <- vif_ar1(rho)

  # ---- Newey-West sandwich --------------------------------------------------
  J <- model_jacobian(cf["r"], cf["K"], t, cf["O2_0"])
  L_rule    <- max(1L, floor(4 * (n / 100)^(2 / 9)))
  L_matched <- nw_bandwidth_matched(n, rho)
  V_nw_rule <- nw_vcov(J, e, L = L_rule)
  V_nw_mtch <- nw_vcov(J, e, L = L_matched)

  # ---- AR(1) GLS refit ------------------------------------------------------
  gls <- fit_ar1_gls(t, y, start)
  se_gls <- rep(NA_real_, 3); names(se_gls) <- c("r", "K", "O2_0")
  rho_gls <- NA_real_; cf_gls <- rep(NA_real_, 3); names(cf_gls) <- names(se_gls)
  if (!is.null(gls)) {
    se_gls  <- tryCatch(summary(gls$fit)$coefficients[, "Std. Error"],
                        error = function(err) se_gls)
    cf_gls  <- stats::coef(gls$fit)
    rho_gls <- gls$rho
  }

  gv <- function(V, p) if (is.null(V)) NA_real_ else sqrt(max(V[p, p], 0))

  rows[[i]] <- tibble::tibble(
    key,
    n = n, span_min = max(t) - min(t), dt_median = stats::median(diff(t)),
    r = unname(cf["r"]), K = unname(cf["K"]), O2_0 = unname(cf["O2_0"]),
    sigma = sqrt(sum(e^2) / (n - 3)),
    acf_lag1 = unname(acf1), rho_ar1 = rho, durbin_watson = dw,
    ljung_box_p = unname(lb$p.value), ljung_box_stat = unname(lb$statistic),
    n_eff = neff, vif_analytic = vif,
    L_nw_rule = L_rule, L_nw_matched = L_matched,
    se_r_naive = unname(se_naive["r"]), se_K_naive = unname(se_naive["K"]),
    se_r_nw_rule = gv(V_nw_rule, "r"), se_K_nw_rule = gv(V_nw_rule, "K"),
    se_r_nw_matched = gv(V_nw_mtch, "r"), se_K_nw_matched = gv(V_nw_mtch, "K"),
    se_r_gls = unname(se_gls["r"]), se_K_gls = unname(se_gls["K"]),
    rho_gls = rho_gls,
    r_gls = unname(cf_gls["r"]), K_gls = unname(cf_gls["K"]),
    vcov_corr_rK = if (is.null(V_naive)) NA_real_ else
      V_naive["r", "K"] / sqrt(V_naive["r", "r"] * V_naive["K", "K"])
  )
  if (i %% 15 == 0) message("    ", i, "/", nrow(fits), " curves")
}

A1 <- dplyr::bind_rows(rows) |>
  dplyr::mutate(
    ratio_r_gls        = se_r_gls / se_r_naive,
    ratio_K_gls        = se_K_gls / se_K_naive,
    ratio_r_nw_matched = se_r_nw_matched / se_r_naive,
    ratio_K_nw_matched = se_K_nw_matched / se_K_naive,
    ratio_r_nw_rule    = se_r_nw_rule / se_r_naive,
    rel_se_r_naive     = se_r_naive / r,
    rel_se_r_gls       = se_r_gls / r,
    rel_se_K_naive     = se_K_naive / K,
    rel_se_K_gls       = se_K_gls / K
  )
d5_write(A1, "A1_per_curve_residual_diagnostics.csv")

# ---- distributions ----------------------------------------------------------
q <- function(x) {
  x <- x[is.finite(x)]
  tibble::tibble(n = length(x), min = min(x), q25 = stats::quantile(x, .25),
                 median = stats::median(x), q75 = stats::quantile(x, .75),
                 max = max(x))
}
A2 <- dplyr::bind_rows(
  dplyr::mutate(q(A1$rho_ar1),       quantity = "rho_ar1",       .before = 1),
  dplyr::mutate(q(A1$acf_lag1),      quantity = "acf_lag1",      .before = 1),
  dplyr::mutate(q(A1$durbin_watson), quantity = "durbin_watson", .before = 1),
  dplyr::mutate(q(A1$n),             quantity = "n_points",      .before = 1),
  dplyr::mutate(q(A1$n_eff),         quantity = "n_eff",         .before = 1),
  dplyr::mutate(q(A1$vif_analytic),  quantity = "vif_analytic",  .before = 1)
)
d5_write(A2, "A2_autocorrelation_summary_overall.csv")

A2t <- A1 |> dplyr::group_by(Taxon) |>
  dplyr::summarise(
    n_curves = dplyr::n(),
    median_n = stats::median(n),
    median_rho = stats::median(rho_ar1),
    min_rho = min(rho_ar1), max_rho = max(rho_ar1),
    median_dw = stats::median(durbin_watson),
    median_n_eff = stats::median(n_eff),
    median_vif = stats::median(vif_analytic),
    .groups = "drop") |>
  dplyr::arrange(dplyr::desc(median_rho))
d5_write(A2t, "A2_autocorrelation_summary_by_taxon.csv")

A3 <- dplyr::bind_rows(
  dplyr::mutate(q(A1$ratio_r_gls),        quantity = "se_r  GLS / naive",        .before = 1),
  dplyr::mutate(q(A1$ratio_K_gls),        quantity = "se_K  GLS / naive",        .before = 1),
  dplyr::mutate(q(A1$ratio_r_nw_matched), quantity = "se_r  NW(matched) / naive", .before = 1),
  dplyr::mutate(q(A1$ratio_K_nw_matched), quantity = "se_K  NW(matched) / naive", .before = 1),
  dplyr::mutate(q(A1$ratio_r_nw_rule),    quantity = "se_r  NW(rule) / naive",    .before = 1)
)
d5_write(A3, "A3_se_inflation_summary.csv")

# ACF profile across all curves, lag by lag. Under white noise each lag should
# sit inside +/- 1.96/sqrt(n); the band is reported alongside so the report can
# show there is no structure at any lag, not merely at lag 1.
A5 <- dplyr::bind_rows(acf_rows) |>
  dplyr::group_by(lag) |>
  dplyr::summarise(
    n_curves = dplyr::n(),
    mean_acf = mean(acf), median_acf = stats::median(acf),
    q05 = stats::quantile(acf, .05), q95 = stats::quantile(acf, .95),
    frac_outside_white_noise_band =
      mean(abs(acf) > 1.96 / sqrt(stats::median(A1$n))),
    .groups = "drop") |>
  dplyr::mutate(white_noise_band = 1.96 / sqrt(stats::median(A1$n)))
d5_write(A5, "A5_residual_acf_by_lag.csv")

A6 <- tibble::tibble(
  quantity = c("curves with Ljung-Box p < 0.05 (lag 10)",
               "curves with Ljung-Box p < 0.01 (lag 10)",
               "curves with |rho| > 0.5",
               "curves with rho > 0.9 (the hypothesised regime)",
               "curves where AR(1) GLS SE exceeds naive by >20%"),
  n_curves = c(sum(A1$ljung_box_p < 0.05, na.rm = TRUE),
               sum(A1$ljung_box_p < 0.01, na.rm = TRUE),
               sum(abs(A1$rho_ar1) > 0.5, na.rm = TRUE),
               sum(A1$rho_ar1 > 0.9, na.rm = TRUE),
               sum(A1$ratio_r_gls > 1.2, na.rm = TRUE)),
  of_total = nrow(A1))
d5_write(A6, "A6_autocorrelation_hypothesis_test.csv")

A4 <- A1 |> dplyr::group_by(Taxon) |>
  dplyr::summarise(
    median_ratio_r_gls = stats::median(ratio_r_gls, na.rm = TRUE),
    median_ratio_K_gls = stats::median(ratio_K_gls, na.rm = TRUE),
    median_rel_se_r_naive = stats::median(rel_se_r_naive, na.rm = TRUE),
    median_rel_se_r_gls   = stats::median(rel_se_r_gls, na.rm = TRUE),
    .groups = "drop")
d5_write(A4, "A4_se_inflation_by_taxon.csv")

message(sprintf(
  "\n  rho: median %.4f (range %.4f-%.4f) | n: median %d -> n_eff median %.1f",
  stats::median(A1$rho_ar1), min(A1$rho_ar1), max(A1$rho_ar1),
  as.integer(stats::median(A1$n)), stats::median(A1$n_eff)))
message(sprintf(
  "  analytic VIF sqrt((1+rho)/(1-rho)): median %.1fx (range %.1f-%.1f)",
  stats::median(A1$vif_analytic), min(A1$vif_analytic), max(A1$vif_analytic)))
message(sprintf(
  "  AR(1) GLS SE / naive SE: r median %.2fx (%.2f-%.2f), K median %.2fx (%.2f-%.2f)",
  stats::median(A1$ratio_r_gls, na.rm = TRUE), min(A1$ratio_r_gls, na.rm = TRUE),
  max(A1$ratio_r_gls, na.rm = TRUE),
  stats::median(A1$ratio_K_gls, na.rm = TRUE), min(A1$ratio_K_gls, na.rm = TRUE),
  max(A1$ratio_K_gls, na.rm = TRUE)))
message("D5 PART A done.")
