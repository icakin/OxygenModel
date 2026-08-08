# =============================================================================
# D6 PART D - the key test: window sweep under each model, plus a power check
# =============================================================================
# Three things:
#
#  (1) POWER CHECK. PART C found no support for lag or saturation. That is only
#      meaningful if the same models COULD detect them when present. This
#      re-uses D5's misspecification generators (scripts/10_simulation_recovery.R)
#      rather than writing new ones, so simulated and observed diagnostics sit on
#      the same footing, and asks: on a curve simulated WITH a 15-min lag or with
#      saturation, does M1/M2g recover it and does AIC prefer it?
#
#  (2) THE WINDOW SWEEP under M0 and under the best alternative. If an extension
#      were capturing something real, window sensitivity should FALL. This is the
#      strongest available test against overfitting.
#
#  (3) Between-taxon ordering of per-cell R, and the Fig 6 growth-respiration
#      slope, under each model.
#
# OUTPUT  runs/D6_analysis/D1_*.csv .. D5_*.csv
# =============================================================================

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d6_common.R"))

set.seed(20260807)
message("\n== D6 PART D: power check, window sweep, downstream ==")

curves <- readr::read_csv(FITCURVES_CSV,     show_col_types = FALSE)
fits   <- readr::read_csv(RESULTS_FIT_CSV,   show_col_types = FALSE)
res    <- readr::read_csv(RESULTS_FINAL_CSV, show_col_types = FALSE)
C1     <- readr::read_csv(d6f("C1_model_fits_per_curve.csv"), show_col_types = FALSE)
ID     <- id_cols(fits)

# ---------------------------------------------------------------------------
# (1) POWER CHECK - reusing D5's generators from 10_simulation_recovery.R
# ---------------------------------------------------------------------------
# D5's arms: "lag" shifts the onset by lag_min; "saturating" decays the growth
# rate logistically over the window. Reproduced here with the same shapes and
# the same trimming rule (stop at 40% of start), then fitted with M0/M1/M2g.
N0_true <- 1e9; O2_ref_true <- 1; O2_0_true <- 1

d6_sim <- function(arm, r_true, R_true, noise_sd, dt, duration,
                   lag_min = 15, sat_frac = 0.5) {
  K_true <- (R_true * N0_true) / O2_ref_true
  t_full <- seq(0, duration, by = dt)
  O2 <- switch(arm,
    "none" = m0_model(r_true, K_true, t_full, O2_0_true),
    "lag"  = m0_model(r_true, K_true, pmax(t_full - lag_min, 0), O2_0_true),
    "saturating" = {
      rt  <- r_true / (1 + exp((t_full - sat_frac * duration) / (0.1 * duration)))
      cum <- cumsum(c(0, diff(t_full)) * rt)
      O2_0_true - (K_true / r_true) * (exp(cum) - 1)
    })
  idx <- which(O2 <= 0.4 * O2_0_true)[1]
  if (!is.na(idx) && idx > 5) { t <- t_full[1:idx]; O2 <- O2[1:idx] } else t <- t_full
  y <- O2 + stats::rnorm(length(t), 0, noise_sd); y[y <= 0] <- 1e-6
  list(t = t, y = y, K_true = K_true)
}

pow_rows <- list()
pgrid <- tidyr::expand_grid(arm = c("none", "lag", "saturating"),
                            r_true = c(0.01, 0.02, 0.035),
                            R_true = c(2e-12, 5e-12, 8e-12),
                            rep_id = 1:12)
for (j in seq_len(nrow(pgrid))) {
  s <- d6_sim(pgrid$arm[j], pgrid$r_true[j], pgrid$R_true[j], 0.005, 1, 120)
  d <- data.frame(t = s$t, y = s$y)
  s0 <- list(r = pgrid$r_true[j], K = s$K_true, O2_0 = 1)
  tlag_max <- 0.5 * (max(s$t) - min(s$t))
  f0 <- fit_m0(s$t, s$y, s0)
  f1 <- fit_try(y ~ m1_model(r, K, t, O2_0, tlag), d, c(s0, tlag = 5),
                c(FIT_LOWER, tlag = 0), c(FIT_UPPER, tlag = tlag_max))
  f2 <- fit_try(y ~ m2g_model(r, K, t, O2_0, Nmax_rel), d, c(s0, Nmax_rel = 50),
                c(FIT_LOWER, Nmax_rel = 1.01), c(FIT_UPPER, Nmax_rel = 1e6))
  gt <- function(f) if (inherits(f, "try-error") || is.null(f)) NA_real_ else unname(stats::AIC(f))
  pow_rows[[j]] <- tibble::tibble(
    arm = pgrid$arm[j], r_true = pgrid$r_true[j], R_true = pgrid$R_true[j],
    AIC_M0 = gt(f0), AIC_M1 = gt(f1), AIC_M2g = gt(f2),
    tlag_hat = if (inherits(f1, "try-error")) NA_real_ else unname(stats::coef(f1)[["tlag"]]),
    Nmax_hat = if (inherits(f2, "try-error")) NA_real_ else unname(stats::coef(f2)[["Nmax_rel"]]))
}
POW <- dplyr::bind_rows(pow_rows) |>
  dplyr::mutate(dAIC_M1 = AIC_M1 - AIC_M0, dAIC_M2g = AIC_M2g - AIC_M0)

D1 <- POW |> dplyr::group_by(arm) |>
  dplyr::summarise(
    n = dplyr::n(),
    pct_M1_preferred = 100 * mean(dAIC_M1 < -2, na.rm = TRUE),
    median_tlag_hat = stats::median(tlag_hat, na.rm = TRUE),
    pct_M2g_preferred = 100 * mean(dAIC_M2g < -2, na.rm = TRUE),
    median_Nmax_hat = stats::median(Nmax_hat, na.rm = TRUE),
    .groups = "drop")
d6_write(D1, "D1_power_check.csv")

# ---------------------------------------------------------------------------
# (2) THE WINDOW SWEEP under M0 and the best alternative
# ---------------------------------------------------------------------------
# Same design as 08: both edges shifted over {-6,-3,0,+3,+6} minutes. 08 itself
# is not modified and its committed output is untouched; this re-runs the sweep
# here so both models are swept identically.
OFF <- c(-6, -3, 0, 3, 6)
dep <- load_depletion_table(); fcf <- load_fc_final()

sweep_rows <- list()
for (i in seq_len(nrow(fits))) {
  key <- fits[i, ID, drop = FALSE]
  gall <- curves |> dplyr::semi_join(key, by = ID) |> dplyr::arrange(Time0_min)
  rr <- res |> dplyr::semi_join(key, by = ID)
  if (!nrow(rr) || nrow(gall) < 20) next
  dr <- dep |> dplyr::filter(Taxon == rr$Taxon[1], Replicate == rr$Replicate[1])
  fr <- fcf |> dplyr::filter(Taxon == rr$Taxon[1], Replicate == rr$Replicate[1])
  if (!nrow(dr) || !nrow(fr)) next
  t_dep <- dr$t_depletion_min[1]; FC <- fr$FC_Final[1]; O2ref <- rr$O2_ref[1]
  t0abs <- rr$fit_start_min[1]
  s0 <- list(r = rr$r_per_minute[1], K = rr$K[1], O2_0 = rr$O2_0[1])
  tmin <- min(gall$Time0_min); tmax <- max(gall$Time0_min)

  for (ds in OFF) for (de in OFF) {
    sel <- gall$Time0_min >= (tmin + ds) & gall$Time0_min <= (tmax + de)
    g <- gall[sel, ]
    if (nrow(g) < 15) next
    t <- g$Time0_min - min(g$Time0_min); y <- g$Oxygen_norm
    d <- data.frame(t = t, y = y)
    start_abs <- t0abs + ds; T_end <- max(t)
    for (nm in c("M0", "M2g")) {
      f <- if (nm == "M0") fit_m0(t, y, s0) else
        fit_try(y ~ m2g_model(r, K, t, O2_0, Nmax_rel), d, c(s0, Nmax_rel = 50),
                c(FIT_LOWER, Nmax_rel = 1.01), c(FIT_UPPER, Nmax_rel = 1e6))
      if (inherits(f, "try-error") || is.null(f)) next
      cf <- stats::coef(f); pred <- y - stats::residuals(f)
      N0s <- FC * FC_TO_CELLS_PER_L * exp(-cf[["r"]] * (t_dep - start_abs))
      Rp <- percell_R(pred[1] - pred[length(pred)], O2ref, N0s, cf[["r"]], T_end)
      sweep_rows[[length(sweep_rows) + 1]] <- tibble::tibble(
        key, model = nm, ds = ds, de = de, is_base = (ds == 0 & de == 0),
        r = cf[["r"]], K = cf[["K"]], Rpc = Rp)
    }
  }
  if (i %% 15 == 0) message("    swept ", i, "/", nrow(fits))
}
SW <- dplyr::bind_rows(sweep_rows)

bs <- SW |> dplyr::filter(is_base) |>
  dplyr::select(dplyr::all_of(ID), model, r0 = r, K0 = K, R0 = Rpc)
SW <- SW |> dplyr::left_join(bs, by = c(ID, "model")) |>
  dplyr::mutate(dK_pct = 100 * (K - K0) / K0,
                dR_pct = 100 * (Rpc - R0) / R0,
                dr_pct = 100 * (r - r0) / r0)
d6_write(SW, "D2_window_sweep_raw.csv")

D3 <- SW |> dplyr::filter(!is_base) |>
  dplyr::group_by(dplyr::across(dplyr::all_of(ID)), model) |>
  dplyr::summarise(mK = max(abs(dK_pct), na.rm = TRUE),
                   mR = max(abs(dR_pct), na.rm = TRUE),
                   mr = max(abs(dr_pct), na.rm = TRUE), .groups = "drop") |>
  dplyr::group_by(model) |>
  dplyr::summarise(n_curves = dplyr::n(),
                   median_max_dK_pct = stats::median(mK, na.rm = TRUE),
                   p90_max_dK_pct = stats::quantile(mK, .9, na.rm = TRUE),
                   median_max_dR_pct = stats::median(mR, na.rm = TRUE),
                   p90_max_dR_pct = stats::quantile(mR, .9, na.rm = TRUE),
                   median_max_dr_pct = stats::median(mr, na.rm = TRUE),
                   .groups = "drop")
d6_write(D3, "D3_window_sweep_by_model.csv")

# ---------------------------------------------------------------------------
# (3) Ordering and the Fig 6 slope under each model
# ---------------------------------------------------------------------------
tx <- C1 |> dplyr::filter(converged, is.finite(R_C_fg_cell_h)) |>
  dplyr::group_by(model, Taxon) |>
  dplyr::summarise(R_taxon = mean(R_C_fg_cell_h), r_taxon = mean(r), .groups = "drop")
ref <- tx |> dplyr::filter(model == "M0") |> dplyr::select(Taxon, R0 = R_taxon)
D4 <- tx |> dplyr::left_join(ref, by = "Taxon") |>
  dplyr::group_by(model) |>
  dplyr::summarise(
    spearman_vs_M0 = stats::cor(R_taxon, R0, method = "spearman"),
    pearson_vs_M0 = stats::cor(R_taxon, R0),
    median_pct_shift_R = 100 * stats::median((R_taxon - R0) / R0),
    .groups = "drop")

# The growth-respiration relation. NOTE: this is the POOLED log-log slope across
# all 75 curves, NOT the published Fig 6 slope, which 06_main_figures.R obtains
# from a random-intercept mixed model on (1 | Taxon) and reports as 0.283. What
# is compared here is the same quantity under each model, so the MODEL EFFECT is
# like-for-like; the absolute value is not the published figure.
slope <- C1 |> dplyr::filter(converged, is.finite(R_C_fg_cell_h), R_C_fg_cell_h > 0) |>
  dplyr::group_by(model) |>
  dplyr::summarise(
    pooled_logR_logr_slope = unname(stats::coef(stats::lm(log(R_C_fg_cell_h) ~ log(r)))[2]),
    pooled_r2 = summary(stats::lm(log(R_C_fg_cell_h) ~ log(r)))$r.squared,
    .groups = "drop")
D4 <- D4 |> dplyr::left_join(slope, by = "model")
d6_write(D4, "D4_ordering_and_fig6.csv")

D5 <- tibble::tibble(
  quantity = c("median max|dK_pct| under M0", "median max|dK_pct| under M2g",
               "median max|dR_pct| under M0", "median max|dR_pct| under M2g",
               "Spearman ordering M2g vs M0", "pooled log R ~ log r slope under M0 (NOT the published Fig 6 slope)",
               "pooled log R ~ log r slope under M2g"),
  value = c(D3$median_max_dK_pct[D3$model == "M0"], D3$median_max_dK_pct[D3$model == "M2g"],
            D3$median_max_dR_pct[D3$model == "M0"], D3$median_max_dR_pct[D3$model == "M2g"],
            D4$spearman_vs_M0[D4$model == "M2g"], D4$pooled_logR_logr_slope[D4$model == "M0"],
            D4$pooled_logR_logr_slope[D4$model == "M2g"]))
d6_write(D5, "D5_headline.csv")

print(as.data.frame(D1), row.names = FALSE, digits = 4); cat("\n")
print(as.data.frame(D3), row.names = FALSE, digits = 4); cat("\n")
print(as.data.frame(D4), row.names = FALSE, digits = 4)
message("D6 PART D done.")
