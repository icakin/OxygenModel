# =============================================================================
# D6 PART C - lag and saturation model variants, as ALTERNATIVE fits
# =============================================================================
# The default estimator is NOT changed. Every model here is an additional fit
# writing to runs/D6_analysis/.
#
#   M0   the pipeline's model, as reference
#   M1   + explicit lag/equilibration (nuisance, absorbing the instrumental
#        stabilisation transient the paper's methods describe - not biology)
#   M2o  + O2 limitation (Michaelis-Menten in O2, Km)
#   M2g  + growth saturation (logistic biomass, Nmax_rel)
#   M3o  M1 + M2o        M3g  M1 + M2g
#
# PART A could not separate the two mechanisms because neither is present, so
# BOTH are implemented and compared, as the brief requires in that case.
#
# nlsLM throughout. Escalation to a Bayesian fit would be for IDENTIFICATION
# only - see the identification diagnostics written here - never for uncertainty,
# which D5 established is not the issue.
#
# OUTPUT  runs/D6_analysis/C1_*.csv .. C4_*.csv
# =============================================================================

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d6_common.R"))

message("\n== D6 PART C: model variants M0-M3 ==")

curves <- readr::read_csv(FITCURVES_CSV,     show_col_types = FALSE)
fits   <- readr::read_csv(RESULTS_FIT_CSV,   show_col_types = FALSE)
res    <- readr::read_csv(RESULTS_FINAL_CSV, show_col_types = FALSE)
ID     <- id_cols(fits)

# Consumed oxygen over the window, from a model's own trajectory: needed so
# per-cell R is computed on the same footing for every variant.
consumed <- function(pred) pred[1] - pred[length(pred)]

rows <- list()
for (i in seq_len(nrow(fits))) {
  key <- fits[i, ID, drop = FALSE]
  g <- curves |> dplyr::semi_join(key, by = ID) |> dplyr::arrange(Time0_min)
  if (nrow(g) < 15) next
  t <- g$Time0_min; y <- g$Oxygen_norm; n <- length(t)
  d <- data.frame(t = t, y = y)
  rr <- res |> dplyr::semi_join(key, by = ID)
  if (!nrow(rr)) next
  O2ref <- rr$O2_ref[1]; T_end <- rr$T_end_min[1]
  N0 <- rr$N0_cells_per_L[1]
  s0 <- list(r = fits$r_per_minute[i], K = fits$K[i], O2_0 = fits$O2_0[i])

  # window length caps the lag: a lag longer than the window is meaningless
  tlag_max <- 0.5 * (max(t) - min(t))

  M <- list()
  M$M0  <- list(fit = fit_m0(t, y, s0), k = 3,
                lower = FIT_LOWER, upper = FIT_UPPER)
  M$M1  <- list(fit = fit_try(y ~ m1_model(r, K, t, O2_0, tlag), d,
                  c(s0, tlag = min(5, tlag_max / 2)),
                  c(FIT_LOWER, tlag = 0), c(FIT_UPPER, tlag = tlag_max)),
                k = 4, lower = c(FIT_LOWER, tlag = 0),
                upper = c(FIT_UPPER, tlag = tlag_max))
  M$M2o <- list(fit = fit_try(y ~ m2o_predict(r, K, t, O2_0, Km), d,
                  c(s0, Km = 0.05), c(FIT_LOWER, Km = 1e-6), c(FIT_UPPER, Km = 10)),
                k = 4, lower = c(FIT_LOWER, Km = 1e-6), upper = c(FIT_UPPER, Km = 10))
  M$M2g <- list(fit = fit_try(y ~ m2g_model(r, K, t, O2_0, Nmax_rel), d,
                  c(s0, Nmax_rel = 50), c(FIT_LOWER, Nmax_rel = 1.01),
                  c(FIT_UPPER, Nmax_rel = 1e6)),
                k = 4, lower = c(FIT_LOWER, Nmax_rel = 1.01),
                upper = c(FIT_UPPER, Nmax_rel = 1e6))
  M$M3o <- list(fit = fit_try(y ~ m3o_predict(r, K, t, O2_0, Km, tlag), d,
                  c(s0, Km = 0.05, tlag = min(5, tlag_max / 2)),
                  c(FIT_LOWER, Km = 1e-6, tlag = 0),
                  c(FIT_UPPER, Km = 10, tlag = tlag_max)),
                k = 5, lower = c(FIT_LOWER, Km = 1e-6, tlag = 0),
                upper = c(FIT_UPPER, Km = 10, tlag = tlag_max))
  M$M3g <- list(fit = fit_try(y ~ m3g_model(r, K, t, O2_0, Nmax_rel, tlag), d,
                  c(s0, Nmax_rel = 50, tlag = min(5, tlag_max / 2)),
                  c(FIT_LOWER, Nmax_rel = 1.01, tlag = 0),
                  c(FIT_UPPER, Nmax_rel = 1e6, tlag = tlag_max)),
                k = 5, lower = c(FIT_LOWER, Nmax_rel = 1.01, tlag = 0),
                upper = c(FIT_UPPER, Nmax_rel = 1e6, tlag = tlag_max))

  for (nm in names(M)) {
    m <- M[[nm]]
    if (inherits(m$fit, "try-error") || is.null(m$fit)) {
      rows[[length(rows) + 1]] <- tibble::tibble(
        key, model = nm, converged = FALSE, k = m$k,
        AIC = NA_real_, BIC = NA_real_, rss = NA_real_, sigma = NA_real_,
        r = NA_real_, K = NA_real_, O2_0 = NA_real_,
        extra1 = NA_real_, extra2 = NA_real_, at_bounds = NA_character_,
        R_C_fg_cell_h = NA_real_, late_mean_sd = NA_real_)
      next
    }
    cf <- stats::coef(m$fit); e <- stats::residuals(m$fit)
    icv <- ic(m$fit, n)
    pred <- y - e
    Rp <- percell_R(consumed(pred), O2ref, N0, cf[["r"]], T_end)
    cut <- t >= (min(t) + (2/3) * (max(t) - min(t)))
    ex <- setdiff(names(cf), c("r", "K", "O2_0"))
    rows[[length(rows) + 1]] <- tibble::tibble(
      key, model = nm, converged = TRUE, k = m$k,
      AIC = icv[["AIC"]], BIC = icv[["BIC"]], rss = icv[["rss"]],
      sigma = sqrt(icv[["rss"]] / (n - m$k)),
      r = cf[["r"]], K = cf[["K"]], O2_0 = cf[["O2_0"]],
      extra1 = if (length(ex) >= 1) unname(cf[[ex[1]]]) else NA_real_,
      extra2 = if (length(ex) >= 2) unname(cf[[ex[2]]]) else NA_real_,
      at_bounds = at_bounds(m$fit, m$lower, m$upper),
      R_C_fg_cell_h = Rp,
      late_mean_sd = mean(e[cut]) / stats::sd(e))
  }
  if (i %% 15 == 0) message("    ", i, "/", nrow(fits), " curves")
}

C1 <- dplyr::bind_rows(rows)
d6_write(C1, "C1_model_fits_per_curve.csv")

# ---- convergence, identification, information criteria ----------------------
C2 <- C1 |>
  dplyr::group_by(model) |>
  dplyr::summarise(
    n_curves = dplyr::n(),
    n_converged = sum(converged),
    n_at_bounds = sum(nzchar(at_bounds) & !is.na(at_bounds)),
    pct_at_bounds = 100 * sum(nzchar(at_bounds) & !is.na(at_bounds)) / dplyr::n(),
    median_AIC = stats::median(AIC, na.rm = TRUE),
    median_BIC = stats::median(BIC, na.rm = TRUE),
    median_sigma = stats::median(sigma, na.rm = TRUE),
    median_late_resid_sd = stats::median(late_mean_sd, na.rm = TRUE),
    .groups = "drop")
d6_write(C2, "C2_model_summary.csv")

# ---- head-to-head against M0 -------------------------------------------------
base <- C1 |> dplyr::filter(model == "M0") |>
  dplyr::select(dplyr::all_of(ID), AIC0 = AIC, BIC0 = BIC, rss0 = rss,
                r0 = r, K0 = K, R0 = R_C_fg_cell_h)
C3 <- C1 |> dplyr::left_join(base, by = ID) |>
  dplyr::filter(model != "M0") |>
  dplyr::mutate(dAIC = AIC - AIC0, dBIC = BIC - BIC0,
                rss_reduction_pct = 100 * (rss0 - rss) / rss0,
                dR_pct = 100 * (R_C_fg_cell_h - R0) / R0,
                dr_pct = 100 * (r - r0) / r0,
                dK_pct = 100 * (K - K0) / K0) |>
  dplyr::group_by(model) |>
  dplyr::summarise(
    n = dplyr::n(),
    curves_AIC_better = sum(dAIC < -2, na.rm = TRUE),
    curves_BIC_better = sum(dBIC < -2, na.rm = TRUE),
    median_dAIC = stats::median(dAIC, na.rm = TRUE),
    median_dBIC = stats::median(dBIC, na.rm = TRUE),
    median_rss_reduction_pct = stats::median(rss_reduction_pct, na.rm = TRUE),
    median_dR_pct = stats::median(dR_pct, na.rm = TRUE),
    median_dr_pct = stats::median(dr_pct, na.rm = TRUE),
    median_dK_pct = stats::median(dK_pct, na.rm = TRUE),
    .groups = "drop")
d6_write(C3, "C3_vs_M0.csv")

# ---- the extra parameters: are they identified? -----------------------------
C4 <- C1 |> dplyr::filter(model != "M0", converged) |>
  dplyr::group_by(model) |>
  dplyr::summarise(
    parameter = dplyr::first(dplyr::case_when(
      model %in% c("M1") ~ "tlag (min)",
      model %in% c("M2o") ~ "Km (normalised O2)",
      model %in% c("M2g") ~ "Nmax_rel",
      model %in% c("M3o") ~ "Km (normalised O2); tlag",
      TRUE ~ "Nmax_rel; tlag")),
    median_extra1 = stats::median(extra1, na.rm = TRUE),
    min_extra1 = min(extra1, na.rm = TRUE), max_extra1 = max(extra1, na.rm = TRUE),
    median_extra2 = stats::median(extra2, na.rm = TRUE),
    pct_at_bounds = 100 * mean(nzchar(at_bounds)),
    .groups = "drop")
d6_write(C4, "C4_extra_parameters.csv")

print(as.data.frame(C2), row.names = FALSE, digits = 4)
cat("\n")
print(as.data.frame(C3), row.names = FALSE, digits = 4)
cat("\n")
print(as.data.frame(C4), row.names = FALSE, digits = 4)
message("D6 PART C done.")
