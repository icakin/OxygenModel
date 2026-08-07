# =============================================================================
# d6_common.R - shared helpers for the D6 model-structure scripts
# =============================================================================
# DIAGNOSTIC + ALTERNATIVE FITS. Nothing here writes to results/, data/,
# scripts/config.R or scripts/original_scripts/. Every output goes to
# runs/D6_analysis/.
#
# The curve key is discovered from the data (Taxon [x Temperature] x Replicate),
# so the same code runs on a Candida-style temperature series unchanged.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(minpack.lm); library(tidyr)
})

.d6_this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE)) else getwd()
})
D6_BASE <- normalizePath(file.path(.d6_this, "..", "..", ".."), mustWork = FALSE)
.this_dir <- file.path(D6_BASE, "scripts")     # so config.R detects the real root
source(file.path(D6_BASE, "scripts", "config.R"))
stopifnot(normalizePath(base_dir) == normalizePath(D6_BASE))

D6_OUT <- file.path(D6_BASE, "runs", "D6_analysis")
dir.create(D6_OUT, showWarnings = FALSE, recursive = TRUE)
d6f <- function(name) file.path(D6_OUT, name)
d6_write <- function(x, name) {
  readr::write_csv(x, d6f(name))
  message(sprintf("  wrote %-54s %4d rows", name, nrow(x)))
  invisible(x)
}

ID_CANDIDATES <- c("Taxon", "Temperature", "Replicate")
id_cols <- function(df) intersect(ID_CANDIDATES, names(df))

# ---- the models --------------------------------------------------------------
# M0 - the pipeline's model, exactly as config.R defines it.
#   O2*(t) = O2_0 + (K/r)(1 - exp(r t))
m0_model <- function(r, K, t, O2_0) O2_0 + (K / r) * (1 - exp(r * t))

# M1 - M0 plus an explicit lag/equilibration offset.
#   NUISANCE, NOT BIOLOGY. The paper's methods describe a ~45 min stabilisation
#   during which the optode settles and the reported value drifts upward. A lag
#   here absorbs that instrumental transient: consumption is taken to begin at
#   t = tlag, and before that the trace is flat at O2_0.
#   O2*(t) = O2_0 + (K/r)(1 - exp(r max(t - tlag, 0)))
m1_model <- function(r, K, t, O2_0, tlag) {
  te <- pmax(t - tlag, 0)
  O2_0 + (K / r) * (1 - exp(r * te))
}

# M2o - M0 plus OXYGEN LIMITATION (Michaelis-Menten in O2).
#   Respiration becomes O2-dependent as O2 falls:
#     dO2/dt = -K exp(r t) * O2 / (Km + O2)
#   No closed form, so it is integrated numerically on the observation grid.
#   Prediction: the deviation scales with ABSOLUTE O2 and appears in any curve
#   that gets low enough, regardless of taxon.
m2o_predict <- function(r, K, t, O2_0, Km) {
  n <- length(t); y <- numeric(n); y[1] <- O2_0
  for (i in seq_len(n - 1)) {
    h <- t[i + 1] - t[i]
    f <- function(tt, o2) -K * exp(r * tt) * max(o2, 0) / (Km + max(o2, 0))
    k1 <- f(t[i],           y[i])
    k2 <- f(t[i] + h / 2,   y[i] + h * k1 / 2)
    k3 <- f(t[i] + h / 2,   y[i] + h * k2 / 2)
    k4 <- f(t[i] + h,       y[i] + h * k3)
    y[i + 1] <- y[i] + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
  }
  y
}

# M2g - M0 plus GROWTH SATURATION (logistic biomass).
#   N(t)/N0 = Nmax_rel / (1 + (Nmax_rel - 1) exp(-r t)), so consumption is
#     dO2/dt = -K * N(t)/N0
#   integrated in closed form:
#     O2(t) = O2_0 - (K Nmax_rel / r) * log( (Nmax_rel - 1 + exp(r t)) / Nmax_rel )
#   Prediction: the deviation scales with ELAPSED GROWTH (r t), not with O2.
m2g_model <- function(r, K, t, O2_0, Nmax_rel) {
  O2_0 - (K * Nmax_rel / r) * log((Nmax_rel - 1 + exp(r * t)) / Nmax_rel)
}

# M3 - lag + the supported saturation mechanism.
m3o_predict <- function(r, K, t, O2_0, Km, tlag)
  m2o_predict(r, K, pmax(t - tlag, 0), O2_0, Km)
m3g_model <- function(r, K, t, O2_0, Nmax_rel, tlag)
  m2g_model(r, K, pmax(t - tlag, 0), O2_0, Nmax_rel)

# ---- fitting ----------------------------------------------------------------
# M0 refit, exactly as 05_oxygen_fits.R does (bounds and control from config.R).
fit_m0 <- function(t, y, start) {
  d <- data.frame(t = t, y = y)
  try(minpack.lm::nlsLM(y ~ m0_model(r, K, t, O2_0), data = d, start = start,
        lower = FIT_LOWER, upper = FIT_UPPER,
        control = minpack.lm::nls.lm.control(maxiter = 1024)), silent = TRUE)
}

# Generic bounded nlsLM wrapper for the extended models.
fit_try <- function(formula, data, start, lower, upper, maxiter = 1024) {
  try(minpack.lm::nlsLM(formula, data = data, start = start,
        lower = lower, upper = upper,
        control = minpack.lm::nls.lm.control(maxiter = maxiter)), silent = TRUE)
}

# Did any parameter finish at (or within tol of) a bound? Identification check.
at_bounds <- function(fit, lower, upper, tol = 1e-6) {
  if (inherits(fit, "try-error") || is.null(fit)) return(NA_character_)
  cf <- stats::coef(fit)
  hit <- character(0)
  for (p in names(cf)) {
    if (p %in% names(lower) && abs(cf[[p]] - lower[[p]]) <= tol * max(1, abs(lower[[p]])))
      hit <- c(hit, paste0(p, "@lower"))
    if (p %in% names(upper) && abs(cf[[p]] - upper[[p]]) <= tol * max(1, abs(upper[[p]])))
      hit <- c(hit, paste0(p, "@upper"))
  }
  if (length(hit)) paste(hit, collapse = ",") else ""
}

ic <- function(fit, n) {
  if (inherits(fit, "try-error") || is.null(fit)) return(c(AIC = NA, BIC = NA, rss = NA, k = NA))
  k <- length(stats::coef(fit))
  rss <- sum(stats::residuals(fit)^2)
  c(AIC = unname(stats::AIC(fit)), BIC = unname(stats::BIC(fit)), rss = rss, k = k)
}

# Per-cell R from a fit, on the SAME footing as 05_oxygen_fits.R:
#   C_tot = (K/r)(e^{rT}-1) O2_ref ; biomass = N0 (e^{rT}-1)/r ; R = C_tot/biomass
# For the saturating model the consumed-oxygen integral differs, so C_tot is
# taken from the model's own cumulative consumption over the window.
percell_R <- function(C_tot_norm, O2_ref, N0, r, T_end) {
  C_tot <- C_tot_norm * O2_ref
  biomass <- N0 * (exp(r * T_end) - 1) / r
  if (!is.finite(C_tot) || !is.finite(biomass) || biomass <= 0) return(NA_real_)
  (C_tot / biomass) * O2_to_C_mass * MG_TO_FG * MIN_TO_H
}
