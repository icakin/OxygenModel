# =============================================================================
# d8_common.R - shared helpers for the D8 balance analysis
# =============================================================================
# Analysis on existing outputs. Fits no new oxygen curves, changes no estimator,
# alters no committed number. Everything writes to runs/D8_analysis/.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(tidyr)
})

.d8_this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE)) else getwd()
})
D8_BASE <- normalizePath(file.path(.d8_this, "..", "..", ".."), mustWork = FALSE)
.this_dir <- file.path(D8_BASE, "scripts")
source(file.path(D8_BASE, "scripts", "config.R"))
stopifnot(normalizePath(base_dir) == normalizePath(D8_BASE))

D8_OUT <- file.path(D8_BASE, "runs", "D8_analysis")
dir.create(D8_OUT, showWarnings = FALSE, recursive = TRUE)
d8f <- function(n) file.path(D8_OUT, n)
d8_write <- function(x, n) {
  readr::write_csv(x, d8f(n))
  message(sprintf("  wrote %-48s %4d rows", n, nrow(x)))
  invisible(x)
}

BOLTZ_K <- 8.617333262e-5     # eV / K
to_K    <- function(Tc) Tc + 273.15
arr_x   <- function(Tc) -1 / (BOLTZ_K * to_K(Tc))   # ln(rate) = ln(a) + E * arr_x

# ---- 11_temperature_cue.R's OWN constants, which override config.R ----------
# Line ~39-40 of that script. config.R sets INOC_DELAY_MIN <- 0; the temperature
# script sets 45. This is a KNOWN INCONSISTENCY and is handled explicitly
# everywhere below - never silently reconciled.
INOC_CELLS_PER_uL_11 <- 550
INOC_DELAY_MIN_11    <- 45
INOC_DELAY_MIN_CFG   <- INOC_DELAY_MIN     # from config.R (0)

# ---- Boltzmann-Arrhenius on the RISING limb ---------------------------------
# Fits ln(rate) = ln(a) + E * (-1/kT) over points up to and including the
# observed peak, which is the standard activation-energy estimate. Returns E in
# eV with a 95% interval.
ba_fit <- function(Tc, rate, label = "", rising_only = TRUE) {
  ok <- is.finite(Tc) & is.finite(rate) & rate > 0
  Tc <- Tc[ok]; rate <- rate[ok]
  if (rising_only) {
    agg <- tapply(rate, Tc, stats::median)
    Tpk <- as.numeric(names(agg))[which.max(agg)]
    keep <- Tc <= Tpk
  } else keep <- rep(TRUE, length(Tc))
  x <- arr_x(Tc[keep]); y <- log(rate[keep])
  m <- stats::lm(y ~ x); ci <- stats::confint(m)
  tibble::tibble(
    quantity = label, n = sum(keep), T_peak_C = if (rising_only) Tpk else NA_real_,
    E_eV = unname(stats::coef(m)[2]), E_lo = ci[2, 1], E_hi = ci[2, 2],
    p = summary(m)$coefficients[2, 4], r_squared = summary(m)$r.squared)
}

# ---- argmin of a ratio against temperature ----------------------------------
# A quadratic in T on the log scale is the cheapest smooth interpolant that has
# a single interior optimum; the vertex gives the argmin. Bootstrapped over
# curves for an interval.
argmin_quad <- function(Tc, ratio, nboot = 2000, seed = 20260808) {
  f <- function(Tv, rv) {
    m <- stats::lm(log(rv) ~ poly(Tv, 2, raw = TRUE))
    b <- unname(stats::coef(m))
    if (length(b) < 3 || !is.finite(b[3]) || b[3] <= 0) return(NA_real_)
    -b[2] / (2 * b[3])
  }
  pt <- f(Tc, ratio)
  set.seed(seed)
  bs <- replicate(nboot, {
    i <- sample(seq_along(Tc), replace = TRUE)
    if (length(unique(Tc[i])) < 3) return(NA_real_)
    f(Tc[i], ratio[i])
  })
  bs <- bs[is.finite(bs)]
  c(argmin = pt,
    lo = if (length(bs)) unname(stats::quantile(bs, .025)) else NA_real_,
    hi = if (length(bs)) unname(stats::quantile(bs, .975)) else NA_real_,
    n_boot_ok = length(bs))
}
