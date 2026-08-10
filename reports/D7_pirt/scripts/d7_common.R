# =============================================================================
# d7_common.R - shared helpers for the D7 Pirt analysis
# =============================================================================
# ANALYSIS ON EXISTING OUTPUTS ONLY. Fits no oxygen curves, changes no estimator,
# alters no committed number. Everything writes to runs/D7_analysis/.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(tidyr)
})

.d7_this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE)) else getwd()
})
D7_BASE <- normalizePath(file.path(.d7_this, "..", "..", ".."), mustWork = FALSE)
.this_dir <- file.path(D7_BASE, "scripts")
source(file.path(D7_BASE, "scripts", "config.R"))
stopifnot(normalizePath(base_dir) == normalizePath(D7_BASE))

D7_OUT <- file.path(D7_BASE, "runs", "D7_analysis")
dir.create(D7_OUT, showWarnings = FALSE, recursive = TRUE)
d7f <- function(n) file.path(D7_OUT, n)
d7_write <- function(x, n) {
  readr::write_csv(x, d7f(n))
  message(sprintf("  wrote %-50s %4d rows", n, nrow(x)))
  invisible(x)
}

# ---- Pirt fit: R = m + G/Y, respiration as RESPONSE, LINEAR axes -------------
# Returns the intercept (maintenance m), slope (1/Y), their 95% intervals, R^2,
# the implied yield Y, and the extrapolation leverage - how far the intercept at
# G = 0 sits outside the observed range of G, in units of that range.
pirt_fit <- function(G, R, label = "") {
  ok <- is.finite(G) & is.finite(R)
  G <- G[ok]; R <- R[ok]
  m <- stats::lm(R ~ G)
  ci <- stats::confint(m)
  s  <- summary(m)
  gmin <- min(G); gmax <- max(G); grange <- gmax - gmin
  slope <- unname(stats::coef(m)[2])
  tibble::tibble(
    dataset = label, n = length(G),
    m_intercept = unname(stats::coef(m)[1]),
    m_lo = ci[1, 1], m_hi = ci[1, 2],
    m_p = s$coefficients[1, 4],
    slope_inv_Y = slope,
    slope_lo = ci[2, 1], slope_hi = ci[2, 2],
    slope_p = s$coefficients[2, 4],
    Y = if (is.finite(slope) && slope > 0) 1 / slope else NA_real_,
    Y_lo = if (is.finite(ci[2,2]) && ci[2,2] > 0) 1 / ci[2, 2] else NA_real_,
    Y_hi = if (is.finite(ci[2,1]) && ci[2,1] > 0) 1 / ci[2, 1] else NA_real_,
    r_squared = s$r.squared,
    G_min = gmin, G_max = gmax,
    # how far below the observed data does G = 0 sit, in observed-range units
    extrapolation_ranges_to_zero = if (grange > 0) gmin / grange else NA_real_,
    G_min_over_G_max = gmin / gmax,
    corr_G_R = stats::cor(G, R))
}
