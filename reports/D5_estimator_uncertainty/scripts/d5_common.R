# =============================================================================
# d5_common.R - shared helpers for the D5 diagnostic scripts
# =============================================================================
# Sourced by every d5_part*.R script. Defines the project paths, the curve-level
# refit used throughout, and the AR(1) machinery.
#
# DIAGNOSTIC ONLY. Nothing here writes to results/, data/ or scripts/. Every
# output goes to runs/D5_analysis/.
#
# MULTI-DATASET NOTE: the grouping key is discovered from the data
# (Taxon [x Temperature] x Replicate), never hard-coded, so the same code runs
# on a Candida-style temperature series without change.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(minpack.lm)
})

.d5_this <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE)) else getwd()
})
D5_BASE <- normalizePath(file.path(.d5_this, "..", "..", ".."), mustWork = FALSE)

# config.R auto-detects the project root as dirname(.this_dir) when .this_dir
# exists, and otherwise from its own call frame - which, sourced from here,
# resolves to this report's folder and would create a stray results/ tree under
# reports/. Point it at the real scripts/ directory first.
.this_dir <- file.path(D5_BASE, "scripts")
source(file.path(D5_BASE, "scripts", "config.R"))
stopifnot(normalizePath(base_dir) == normalizePath(D5_BASE))

D5_OUT <- file.path(D5_BASE, "runs", "D5_analysis")
dir.create(D5_OUT, showWarnings = FALSE, recursive = TRUE)
d5f <- function(name) file.path(D5_OUT, name)

d5_write <- function(x, name) {
  readr::write_csv(x, d5f(name))
  message(sprintf("  wrote %-52s %4d rows", name, nrow(x)))
  invisible(x)
}

# ---- grouping key ------------------------------------------------------------
# The curve identifier. Taxon x Replicate here; Taxon x Temperature x Replicate
# in a temperature-gradient dataset. Discovered, not assumed.
ID_CANDIDATES <- c("Taxon", "Temperature", "Replicate")
id_cols <- function(df) intersect(ID_CANDIDATES, names(df))
curve_id <- function(df, cols = id_cols(df)) do.call(paste, c(df[cols], sep = " | "))

# ---- the model ---------------------------------------------------------------
# Identical to config.R's resp_model; restated so a reader of this script does
# not have to hold config.R in their head.
#   O2_norm(t) = O2_0 + (K / r) * (1 - exp(r * t))
d5_model <- function(r, K, t, O2_0) O2_0 + (K / r) * (1 - exp(r * t))

# ---- ordinary (naive) refit, exactly as 05_oxygen_fits.R does ----------------
# 05 starts from a heuristic; we start from its converged answer, which lands on
# the same optimum and avoids re-deriving the heuristic. Bounds, control and the
# model are 05's.
fit_ols <- function(t, y, start) {
  d <- data.frame(t = t, y = y)
  try(minpack.lm::nlsLM(
        y ~ d5_model(r, K, t, O2_0), data = d,
        start = start, lower = FIT_LOWER, upper = FIT_UPPER,
        control = minpack.lm::nls.lm.control(maxiter = 1024)),
      silent = TRUE)
}

# ---- AR(1) helpers -----------------------------------------------------------
# rho by lag-1 regression through the origin on the residual series, which is
# the standard AR(1) moment estimator and is what the variance-inflation
# expression below assumes.
ar1_rho <- function(e) {
  n <- length(e)
  if (n < 3) return(NA_real_)
  num <- sum(e[-1] * e[-n]); den <- sum(e[-n]^2)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  max(min(num / den, 0.9999), -0.9999)
}

durbin_watson <- function(e) sum(diff(e)^2) / sum(e^2)

# Effective sample size and the variance-inflation factor for a mean-like
# statistic under AR(1): sd inflates by sqrt((1+rho)/(1-rho)).
n_eff_ar1 <- function(n, rho) n * (1 - rho) / (1 + rho)
vif_ar1   <- function(rho) sqrt((1 + rho) / (1 - rho))

# ---- AR(1) generalised least squares refit (Prais-Winsten) -------------------
# WHY THIS RATHER THAN Newey-West AS THE HEADLINE CORRECTION: with rho near 1 a
# HAC estimator needs a bandwidth comparable to the correlation length, and the
# Bartlett kernel converges slowly there - at rho = 0.98 the correlation length
# is ~50 points against n ~ 110, so a rule-of-thumb bandwidth of 4 corrects
# almost nothing. Modelling the AR(1) structure explicitly uses it directly.
# Newey-West is computed alongside, at a stated bandwidth, as a cross-check.
#
# The transform: y*_1 = sqrt(1-rho^2) y_1, y*_t = y_t - rho y_{t-1}, and the
# mean function is transformed the same way, so the transformed problem has
# independent errors and its vcov IS the GLS covariance.
fit_ar1_gls <- function(t, y, start, rho_init = NULL, iters = 3L) {
  n <- length(y)
  f0 <- fit_ols(t, y, start)
  if (inherits(f0, "try-error")) return(NULL)
  e  <- y - stats::predict(f0)
  rho <- if (is.null(rho_init)) ar1_rho(e) else rho_init
  fit <- f0
  for (i in seq_len(iters)) {
    if (!is.finite(rho) || abs(rho) >= 0.9999) break
    w1 <- sqrt(1 - rho^2)
    yt <- c(w1 * y[1], y[-1] - rho * y[-n])
    modf <- function(r, K, O2_0) {
      f <- d5_model(r, K, t, O2_0)
      c(w1 * f[1], f[-1] - rho * f[-n])
    }
    d <- data.frame(yt = yt)
    ft <- try(minpack.lm::nlsLM(
                yt ~ modf(r, K, O2_0), data = d,
                start = as.list(stats::coef(fit)),
                lower = FIT_LOWER, upper = FIT_UPPER,
                control = minpack.lm::nls.lm.control(maxiter = 1024)),
              silent = TRUE)
    if (inherits(ft, "try-error")) break
    fit <- ft
    e   <- y - d5_model(coef(fit)["r"], coef(fit)["K"], t, coef(fit)["O2_0"])
    rho <- ar1_rho(e)
  }
  list(fit = fit, rho = rho, resid = e)
}

# ---- Newey-West (HAC) sandwich for a nonlinear least-squares fit -------------
# Var = (J'J)^-1 [ S_0 + sum_l w_l (S_l + S_l') ] (J'J)^-1, Bartlett weights
# w_l = 1 - l/(L+1), score contributions u_t = J_t e_t.
# Default bandwidth is the usual Newey-West rule L = floor(4 (n/100)^(2/9));
# callers may pass a larger, correlation-length-matched bandwidth.
nw_vcov <- function(J, e, L = NULL) {
  n <- nrow(J); p <- ncol(J)
  if (is.null(L)) L <- max(1L, floor(4 * (n / 100)^(2 / 9)))
  u <- J * e
  meat <- crossprod(u)
  if (L >= 1) for (l in seq_len(L)) {
    w <- 1 - l / (L + 1)
    S <- crossprod(u[-seq_len(l), , drop = FALSE], u[seq_len(n - l), , drop = FALSE])
    meat <- meat + w * (S + t(S))
  }
  bread <- tryCatch(solve(crossprod(J)), error = function(err) NULL)
  if (is.null(bread)) return(NULL)
  V <- bread %*% meat %*% bread
  dimnames(V) <- list(colnames(J), colnames(J))
  V
}

# Analytic Jacobian of d5_model wrt (r, K, O2_0).
model_jacobian <- function(r, K, t, O2_0) {
  E <- exp(r * t)
  cbind(r    = -(K / r^2) * (1 - E) - (K / r) * t * E,
        K    = (1 / r) * (1 - E),
        O2_0 = rep(1, length(t)))
}
