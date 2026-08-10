# =============================================================================
# D7 PARTS C + D - reconcile with Figure 6, and quantify identifiability
# =============================================================================
# The gate in PART A failed, so no maintenance or yield estimate is carried
# forward. Two things remain worth establishing:
#
#  (C) How the published log-log slope of 0.283 relates to the Pirt
#      parameterisation - i.e. what Fig 6 already implies about maintenance, and
#      whether that agrees with the direct fit. If the two parameterisations of
#      the same data disagree, that locates the problem.
#
#  (D) Whether the split is identifiable AT ALL at this design, independent of
#      the shared-r artefact: simulate from a KNOWN Pirt relationship at the
#      observed n, spread of G and residual scale, and check recovery and the
#      intercept-slope sampling correlation. Then report the spread in G that
#      WOULD be needed.
#
# OUTPUT  runs/D7_analysis/C1_*.csv, D1_*.csv .. D3_*.csv
# =============================================================================

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d7_common.R"))
suppressPackageStartupMessages({ library(lme4); library(lmerTest) })

set.seed(20260808)
message("\n== D7 PARTS C + D ==")

res <- readr::read_csv(RESULTS_FINAL_CSV, show_col_types = FALSE) |>
  dplyr::filter(fit_ok, is.finite(G_C_fg_cell_h), G_C_fg_cell_h > 0,
                is.finite(R_C_fg_cell_h), R_C_fg_cell_h > 0)

# ---------------------------------------------------------------------------
# PART C - what the 0.283 log-log slope implies about maintenance
# ---------------------------------------------------------------------------
# Fig 6 fits log_G ~ log_Rc with (1 | Taxon) - GROWTH as the response. Refit here
# with the SAME random-effects structure so only the parameterisation differs,
# then invert.
rm_df <- res |>
  dplyr::mutate(Taxon = factor(Taxon),
                log_G = log10(G_C_fg_cell_h), log_R = log10(R_C_fg_cell_h))
x0 <- mean(rm_df$log_R)
rm_df$log_Rc <- rm_df$log_R - x0
mm_pub <- lmerTest::lmer(log_G ~ log_Rc + (1 | Taxon), data = rm_df, REML = FALSE)
slope_pub <- unname(lme4::fixef(mm_pub)[2])

# The reverse regression, same hierarchy - this is the elasticity Pirt speaks to.
y0 <- mean(rm_df$log_G)
rm_df$log_Gc <- rm_df$log_G - y0
mm_rev <- lmerTest::lmer(log_R ~ log_Gc + (1 | Taxon), data = rm_df, REML = FALSE)
slope_rev <- unname(lme4::fixef(mm_rev)[2])

# Under Pirt R = m + G/Y, the elasticity is
#   d log R / d log G = (G/Y) / (m + G/Y)
# which lies in (0, 1] for ANY m >= 0, approaching 1 as m -> 0. An elasticity
# ABOVE 1 is impossible with non-negative maintenance, and implies
#   m = (G/Y) (1/elasticity - 1) < 0.
elasticity_from_pub <- 1 / slope_pub
Gbar <- mean(res$G_C_fg_cell_h); Rbar <- mean(res$R_C_fg_cell_h)
# Solve directly: elasticity e = (G/Y)/R at the mean, so G/Y = e*R and
# m = R - G/Y = R(1 - e).
implied_m_pub <- Rbar * (1 - elasticity_from_pub)
implied_m_rev <- Rbar * (1 - slope_rev)

A3 <- readr::read_csv(d7f("A3_gate_pirt_fits.csv"), show_col_types = FALSE)
m_direct <- A3$m_intercept[A3$dataset == "N0 default (r-dependent)"]

C1 <- tibble::tibble(
  quantity = c("published Fig 6 slope, log_G ~ log_Rc, (1|Taxon)",
               "refit here with the same structure",
               "implied elasticity d log R / d log G = 1/slope",
               "reverse regression log_R ~ log_Gc, (1|Taxon)",
               "max elasticity possible with m >= 0",
               "mean per-cell R (fg C / cell / h)",
               "implied maintenance m from the Fig 6 elasticity",
               "implied maintenance m from the reverse regression",
               "maintenance m from the DIRECT linear Pirt fit (PART A)"),
  value = c(0.283, slope_pub, elasticity_from_pub, slope_rev, 1.0,
            Rbar, implied_m_pub, implied_m_rev, m_direct))
d7_write(C1, "C1_fig6_reconciliation.csv")

# ---------------------------------------------------------------------------
# PART D - identifiability at this design, setting the artefact aside
# ---------------------------------------------------------------------------
# Simulate from a KNOWN Pirt relationship using the OBSERVED G values and the
# observed residual scale, then refit. This asks the cleanest possible question:
# even if R on G were pure physiology, could m and Y be recovered here?
G_obs <- res$G_C_fg_cell_h
fit_obs <- stats::lm(R_C_fg_cell_h ~ G_C_fg_cell_h, data = res)
sigma_obs <- stats::sigma(fit_obs)

sim_recover <- function(G, m_true, Y_true, sigma, nsim = 2000) {
  slope_true <- 1 / Y_true
  out <- matrix(NA_real_, nsim, 2)
  for (b in seq_len(nsim)) {
    R <- m_true + G * slope_true + stats::rnorm(length(G), 0, sigma)
    cf <- stats::coef(stats::lm(R ~ G))
    out[b, ] <- cf
  }
  tibble::tibble(
    n = length(G),
    G_min = min(G), G_max = max(G), G_spread_ratio = max(G) / min(G),
    m_true = m_true, Y_true = Y_true,
    m_mean = mean(out[, 1]), m_sd = stats::sd(out[, 1]),
    m_bias = mean(out[, 1]) - m_true,
    slope_mean = mean(out[, 2]), slope_sd = stats::sd(out[, 2]),
    corr_intercept_slope = stats::cor(out[, 1], out[, 2]),
    pct_m_negative = 100 * mean(out[, 1] < 0),
    m_rel_precision = stats::sd(out[, 1]) / abs(m_true))
}

# A plausible truth: maintenance ~20% of mean observed R, yield 0.3.
m_true <- 0.2 * Rbar; Y_true <- 0.3
D1 <- sim_recover(G_obs, m_true, Y_true, sigma_obs)
D1$scenario <- "observed design (15 taxa, 30 C)"

# What spread in G would be needed? Rescale the observed G about its mean to
# widen the range, holding n and sigma fixed.
widen <- function(G, factor) {
  mu <- mean(G); mu + (G - mu) * factor
}
D2 <- dplyr::bind_rows(lapply(c(1, 2, 4, 8, 16), function(f) {
  Gw <- widen(G_obs, f)
  Gw <- Gw - min(Gw) + max(min(G_obs) / f, 1e-6)   # keep strictly positive
  s <- sim_recover(Gw, m_true, Y_true, sigma_obs, nsim = 1000)
  s$widen_factor <- f
  s
}))
d7_write(D2, "D2_identifiability_vs_spread.csv")
d7_write(D1, "D1_identifiability_observed.csv")

# How much MORE data would the observed design need? sd(m) scales as 1/sqrt(n),
# so to reach a target relative precision the required n scales as the square of
# the ratio. Reported as a concrete, quotable number.
target_rel <- 0.20
n_needed <- ceiling(D1$n * (D1$m_rel_precision / target_rel)^2)

D3 <- tibble::tibble(
  quantity = c("observed G spread, max/min",
               "intercept-slope sampling correlation at the observed design",
               "relative precision of m at the observed design (sd/|m|)",
               "% of simulated datasets giving a NEGATIVE m",
               "residual sd used (from the observed fit)",
               "assumed true m (20% of mean R)",
               "assumed true Y",
               "n needed for 20% relative precision on m, at this noise level",
               "n actually available"),
  value = c(D1$G_spread_ratio, D1$corr_intercept_slope, D1$m_rel_precision,
            D1$pct_m_negative, sigma_obs, m_true, Y_true, n_needed, D1$n))
d7_write(D3, "D3_identifiability_headline.csv")

message("\n  PART C - Fig 6 reconciliation:")
print(as.data.frame(C1), row.names = FALSE, digits = 4)
message("\n  PART D - identifiability at the observed design:")
print(as.data.frame(D1[, c("n", "G_spread_ratio", "corr_intercept_slope",
                           "m_bias", "m_sd", "m_rel_precision", "pct_m_negative")]),
      row.names = FALSE, digits = 4)
message("\n  PART D - how it improves with spread in G:")
print(as.data.frame(D2[, c("widen_factor", "G_spread_ratio", "corr_intercept_slope",
                           "m_rel_precision", "pct_m_negative")]),
      row.names = FALSE, digits = 4)
message("\nD7 PARTS C+D done.")
