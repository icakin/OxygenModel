# =============================================================================
# 14_joint_hierarchical_prototype.R   (v2 - reparameterised)
#   PROTOTYPE: one hierarchical model fitted to the oxygen traces directly,
#   instead of the current two-stage "fit each curve, then summarise".
# =============================================================================
# WHY THIS EXISTS
#   The pipeline fits every curve on its own with nlsLM, keeps the point
#   estimates of r and K, and carries those into a separate second model.
#   Within a single fit r and K are correlated at about -0.985 (PART 1 below):
#   the data constrains their COMBINATION far better than either alone. Stage
#   one keeps only coefficients and SEs and discards vcov(), so that structure
#   never reaches R = K * O2_ref / N0, which depends on both.
#
# WHY v2 IS PARAMETERISED DIFFERENTLY
#   v1 estimated (r, K) directly and would not converge: R-hat 3.7, with
#   hundreds of max-treedepth hits. That is the sampler crawling along the same
#   near-degenerate r-K ridge described above. The cure is to estimate a pair
#   that is NOT degenerate:
#
#       O2(t) = O2_0 - D * (exp(r*t) - 1) / (exp(r*T) - 1)
#
#   where T is the curve's own window length and D is the TOTAL drawdown over
#   that window. D is read almost directly off the data (median 0.56, log-sd
#   0.28) and r is then just the curve shape, so the two are close to
#   independent. This is algebraically the same model: afterwards
#
#       K = D * r / (exp(r*T) - 1)
#
#   recovers the original parameter exactly, which PART 3 does on the posterior.
#
# WHAT IT DOES
#   PART 1  diagnostic: per-curve r-K correlation, and fit-SE vs replicate scatter
#   PART 2  the joint hierarchical fit (brms), reparameterised
#   PART 3  comparison: joint vs two-stage uncertainty on taxon-level r and K
#
# INPUT   results/tables/oxygen_fit_curves.csv     (from 05_oxygen_fits.R)
#         results/tables/oxygen_results_with_R.csv
# OUTPUT  results/rds/joint_hierarchical_fit.rds
#         results/tables/joint_vs_twostage_uncertainty.csv
#         results/tables/per_curve_rK_covariance.csv
#
# READ THE CONVERGENCE LINE BEFORE THE RATIOS. If max R-hat > 1.01 the
# comparison is meaningless no matter how sensible the numbers look.
#
# REQUIRES brms + Stan backend, minpack.lm, dplyr, posterior.
# RUN      Rscript scripts/14_joint_hierarchical_prototype.R
# =============================================================================

QUICK <- TRUE     # TRUE = 4 chains x 1500. FALSE = 4 chains x 4000.

# ---- locate scripts/ and shared config --------------------------------------
.this_dir <- if (
  requireNamespace("rstudioapi", quietly = TRUE) &&
  rstudioapi::isAvailable() &&
  nzchar(rstudioapi::getActiveDocumentContext()$path)
) {
  dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1])))
  else tryCatch(dirname(sys.frame(1)$ofile), error = function(e) getwd())
}
source(file.path(.this_dir, "config.R"))
suppressPackageStartupMessages({
  library(dplyr); library(brms); library(minpack.lm); library(posterior)
})

curves <- readr::read_csv(tbl("oxygen_fit_curves.csv"),     show_col_types = FALSE)
res    <- readr::read_csv(tbl("oxygen_results_with_R.csv"), show_col_types = FALSE)

# =============================================================================
# PART 1 - what does the two-stage approach discard?
# =============================================================================
message("\n[1/3] Per-curve r-K covariance diagnostic ...")
resp_model <- function(r, K, t, O2_0) O2_0 + (K / r) * (1 - exp(r * t))

cov_rows <- list()
for (i in seq_len(nrow(res))) {
  tax <- res$Taxon[i]; rp <- res$Replicate[i]
  df <- curves %>% filter(Taxon == tax, Replicate == rp) %>% arrange(Time0_min)
  if (nrow(df) < 6) next
  fit <- try(nlsLM(Oxygen_norm ~ resp_model(r, K, Time0_min, O2_0), data = df,
                   start = list(r = res$r_per_minute[i], K = res$K[i], O2_0 = res$O2_0[i]),
                   lower = FIT_LOWER, upper = FIT_UPPER,
                   control = nls.lm.control(maxiter = 1024)), silent = TRUE)
  if (inherits(fit, "try-error")) next
  V <- try(vcov(fit), silent = TRUE)
  if (inherits(V, "try-error") || !all(c("r", "K") %in% rownames(V))) next
  cf <- coef(fit); se_r <- sqrt(V["r","r"]); se_K <- sqrt(V["K","K"])
  cov_rows[[length(cov_rows)+1]] <- data.frame(
    Taxon = tax, Replicate = rp, r = cf["r"], K = cf["K"],
    se_r = se_r, se_K = se_K, rho_rK = V["r","K"]/(se_r*se_K))
}
cov_df <- bind_rows(cov_rows)
readr::write_csv(cov_df, tbl("per_curve_rK_covariance.csv"))
message(sprintf("      within-fit cor(r, K): median %.3f  (n = %d curves)",
                median(cov_df$rho_rK, na.rm = TRUE), nrow(cov_df)))

# =============================================================================
# PART 2 - joint hierarchical fit, reparameterised as (r, D)
# =============================================================================
message("\n[2/3] Fitting joint hierarchical model ...")
d <- curves %>%
  filter(is.finite(Time0_min), is.finite(Oxygen_norm)) %>%
  group_by(Taxon, Replicate) %>%
  mutate(Twin = max(Time0_min, na.rm = TRUE)) %>%   # per-curve window length
  ungroup() %>%
  transmute(Taxon = factor(Taxon), Replicate,
            curve = factor(paste(Taxon, Replicate, sep = "_")),
            Time0_min, Twin, Oxygen_norm)

# O2(t) = O20 - D * (exp(r t) - 1)/(exp(r T) - 1)
f <- bf(Oxygen_norm ~ O20 - exp(logD) *
          (exp(exp(logr)*Time0_min) - 1) / (exp(exp(logr)*Twin) - 1),
        O20  ~ 1,                                   # O2_0 is ~1.00 in all curves
        logr ~ 1 + (1 | Taxon) + (1 | curve),
        logD ~ 1 + (1 | Taxon) + (1 | curve),
        nl = TRUE)

pr <- c(prior(normal(-3.8, 0.5),  nlpar = "logr"),
        prior(normal(-0.58, 0.4), nlpar = "logD"),
        prior(normal(1, 0.02),    nlpar = "O20"),
        prior(exponential(10),    class = "sd", nlpar = "logr"),
        prior(exponential(10),    class = "sd", nlpar = "logD"),
        prior(exponential(1000),  class = "sigma"))  # data wants sigma ~0.001

ITER <- if (QUICK) 1500 else 4000
CH   <- 4

t0 <- Sys.time()
fit <- brm(f, data = d, prior = pr,
           chains = CH, iter = ITER, warmup = ITER/2,
           cores = min(CH, parallel::detectCores()),
           control = list(adapt_delta = 0.90, max_treedepth = 10),
           init = function() list(b_logr = as.array(-3.8),
                                  b_logD = as.array(-0.58),
                                  b_O20  = as.array(1),
                                  sigma  = 0.002),   # keep chains out of the
                                                     # "noise explains everything" basin

           init_r = 0.2,      # keep ALL chains near sensible starts. brms
                              # scatters unnamed params widely by default, which
                              # is what let one chain start in a hopeless region
                              # and never recover (sigma 0.184 vs 0.001).
           refresh = max(1, ITER/20), seed = 1)
message(sprintf("      done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
saveRDS(fit, file.path(rds_dir, "joint_hierarchical_fit.rds"))

rh  <- max(rhat(fit), na.rm = TRUE)
div <- sum(subset(nuts_params(fit), Parameter == "divergent__")$Value)
message(sprintf("      max R-hat = %.4f | divergences = %d", rh, div))
if (rh > 1.01) {
  message("\n  *** NOT CONVERGED (R-hat > 1.01). The ratios below are NOT usable. ***")
  message("  Try QUICK <- FALSE, and if it still fails send me the R-hat and treedepth warnings.\n")
} else {
  message("      converged.\n")
}

# =============================================================================
# PART 3 - joint vs two-stage uncertainty, per taxon
# =============================================================================
message("[3/3] Comparing joint vs two-stage uncertainty ...")
Tref <- median(d$Twin, na.rm = TRUE)

dr   <- as_draws_df(fit)
taxa <- levels(d$Taxon)
get <- function(p, tx) {
  cn <- sprintf("r_Taxon__%s[%s,Intercept]", p, tx)
  if (cn %in% names(dr)) dr[[cn]] else 0
}
joint <- lapply(taxa, function(tx) {
  logr <- dr[["b_logr_Intercept"]] + get("logr", tx)
  logD <- dr[["b_logD_Intercept"]] + get("logD", tx)
  r <- exp(logr); D <- exp(logD)
  K <- D * r / (exp(r * Tref) - 1)          # back to the original parameter
  data.frame(Taxon = tx, sd_logr_joint = sd(logr), sd_logK_joint = sd(log(K)))
}) %>% bind_rows()

twostage <- res %>% group_by(Taxon) %>%
  summarise(sd_logr_twostage = sd(log(r_per_minute), na.rm = TRUE)/sqrt(n()),
            sd_logK_twostage = sd(log(K), na.rm = TRUE)/sqrt(n()), .groups = "drop")

cmp <- left_join(twostage, joint, by = "Taxon") %>%
  mutate(ratio_logr = sd_logr_joint/sd_logr_twostage,
         ratio_logK = sd_logK_joint/sd_logK_twostage)
readr::write_csv(cmp, tbl("joint_vs_twostage_uncertainty.csv"))
print(as.data.frame(cmp), row.names = FALSE, digits = 3)
message(sprintf(
  "\n  median uncertainty ratio (joint / two-stage): log r %.2fx, log K %.2fx",
  median(cmp$ratio_logr, na.rm = TRUE), median(cmp$ratio_logK, na.rm = TRUE)))
message("  >1 means the two-stage approach was over-confident.")
message(sprintf("  [valid only if max R-hat <= 1.01; this run: %.3f]\n", rh))
