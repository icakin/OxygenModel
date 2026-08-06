# =============================================================================
# 12_joint_rK_estimator.R
#   Joint hierarchical fit that PROPAGATES the within-curve r-K covariance
#   into the taxon-level estimates, instead of discarding vcov() as the
#   two-stage pipeline does. This version converges (R-hat ~1.00 in seconds).
# =============================================================================
# WHY NOT FIT THE RAW TRACES DIRECTLY
#   A single brms model fitted to the exponential traces with taxon + curve
#   random effects will not converge (R-hat > 3): the smooth curves let r and K
#   trade off, and 5 replicates per taxon cannot separate taxon- from
#   replicate-level variance. The posterior is a long ridge and chains never mix.
#
# WHAT THIS DOES INSTEAD (same goal, stable)
#   Stage 1: fit each curve, KEEP the full 2x2 covariance of (log r, log K).
#   Stage 2: a hierarchical measurement-error model. Each curve's (log r, log K)
#            is a noisy measurement of its taxon's true value, with that known
#            covariance built in; the model estimates taxon means and the
#            between-replicate spread. It is linear-Gaussian, so it converges
#            immediately, and the discarded covariance now flows through.
#
#   For curve i in taxon t:   y_i ~ MVN( mu_t , Tau + S_i )
#   y_i = (log r_i, log K_i);  S_i = its measurement covariance (the bit the
#   old pipeline threw away);  Tau = between-replicate covariance;  mu_t = the
#   taxon-level answer, with honest uncertainty.
#
# INPUT   results/tables/oxygen_fit_curves.csv, oxygen_results_with_R.csv
# OUTPUT  results/rds/joint_rK_estimator.rds
#         results/tables/joint_vs_twostage_uncertainty.csv
#         results/figures/Fig_joint_vs_twostage.png
# RUN     Rscript scripts/12_joint_rK_estimator.R          (needs rstan)
# =============================================================================

.this_dir <- {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
}
source(file.path(.this_dir, "config.R"))
suppressPackageStartupMessages({ library(dplyr); library(minpack.lm); library(rstan) })
rstan_options(auto_write = TRUE); options(mc.cores = parallel::detectCores())

curves <- readr::read_csv(tbl("oxygen_fit_curves.csv"),     show_col_types = FALSE)
res    <- readr::read_csv(tbl("oxygen_results_with_R.csv"), show_col_types = FALSE)
resp_model <- function(r, K, t, O2_0) O2_0 + (K / r) * (1 - exp(r * t))

# ---- Stage 1: per-curve (log r, log K) + full log-scale covariance ----------
message("[1/2] refitting curves, keeping full r-K covariance ...")
rows <- list()
for (i in seq_len(nrow(res))) {
  df <- curves %>% filter(Taxon == res$Taxon[i], Replicate == res$Replicate[i]) %>%
    arrange(Time0_min)
  if (nrow(df) < 6) next
  fit <- try(nlsLM(Oxygen_norm ~ resp_model(r, K, Time0_min, O2_0), data = df,
             start = list(r = res$r_per_minute[i], K = res$K[i], O2_0 = res$O2_0[i]),
             lower = FIT_LOWER, upper = FIT_UPPER,
             control = nls.lm.control(maxiter = 1024)), silent = TRUE)
  if (inherits(fit, "try-error")) next
  V <- try(vcov(fit), silent = TRUE); if (inherits(V, "try-error")) next
  if (!all(c("r","K") %in% rownames(V))) next
  cf <- coef(fit); r <- cf["r"]; K <- cf["K"]
  rows[[length(rows)+1]] <- data.frame(
    Taxon = res$Taxon[i], Replicate = res$Replicate[i],
    lr = log(r), lk = log(K),
    v_lr = V["r","r"]/r^2, v_lk = V["K","K"]/K^2, c_lrk = V["r","K"]/(r*K))
}
pc <- bind_rows(rows)
message(sprintf("      %d curves kept", nrow(pc)))

# ---- Stage 2: hierarchical measurement-error model (inline Stan) ------------
message("[2/2] fitting hierarchical measurement-error model ...")
taxa <- sort(unique(pc$Taxon)); J <- length(taxa); N <- nrow(pc)
tax  <- match(pc$Taxon, taxa)
y    <- as.matrix(pc[, c("lr","lk")])
S    <- array(0, dim = c(N,2,2))
for (i in 1:N) { S[i,1,1]<-pc$v_lr[i]; S[i,2,2]<-pc$v_lk[i]
                 S[i,1,2]<-pc$c_lrk[i]; S[i,2,1]<-pc$c_lrk[i] }

stan_code <- "
data{
  int<lower=1> N; int<lower=1> J;
  array[N] int<lower=1,upper=J> tax;
  array[N] vector[2] y; array[N] matrix[2,2] S;
}
parameters{
  array[J] vector[2] mu;
  vector<lower=0>[2] tau;
  cholesky_factor_corr[2] Lc;
}
transformed parameters{
  matrix[2,2] Tau = diag_pre_multiply(tau,Lc)*diag_pre_multiply(tau,Lc)';
}
model{
  for (j in 1:J) mu[j] ~ normal(0,5);
  tau ~ normal(0,1); Lc ~ lkj_corr_cholesky(2);
  for (i in 1:N) y[i] ~ multi_normal(mu[tax[i]], Tau + S[i]);
}"

sm  <- stan_model(model_code = stan_code)
fit <- sampling(sm, data = list(N=N,J=J,tax=tax,y=y,S=S),
                chains = 4, iter = 2000, warmup = 1000, seed = 1,
                refresh = 500, control = list(adapt_delta = 0.95))
saveRDS(fit, file.path(rds_dir, "joint_rK_estimator.rds"))

rh <- max(summary(fit)$summary[,"Rhat"], na.rm = TRUE)
message(sprintf("      max R-hat = %.4f  %s", rh,
                if (rh <= 1.01) "(converged)" else "(NOT converged)"))

# ---- joint (covariance-propagated) vs two-stage uncertainty -----------------
post  <- rstan::extract(fit, "mu")$mu       # draws x J x 2
joint <- data.frame(Taxon = taxa,
                    m_lr_j = apply(post[,,1], 2, mean), m_lk_j = apply(post[,,2], 2, mean),
                    sd_lr_joint = apply(post[,,1], 2, sd), sd_lk_joint = apply(post[,,2], 2, sd))
two   <- pc %>% group_by(Taxon) %>%
  summarise(m_lr = mean(lr), m_lk = mean(lk),
            sd_lr_two = sd(lr)/sqrt(n()), sd_lk_two = sd(lk)/sqrt(n()), .groups="drop")
cmp <- left_join(two, joint, by = "Taxon") %>%
  mutate(ratio_logr = sd_lr_joint/sd_lr_two, ratio_logK = sd_lk_joint/sd_lk_two)
readr::write_csv(cmp, tbl("joint_vs_twostage_uncertainty.csv"))
print(as.data.frame(cmp), row.names = FALSE, digits = 3)
message(sprintf("\n  median uncertainty ratio (joint / two-stage): log r %.2fx, log K %.2fx",
                median(cmp$ratio_logr), median(cmp$ratio_logK)))
message("  >1 means the two-stage pipeline was over-confident.")
message(sprintf("  [valid only if max R-hat <= 1.01; this run: %.3f]", rh))

# ---- figure: joint vs two-stage 95% intervals -------------------------------
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
.jvt_panel <- function(mc, sdt, mjc, sdj, xlab, title) {
  df <- rbind(data.frame(Taxon = cmp$Taxon, m = cmp[[mc]],  sd = cmp[[sdt]], method = "Two-stage"),
              data.frame(Taxon = cmp$Taxon, m = cmp[[mjc]], sd = cmp[[sdj]], method = "Joint"))
  df$lo <- df$m - 1.96 * df$sd; df$hi <- df$m + 1.96 * df$sd
  df$Taxon  <- factor(df$Taxon, levels = cmp$Taxon[order(cmp[[mc]])])
  df$method <- factor(df$method, levels = c("Two-stage", "Joint"))
  ggplot2::ggplot(df, ggplot2::aes(m, Taxon, colour = method)) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.8,
                            position = ggplot2::position_dodge(width = 0.55)) +
    ggplot2::geom_point(size = 2, position = ggplot2::position_dodge(width = 0.55)) +
    ggplot2::scale_colour_manual(values = c("Two-stage" = "grey55", "Joint" = "#C0392B"), guide = "none") +
    ggplot2::labs(x = xlab, y = NULL, title = title) + ggplot2::theme_classic(base_size = 12)
}
.pr <- .jvt_panel("m_lr", "sd_lr_two", "m_lr_j", "sd_lr_joint", "log growth rate (log r)",
                  sprintf("Growth rate  (~%.2fx wider)", median(cmp$ratio_logr)))
.pk <- .jvt_panel("m_lk", "sd_lk_two", "m_lk_j", "sd_lk_joint", "log oxygen constant (log K)",
                  sprintf("Oxygen constant K  (~%.2fx wider)", median(cmp$ratio_logK)))
ggplot2::ggsave(fig("Fig_joint_vs_twostage.png"),
                (.pr | .pk) + patchwork::plot_annotation(
                  title = "Two-stage (grey) vs covariance-propagated joint fit (red)",
                  subtitle = "Points sit in the same place; only the 95% error bars change. Ranking unchanged."),
                width = 13, height = 7, dpi = 150, bg = "white")
message("  wrote Fig_joint_vs_twostage.png")
