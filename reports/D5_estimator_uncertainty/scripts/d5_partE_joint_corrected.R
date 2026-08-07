# =============================================================================
# D5 PART E - the joint estimator under corrected measurement covariance
# =============================================================================
# 12_joint_rK_estimator.R is NOT modified. This is a copy of its two stages that
# swaps the per-curve measurement covariance S_i and writes to new files.
#
# Four variants of S_i, the same Stan model, data, seed and control settings:
#
#   naive       vcov(nls), exactly what 12 uses          -> reproduces 12
#   ar1         the AR(1) GLS covariance from PART A     -> the brief's ask
#   inflate10   naive x 10                                \  counterfactual: how
#   inflate100  naive x 100                               /  large would S_i have
#                                                            to be to matter?
#
# PART A found rho ~ 0, so `ar1` is expected to be almost identical to `naive`.
# The inflation arms are what actually answer the question D6 needs answered -
# whether propagating the within-curve covariance can EVER matter at this design,
# or whether the between-replicate variance Tau dominates regardless.
#
# OUTPUT  runs/D5_analysis/E1_*.csv .. E3_*.csv
# =============================================================================

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d5_common.R"))
suppressPackageStartupMessages(library(rstan))
rstan_options(auto_write = TRUE); options(mc.cores = parallel::detectCores())

message("\n== D5 PART E: joint estimator under four S_i variants ==")

curves <- readr::read_csv(FITCURVES_CSV,     show_col_types = FALSE)
res    <- readr::read_csv(RESULTS_FINAL_CSV, show_col_types = FALSE)
A1     <- readr::read_csv(d5f("A1_per_curve_residual_diagnostics.csv"), show_col_types = FALSE)
ID     <- id_cols(res)

# ---- Stage 1: per-curve (log r, log K) with BOTH covariances ----------------
message("[1/2] refitting curves, keeping naive and AR(1) covariances ...")
rows <- list()
for (i in seq_len(nrow(res))) {
  key <- res[i, ID, drop = FALSE]
  df <- curves |> dplyr::semi_join(key, by = ID) |> dplyr::arrange(Time0_min)
  if (nrow(df) < 6) next
  start <- list(r = res$r_per_minute[i], K = res$K[i], O2_0 = res$O2_0[i])

  f <- fit_ols(df$Time0_min, df$Oxygen_norm, start)
  if (inherits(f, "try-error")) next
  V <- try(stats::vcov(f), silent = TRUE); if (inherits(V, "try-error")) next
  if (!all(c("r", "K") %in% rownames(V))) next
  cf <- stats::coef(f); r <- cf["r"]; K <- cf["K"]

  a <- A1 |> dplyr::semi_join(key, by = ID)
  g <- fit_ar1_gls(df$Time0_min, df$Oxygen_norm, start,
                   rho_init = if (nrow(a)) a$rho_ar1[1] else NULL)
  Vg <- if (is.null(g)) NULL else try(stats::vcov(g$fit), silent = TRUE)
  if (is.null(Vg) || inherits(Vg, "try-error")) Vg <- V
  cg <- if (is.null(g)) cf else stats::coef(g$fit)

  rows[[length(rows) + 1]] <- data.frame(
    key,
    lr = log(r), lk = log(K),
    v_lr = V["r","r"]/r^2, v_lk = V["K","K"]/K^2, c_lrk = V["r","K"]/(r*K),
    v_lr_ar1 = Vg["r","r"]/cg["r"]^2, v_lk_ar1 = Vg["K","K"]/cg["K"]^2,
    c_lrk_ar1 = Vg["r","K"]/(cg["r"]*cg["K"]))
}
pc <- dplyr::bind_rows(rows)
message(sprintf("      %d curves kept", nrow(pc)))
d5_write(tibble::as_tibble(pc), "E0_per_curve_covariances.csv")

# ---- Stage 2: the same Stan model, four times -------------------------------
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
sm <- stan_model(model_code = stan_code)

taxa <- sort(unique(pc$Taxon)); J <- length(taxa); N <- nrow(pc)
tax  <- match(pc$Taxon, taxa)
y    <- as.matrix(pc[, c("lr", "lk")])

make_S <- function(variant) {
  S <- array(0, dim = c(N, 2, 2))
  for (i in seq_len(N)) {
    if (variant == "ar1") {
      S[i,1,1] <- pc$v_lr_ar1[i]; S[i,2,2] <- pc$v_lk_ar1[i]
      S[i,1,2] <- pc$c_lrk_ar1[i]; S[i,2,1] <- pc$c_lrk_ar1[i]
    } else {
      m <- switch(variant, naive = 1, inflate10 = 10, inflate100 = 100)
      S[i,1,1] <- m * pc$v_lr[i]; S[i,2,2] <- m * pc$v_lk[i]
      S[i,1,2] <- m * pc$c_lrk[i]; S[i,2,1] <- m * pc$c_lrk[i]
    }
  }
  S
}

VARIANTS <- c("naive", "ar1", "inflate10", "inflate100")
out <- list(); diag_rows <- list()

for (v in VARIANTS) {
  message("  fitting variant: ", v)
  fit <- sampling(sm, data = list(N = N, J = J, tax = tax, y = y, S = make_S(v)),
                  chains = 4, iter = 2000, warmup = 1000, seed = 1,
                  refresh = 0, control = list(adapt_delta = 0.95))
  post <- rstan::extract(fit, "mu")$mu
  rh <- max(summary(fit)$summary[, "Rhat"], na.rm = TRUE)
  tau_post <- rstan::extract(fit, "tau")$tau
  out[[v]] <- tibble::tibble(
    variant = v, Taxon = taxa,
    m_lr = apply(post[,,1], 2, mean), m_lk = apply(post[,,2], 2, mean),
    sd_lr = apply(post[,,1], 2, stats::sd), sd_lk = apply(post[,,2], 2, stats::sd),
    ci_width_lr = 2 * 1.96 * apply(post[,,1], 2, stats::sd),
    ci_width_lk = 2 * 1.96 * apply(post[,,2], 2, stats::sd))
  diag_rows[[v]] <- tibble::tibble(
    variant = v, max_Rhat = rh,
    sd_lr_min = min(out[[v]]$sd_lr), sd_lr_max = max(out[[v]]$sd_lr),
    sd_lr_spread_pct = 100 * (max(out[[v]]$sd_lr) / min(out[[v]]$sd_lr) - 1),
    sd_lk_min = min(out[[v]]$sd_lk), sd_lk_max = max(out[[v]]$sd_lk),
    sd_lk_spread_pct = 100 * (max(out[[v]]$sd_lk) / min(out[[v]]$sd_lk) - 1),
    tau_lr_mean = mean(tau_post[,1]), tau_lk_mean = mean(tau_post[,2]))
  message(sprintf("    max R-hat %.4f | sd(log r) %.5f-%.5f | tau_lr %.4f",
                  rh, min(out[[v]]$sd_lr), max(out[[v]]$sd_lr), mean(tau_post[,1])))
}

E1 <- dplyr::bind_rows(out)
d5_write(E1, "E1_joint_posterior_by_variant.csv")
E2 <- dplyr::bind_rows(diag_rows)
d5_write(E2, "E2_joint_variant_diagnostics.csv")

base <- dplyr::filter(E1, variant == "naive") |>
  dplyr::select(Taxon, sd_lr_base = sd_lr, sd_lk_base = sd_lk)
E3 <- E1 |> dplyr::left_join(base, by = "Taxon") |>
  dplyr::group_by(variant) |>
  dplyr::summarise(
    median_widening_lr = stats::median(sd_lr / sd_lr_base),
    median_widening_lk = stats::median(sd_lk / sd_lk_base),
    max_widening_lr = max(sd_lr / sd_lr_base),
    .groups = "drop")
d5_write(E3, "E3_joint_interval_widening.csv")

message("D5 PART E done.")
