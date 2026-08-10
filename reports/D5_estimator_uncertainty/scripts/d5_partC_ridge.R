# =============================================================================
# D5 PART C - the r-K ridge, measured three ways
# =============================================================================
#   (i)   corr(log r, log K) from vcov(nls) per curve - what the pipeline sees
#   (ii)  a parametric bootstrap from each fitted curve at its own residual
#         scale, ONCE with i.i.d. noise and ONCE with AR(1) noise
#   (iii) the observed within-taxon correlation across the five real replicates
#
# Plus the rotated basis: the sampling variance of (log r + log K) against
# (log K - log r), to say which combination is determined and which is not.
#
# vcov is on the natural scale; the delta method takes it to the log scale:
#   cov(log r, log K) = cov(r, K) / (r K), so the CORRELATION is unchanged.
# The rotated-basis variances, however, need the log-scale covariance, so that
# is what is used there.
#
# INPUT   results/tables/oxygen_fit_curves.csv, oxygen_fit_results.csv
#         runs/D5_analysis/A1_per_curve_residual_diagnostics.csv
# OUTPUT  runs/D5_analysis/C1_*.csv .. C5_*.csv
# =============================================================================

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d5_common.R"))

set.seed(20260807)          # bootstrap reproducibility
N_BOOT <- 300L              # stated in the report; a diagnostic, not a precision exercise

message("\n== D5 PART C: the r-K ridge ==")
message("  parametric bootstrap replicates: ", N_BOOT, " (seed 20260807)")

curves <- readr::read_csv(FITCURVES_CSV,   show_col_types = FALSE)
fits   <- readr::read_csv(RESULTS_FIT_CSV, show_col_types = FALSE)
A1     <- readr::read_csv(d5f("A1_per_curve_residual_diagnostics.csv"), show_col_types = FALSE)
ID     <- id_cols(fits)

# Simulate an AR(1) series with marginal sd sigma.
ar1_noise <- function(n, rho, sigma) {
  if (!is.finite(rho) || abs(rho) < 1e-8) return(stats::rnorm(n, 0, sigma))
  as.numeric(stats::arima.sim(list(ar = rho), n = n,
                              sd = sigma * sqrt(1 - rho^2)))
}

boot_rows <- vector("list", nrow(fits))
per_curve <- vector("list", nrow(fits))

for (i in seq_len(nrow(fits))) {
  key <- fits[i, ID, drop = FALSE]
  g <- curves |> dplyr::semi_join(key, by = ID) |> dplyr::arrange(Time0_min)
  if (nrow(g) < 10) next
  t <- g$Time0_min; n <- length(t)
  a <- A1 |> dplyr::semi_join(key, by = ID)
  if (!nrow(a)) next
  r0 <- a$r[1]; K0 <- a$K[1]; O0 <- a$O2_0[1]
  sig <- a$sigma[1]; rho <- a$rho_ar1[1]
  start <- list(r = r0, K = K0, O2_0 = O0)
  mu <- d5_model(r0, K0, t, O0)

  # ---- (i) vcov correlation, and the log-scale covariance -----------------
  f <- fit_ols(t, g$Oxygen_norm, start)
  V <- if (inherits(f, "try-error")) NULL else tryCatch(stats::vcov(f), error = function(e) NULL)
  corr_vcov <- NA_real_; var_sum <- NA_real_; var_dif <- NA_real_
  if (!is.null(V) && all(c("r", "K") %in% rownames(V))) {
    corr_vcov <- V["r", "K"] / sqrt(V["r", "r"] * V["K", "K"])
    # delta method to log scale
    Vl <- matrix(c(V["r","r"]/r0^2, V["r","K"]/(r0*K0),
                   V["r","K"]/(r0*K0), V["K","K"]/K0^2), 2, 2)
    var_sum <- Vl[1,1] + Vl[2,2] + 2*Vl[1,2]   # var(log r + log K)
    var_dif <- Vl[1,1] + Vl[2,2] - 2*Vl[1,2]   # var(log K - log r)
  }

  # ---- (ii) parametric bootstrap, iid and AR(1) ---------------------------
  bs <- function(noise_fun) {
    lr <- numeric(0); lk <- numeric(0)
    for (b in seq_len(N_BOOT)) {
      yb <- mu + noise_fun()
      fb <- fit_ols(t, yb, start)
      if (inherits(fb, "try-error")) next
      cb <- stats::coef(fb)
      if (!all(is.finite(cb)) || cb["r"] <= 0 || cb["K"] <= 0) next
      lr <- c(lr, log(cb["r"])); lk <- c(lk, log(cb["K"]))
    }
    if (length(lr) < 20) return(c(NA, NA, NA, NA))
    c(stats::cor(lr, lk), stats::var(lr + lk), stats::var(lk - lr), length(lr))
  }
  b_iid <- bs(function() stats::rnorm(n, 0, sig))
  b_ar1 <- bs(function() ar1_noise(n, rho, sig))

  per_curve[[i]] <- tibble::tibble(
    key, n = n, sigma = sig, rho_ar1 = rho,
    corr_vcov = corr_vcov,
    corr_boot_iid = b_iid[1], corr_boot_ar1 = b_ar1[1],
    var_logsum_vcov = var_sum, var_logdif_vcov = var_dif,
    var_logsum_boot_iid = b_iid[2], var_logdif_boot_iid = b_iid[3],
    var_logsum_boot_ar1 = b_ar1[2], var_logdif_boot_ar1 = b_ar1[3],
    n_boot_ok_iid = b_iid[4], n_boot_ok_ar1 = b_ar1[4])
  if (i %% 15 == 0) message("    ", i, "/", nrow(fits), " curves bootstrapped")
}

C1 <- dplyr::bind_rows(per_curve) |>
  dplyr::mutate(ridge_ratio_vcov = var_logdif_vcov / var_logsum_vcov,
                ridge_ratio_boot_iid = var_logdif_boot_iid / var_logsum_boot_iid,
                ridge_ratio_boot_ar1 = var_logdif_boot_ar1 / var_logsum_boot_ar1)
d5_write(C1, "C1_ridge_per_curve.csv")

# ---- (iii) observed within-taxon correlation across real replicates --------
GRP <- setdiff(ID, "Replicate")
C2 <- A1 |>
  dplyr::group_by(dplyr::across(dplyr::all_of(GRP))) |>
  dplyr::summarise(
    n_rep = dplyr::n(),
    corr_observed_replicates = if (dplyr::n() >= 3)
      stats::cor(log(r), log(K)) else NA_real_,
    var_logsum_replicates = stats::var(log(r) + log(K)),
    var_logdif_replicates = stats::var(log(K) - log(r)),
    .groups = "drop") |>
  dplyr::mutate(ridge_ratio_replicates = var_logdif_replicates / var_logsum_replicates)
d5_write(C2, "C2_ridge_observed_replicates.csv")

q <- function(x) { x <- x[is.finite(x)]
  c(n = length(x), min = min(x), median = stats::median(x), max = max(x)) }
C3 <- dplyr::bind_rows(
  tibble::tibble(method = "(i)   vcov(nls) per curve",              !!!q(C1$corr_vcov)),
  tibble::tibble(method = "(ii)  parametric bootstrap, iid noise",  !!!q(C1$corr_boot_iid)),
  tibble::tibble(method = "(ii)  parametric bootstrap, AR(1) noise",!!!q(C1$corr_boot_ar1)),
  tibble::tibble(method = "(iii) observed across real replicates",  !!!q(C2$corr_observed_replicates))
)
d5_write(C3, "C3_ridge_correlation_three_ways.csv")

C4 <- dplyr::bind_rows(
  tibble::tibble(basis = "var(log r + log K)", source = "vcov",            !!!q(C1$var_logsum_vcov)),
  tibble::tibble(basis = "var(log K - log r)", source = "vcov",            !!!q(C1$var_logdif_vcov)),
  tibble::tibble(basis = "var(log r + log K)", source = "bootstrap iid",   !!!q(C1$var_logsum_boot_iid)),
  tibble::tibble(basis = "var(log K - log r)", source = "bootstrap iid",   !!!q(C1$var_logdif_boot_iid)),
  tibble::tibble(basis = "var(log r + log K)", source = "bootstrap AR(1)", !!!q(C1$var_logsum_boot_ar1)),
  tibble::tibble(basis = "var(log K - log r)", source = "bootstrap AR(1)", !!!q(C1$var_logdif_boot_ar1)),
  tibble::tibble(basis = "var(log r + log K)", source = "real replicates", !!!q(C2$var_logsum_replicates)),
  tibble::tibble(basis = "var(log K - log r)", source = "real replicates", !!!q(C2$var_logdif_replicates))
)
d5_write(C4, "C4_rotated_basis_variances.csv")

C5 <- tibble::tibble(
  quantity = c("median corr(log r, log K), vcov",
               "median corr(log r, log K), bootstrap iid",
               "median corr(log r, log K), bootstrap AR(1)",
               "median corr(log r, log K), real replicates",
               "median var(log K - log r) / var(log r + log K), vcov",
               "median var(log K - log r) / var(log r + log K), bootstrap iid",
               "median var(log K - log r) / var(log r + log K), real replicates",
               "bootstrap replicates per curve"),
  value = c(stats::median(C1$corr_vcov, na.rm = TRUE),
            stats::median(C1$corr_boot_iid, na.rm = TRUE),
            stats::median(C1$corr_boot_ar1, na.rm = TRUE),
            stats::median(C2$corr_observed_replicates, na.rm = TRUE),
            stats::median(C1$ridge_ratio_vcov, na.rm = TRUE),
            stats::median(C1$ridge_ratio_boot_iid, na.rm = TRUE),
            stats::median(C2$ridge_ratio_replicates, na.rm = TRUE),
            N_BOOT))
d5_write(C5, "C5_ridge_headline.csv")

message(sprintf("\n  corr(log r, log K): vcov %.3f | boot-iid %.3f | boot-AR1 %.3f | replicates %.3f",
  stats::median(C1$corr_vcov, na.rm = TRUE), stats::median(C1$corr_boot_iid, na.rm = TRUE),
  stats::median(C1$corr_boot_ar1, na.rm = TRUE), stats::median(C2$corr_observed_replicates, na.rm = TRUE)))
message(sprintf("  ridge ratio var(logK-logr)/var(logr+logK): vcov %.2f | boot %.2f | replicates %.2f",
  stats::median(C1$ridge_ratio_vcov, na.rm = TRUE),
  stats::median(C1$ridge_ratio_boot_iid, na.rm = TRUE),
  stats::median(C2$ridge_ratio_replicates, na.rm = TRUE)))
message("D5 PART C done.")
