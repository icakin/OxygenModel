# =============================================================================
# 08_window_sensitivity.R - Reviewer robustness test: how much do K, R and
# per-cell respiration depend on the manual fit window?
# =============================================================================
# manual_fit_windows.csv sets both edges on all 75 curves and the choice moves
# absolute K and R. This refits every curve over a grid of edge shifts and reports
# (a) the per-curve spread in K, R and r, and (b) whether a systematic shift
# changes the between-taxon ordering. Reuses the pipeline model/constants via
# config.R. Non-destructive: writes its own table + figure, changes no pipeline
# number.
#   INPUT   results/tables/Oxygen_Data_Filtered.csv, manual_fit_windows.csv, data/Ninoc.csv
#   OUTPUT  results/tables/window_sensitivity_percurve.csv
#           results/figures/Fig_window_sensitivity.png
#   RUN     Rscript scripts/08_window_sensitivity.R
# =============================================================================

# ---- locate scripts/ dir and source shared config ---------------------------
.this_dir <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  d <- if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE)) else
         tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (is.null(d) || is.na(d) || !nzchar(d)) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
        nzchar(rstudioapi::getActiveDocumentContext()$path))
      d <- dirname(rstudioapi::getActiveDocumentContext()$path) else d <- getwd()
  }
  d
})
source(file.path(.this_dir, "config.R"))
suppressPackageStartupMessages({ library(dplyr); library(minpack.lm); library(ggplot2)
                                 library(tidyr); library(patchwork) })

ox  <- readr::read_csv(FILTERED_CSV, show_col_types = FALSE) |>
         mutate(Taxon = as.character(Taxon), Replicate = as.character(Replicate))
win <- readr::read_csv(file.path(tables_dir, "manual_fit_windows.csv"), show_col_types = FALSE)
nin <- readr::read_csv(NINOC_CSV, show_col_types = FALSE)

# fit one curve over [start,end]; return r, K, C_tot, T_end, O0, QC flag
fit_one <- function(df0, start, end) {
  df <- df0 |> arrange(Time) |> filter(Time >= start, Time <= end)
  if (nrow(df) < 5) return(NULL)
  df$Time0 <- df$Time - min(df$Time)
  O0 <- mean(head(df$Oxygen, 3)); if (!is.finite(O0) || O0 <= 0) return(NULL)
  df$Oxygen_norm <- df$Oxygen / O0
  T_end <- max(df$Time0); if (!is.finite(T_end) || T_end <= 0) return(NULL)
  seg <- head(df, max(3, floor(0.3 * nrow(df))))
  r0  <- min(max(max(abs(diff(log(pmax(seg$Oxygen_norm,1e-6)))/diff(seg$Time0)),na.rm=TRUE),1e-4),5e-2)
  k0  <- min(max(abs(suppressWarnings(min(diff(df$Oxygen_norm)/diff(df$Time0),na.rm=TRUE))),1e-5),0.1)
  fit <- tryCatch(nlsLM(Oxygen_norm ~ resp_model(r,K,Time0,O2_0), data = df,
                        start = list(r=r0,K=k0,O2_0=1), lower=FIT_LOWER, upper=FIT_UPPER,
                        control = nls.lm.control(maxiter=300)), error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  p <- summary(fit)$coefficients
  r_est <- p["r","Estimate"]; K_est <- p["K","Estimate"]
  pred <- predict(fit, df); resid <- df$Oxygen_norm - pred
  R2 <- 1 - sum(resid^2)/sum((df$Oxygen_norm-mean(df$Oxygen_norm))^2)
  ok <- all(abs(p[,"Std. Error"]/p[,"Estimate"]) < REL_SE_THRESHOLD, R2 >= R2_THRESHOLD,
            is.finite(K_est), K_est>0, K_est<0.5, is.finite(r_est), r_est>=1e-4, r_est<=0.1)
  list(r=r_est, K=K_est, C_tot=(K_est/r_est)*(exp(r_est*T_end)-1)*O0, T_end=T_end, O0=O0, ok=ok)
}
percell_R <- function(f, N_inoc, delta) {           # fg C / cell / h
  N0 <- N_inoc * exp(f$r * delta)
  (f$C_tot / (N0 * (exp(f$r * f$T_end) - 1) / f$r)) * O2_to_C_mass * MG_TO_FG * MIN_TO_H
}

OFF <- c(-6,-3,0,3,6)                                 # minute offsets per edge
rows <- list()
for (i in seq_len(nrow(win))) {
  Tax <- win$Taxon[i]; Rep <- win$Replicate[i]
  df0 <- ox  |> filter(Taxon==Tax, Replicate==Rep)
  nr  <- nin |> filter(Taxon==Tax, Replicate==Rep)
  if (nrow(df0) < 5 || nrow(nr) == 0) next
  N_inoc <- nr$N_inoculation_cells_per_L[1]; delta <- nr$delta_Ninoc_to_N0_min[1]
  base <- fit_one(df0, win$fit_start[i], win$fit_end[i]); if (is.null(base)) next
  baseR <- percell_R(base, N_inoc, delta)
  for (ds in OFF) for (de in OFF) {
    f <- fit_one(df0, win$fit_start[i]+ds, win$fit_end[i]+de); if (is.null(f)) next
    rows[[length(rows)+1]] <- tibble::tibble(
      Taxon=Tax, Replicate=Rep, ds=ds, de=de, is_base=(ds==0 & de==0),
      r=f$r, K=f$K, Rpc=percell_R(f,N_inoc,delta), ok=f$ok,
      dK_pct=100*(f$K-base$K)/base$K, dR_pct=100*(percell_R(f,N_inoc,delta)-baseR)/baseR,
      dr_pct=100*(f$r-base$r)/base$r)
  }
}
res <- dplyr::bind_rows(rows)
readr::write_csv(res, file.path(tables_dir, "window_sensitivity_percurve.csv"))

# ---- (a) per-curve jitter ----------------------------------------------------
pc <- res |> group_by(Taxon,Replicate) |>
  summarise(r=max(abs(dr_pct)), K=max(abs(dK_pct)), R=max(abs(dR_pct)), fok=mean(ok), .groups="drop")
qs <- function(x) sprintf("median %.1f%% | 90th %.1f%% | max %.1f%%", median(x), quantile(x,.9), max(x))
cat("\n===== PER-CURVE WINDOW SENSITIVITY (", nrow(pc), "curves, both edges +/-6 min ) =====\n")
cat("  K          :", qs(pc$K), "\n"); cat("  per-cell R :", qs(pc$R), "\n")
cat("  r          :", qs(pc$r), "\n")
cat(sprintf("  QC pass rate across perturbed windows: median %.0f%%\n", 100*median(pc$fok)))

# ---- (b) systematic shift ----------------------------------------------------
base_tax <- res |> filter(is_base) |> group_by(Taxon) |> summarise(base=mean(Rpc), .groups="drop")
cat("\n===== SYSTEMATIC SHIFT (all windows moved together) =====\n")
for (cc in list(c(-6,-6),c(6,6),c(-6,6),c(6,-6))) {
  s <- res |> filter(ds==cc[1], de==cc[2]) |> group_by(Taxon) |> summarise(v=mean(Rpc), .groups="drop") |>
       inner_join(base_tax, by="Taxon")
  cat(sprintf("  start%+d end%+d: taxon-mean R shifts median %.1f%%, max %.1f%%; ordering rho=%.3f\n",
              cc[1], cc[2], median(abs(100*(s$v-s$base)/s$base)), max(abs(100*(s$v-s$base)/s$base)),
              suppressWarnings(cor(s$v, s$base, method="spearman"))))
}

# ---- figure ------------------------------------------------------------------
pcl <- pc |> dplyr::select(Taxon,Replicate,r,K,R) |> pivot_longer(c(r,K,R), names_to="param", values_to="pct")
pcl$param <- factor(pcl$param, levels=c("r","K","R"), labels=c("r","K","per-cell R"))
pA <- ggplot(pcl, aes(param, pct, fill=param)) + geom_boxplot(width=.55, outlier.size=.7, alpha=.85) +
  scale_fill_manual(values=c("r"="#4C9F70","K"="#E1A730","per-cell R"="#C0504D"), guide="none") +
  labs(title="A  Per-curve window jitter (both edges +/-6 min)", subtitle="max |% change| per curve, n=75",
       x=NULL, y="max |% change| vs chosen window") + theme_classic(base_size=12)
tj <- inner_join(base_tax, res |> filter(ds==6,de==6) |> group_by(Taxon) |> summarise(shifted=mean(Rpc),.groups="drop"), by="Taxon")
pB <- ggplot(tj, aes(base, shifted)) + geom_abline(slope=1,intercept=0,linetype=2,colour="grey60") +
  geom_point(size=2.6, colour="#31688E", alpha=.85) +
  labs(title="B  Between-taxon ordering under a systematic +6/+6 shift",
       subtitle=sprintf("taxon-mean per-cell R; Spearman rho = %.3f", cor(tj$base,tj$shifted,method="spearman")),
       x="baseline taxon-mean per-cell R (fg C / cell / h)", y="shifted") + theme_classic(base_size=12)
ggsave(file.path(figures_dir, "Fig_window_sensitivity.png"), pA/pB, width=8.2, height=8.8, dpi=150, bg="white")
cat("\nWrote window_sensitivity_percurve.csv and Fig_window_sensitivity.png\n")
