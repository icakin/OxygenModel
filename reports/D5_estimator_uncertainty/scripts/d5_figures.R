# =============================================================================
# D5 figures - built ONLY from the artefact CSVs in runs/D5_analysis/
# =============================================================================
# Nothing here recomputes an analysis; each panel renders a CSV, so the report
# and the artefacts cannot drift apart.
# OUTPUT  reports/D5_estimator_uncertainty/figures/*.png|pdf
# =============================================================================

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d5_common.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(tidyr) })

FIGD <- file.path(D5_BASE, "reports", "D5_estimator_uncertainty", "figures")
dir.create(FIGD, showWarnings = FALSE, recursive = TRUE)
sv <- function(p, name, w, h) {
  ggplot2::ggsave(file.path(FIGD, paste0(name, ".png")), p, width = w, height = h,
                  dpi = 200, bg = "white")
  message("  wrote figures/", name, ".png")
}
th <- ggplot2::theme_classic(base_size = 11)

A1 <- readr::read_csv(d5f("A1_per_curve_residual_diagnostics.csv"), show_col_types = FALSE)
A5 <- readr::read_csv(d5f("A5_residual_acf_by_lag.csv"),            show_col_types = FALSE)
B1 <- readr::read_csv(d5f("B1_within_se_vs_between_replicate_by_taxon.csv"), show_col_types = FALSE)
C1 <- readr::read_csv(d5f("C1_ridge_per_curve.csv"),                show_col_types = FALSE)
C2 <- readr::read_csv(d5f("C2_ridge_observed_replicates.csv"),      show_col_types = FALSE)
C6 <- readr::read_csv(d5f("C6_R_propagation_per_curve.csv"),        show_col_types = FALSE)
D2 <- readr::read_csv(d5f("D2_misspecification_summary.csv"),       show_col_types = FALSE)
E1 <- readr::read_csv(d5f("E1_joint_posterior_by_variant.csv"),     show_col_types = FALSE)
F3 <- readr::read_csv(d5f("F3_contrast_counts.csv"),                show_col_types = FALSE)

# ---- Fig 1: the autocorrelation that is not there ---------------------------
p1a <- ggplot(A1, aes(rho_ar1)) +
  geom_histogram(bins = 25, fill = "grey35") +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0.98, colour = "#C0392B", linewidth = 0.9) +
  annotate("text", x = 0.98, y = Inf, label = " hypothesised rho = 0.98",
           colour = "#C0392B", hjust = 0, vjust = 1.6, size = 3.2) +
  coord_cartesian(xlim = c(-0.35, 1.05)) +
  labs(x = expression("AR(1) "*rho*" of the fit residuals"), y = "curves",
       title = "A. Residual autocorrelation, all 75 curves") + th
p1b <- ggplot(A5, aes(lag, median_acf)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "grey80") +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_hline(aes(yintercept =  white_noise_band), colour = "#2980B9", linetype = 3) +
  geom_hline(aes(yintercept = -white_noise_band), colour = "#2980B9", linetype = 3) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.4) +
  labs(x = "lag (minutes)", y = "residual ACF",
       title = "B. ACF by lag (median, 5-95%)",
       subtitle = "dotted = white-noise band") + th
sv(p1a | p1b, "fig1_autocorrelation", 10, 4)

# ---- Fig 2: within-curve SE against replicate spread ------------------------
b <- B1 |> dplyr::mutate(Taxon = forcats::fct_reorder(Taxon, between_sd_logr))
p2 <- ggplot(b) +
  geom_segment(aes(y = Taxon, yend = Taxon, x = within_se_logr_naive,
                   xend = between_sd_logr), colour = "grey70", linewidth = 1) +
  geom_point(aes(y = Taxon, x = within_se_logr_naive, colour = "within-curve se(log r)"), size = 2.4) +
  geom_point(aes(y = Taxon, x = between_sd_logr, colour = "between-replicate sd(log r)"), size = 2.4) +
  scale_x_log10() +
  scale_colour_manual(values = c("within-curve se(log r)" = "#2980B9",
                                 "between-replicate sd(log r)" = "#C0392B"), name = NULL) +
  labs(x = "log scale (log10 axis)", y = NULL,
       title = "Within-curve standard error vs replicate-to-replicate spread",
       subtitle = "median ratio 37x; measurement noise is ~0.07% of replicate variance") +
  th + theme(legend.position = "top")
sv(p2, "fig2_se_vs_reproducibility", 8.5, 5)

# ---- Fig 3: the ridge, three ways, and what it does to R --------------------
rid <- dplyr::bind_rows(
  tibble::tibble(method = "vcov(nls)",          corr = C1$corr_vcov),
  tibble::tibble(method = "bootstrap (iid)",    corr = C1$corr_boot_iid),
  tibble::tibble(method = "bootstrap (AR(1))",  corr = C1$corr_boot_ar1),
  tibble::tibble(method = "real replicates",    corr = C2$corr_observed_replicates)) |>
  dplyr::filter(is.finite(corr)) |>
  dplyr::mutate(method = factor(method, levels = c("vcov(nls)", "bootstrap (iid)",
                                                   "bootstrap (AR(1))", "real replicates")))
p3a <- ggplot(rid, aes(method, corr)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.size = 0.7, fill = "grey90") +
  labs(x = NULL, y = expression(corr(log~r, log~K)),
       title = "A. The r-K ridge, measured three ways") +
  th + theme(axis.text.x = element_text(angle = 20, hjust = 1))
p3b <- ggplot(C6, aes(sd_logR_if_independent, sd_logR)) +
  geom_abline(slope = 1, linetype = 2, colour = "grey50") +
  geom_point(alpha = 0.7, colour = "#C0392B") +
  labs(x = "sd(log R) if r and K were independent", y = "actual sd(log R)",
       title = "B. The ridge PROTECTS per-cell R",
       subtitle = "every curve below the 1:1 line; median ratio 0.39") + th
sv(p3a | p3b, "fig3_ridge", 10, 4.6)

# ---- Fig 4: misspecification dwarfs everything else -------------------------
d <- D2 |> dplyr::mutate(arm = factor(arm, levels = c("none","drift","lag","saturating")))
p4 <- ggplot(d, aes(arm, median_pct_bias_R, fill = arm)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = sprintf("%+.1f%%", median_pct_bias_R),
                vjust = ifelse(median_pct_bias_R >= 0, -0.4, 1.3)), size = 3.4) +
  scale_fill_manual(values = c(none = "grey60", drift = "#2980B9",
                               lag = "#E67E22", saturating = "#C0392B"), guide = "none") +
  labs(x = NULL, y = "median bias in per-cell R (%)",
       title = "Model misspecification, 810 fits per arm",
       subtitle = "recovery under the true model is +0.09%; a lag or a saturating phase is not") +
  th
sv(p4, "fig4_misspecification", 8, 4.4)

# ---- Fig 5: joint posterior SDs are flat across taxa, under every variant ---
e <- E1 |> dplyr::mutate(variant = factor(variant,
        levels = c("naive","ar1","inflate10","inflate100")))
p5 <- ggplot(e, aes(sd_lr, Taxon, colour = variant)) +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_colour_manual(values = c(naive = "grey30", ar1 = "#2980B9",
                                 inflate10 = "#E67E22", inflate100 = "#C0392B"),
                      name = expression(S[i]*" variant")) +
  labs(x = expression("posterior sd of taxon mean "*log~r), y = NULL,
       title = "Taxon-level posterior SDs stay flat under every measurement covariance",
       subtitle = "even inflating S_i by 100x does not widen them") +
  th + theme(legend.position = "top")
sv(p5, "fig5_joint_variants", 9, 5.5)

# ---- Fig 6: the downstream count --------------------------------------------
f <- F3 |> dplyr::mutate(model = factor(model, levels = model[order(resolved_holm)]))
p6 <- ggplot(f, aes(resolved_holm, model)) +
  geom_col(fill = "grey40", width = 0.6) +
  geom_text(aes(label = sprintf("%d / %d  (%.0f%%)", resolved_holm, n_contrasts,
                                pct_resolved_holm)), hjust = -0.08, size = 3.4) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.32)), limits = c(0, 105)) +
  labs(x = "pairwise taxon contrasts resolved (Holm, 95%)", y = NULL,
       title = "Which contrasts survive, under each uncertainty model") + th
sv(p6, "fig6_contrasts", 9, 3.2)

message("D5 figures done.")
