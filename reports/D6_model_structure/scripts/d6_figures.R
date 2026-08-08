# =============================================================================
# D6 figures - built ONLY from the artefact CSVs in runs/D6_analysis/
# =============================================================================
source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d6_common.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

FIGD <- file.path(D6_BASE, "reports", "D6_model_structure", "figures")
dir.create(FIGD, showWarnings = FALSE, recursive = TRUE)
sv <- function(p, name, w, h) {
  ggplot2::ggsave(file.path(FIGD, paste0(name, ".png")), p, width = w, height = h,
                  dpi = 200, bg = "white")
  message("  wrote figures/", name, ".png")
}
th <- ggplot2::theme_classic(base_size = 11)

A1 <- readr::read_csv(d6f("A1_residual_structure_per_curve.csv"), show_col_types = FALSE)
A2 <- readr::read_csv(d6f("A2_residual_profile_by_window_position.csv"), show_col_types = FALSE)
A3 <- readr::read_csv(d6f("A3_residual_profile_by_absolute_O2.csv"), show_col_types = FALSE)
B1 <- readr::read_csv(d6f("B1_window_sensitivity_vs_saturation.csv"), show_col_types = FALSE)
B3 <- readr::read_csv(d6f("B3_variance_explained.csv"), show_col_types = FALSE)
C3 <- readr::read_csv(d6f("C3_vs_M0.csv"), show_col_types = FALSE)
D1 <- readr::read_csv(d6f("D1_power_check.csv"), show_col_types = FALSE)
D3 <- readr::read_csv(d6f("D3_window_sweep_by_model.csv"), show_col_types = FALSE)

# ---- Fig 1: no residual structure -------------------------------------------
p1a <- ggplot(A2, aes(frac_window_mid, mean_resid_sd)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = "position within the fit window (0 = start, 1 = end)",
       y = "mean residual (in residual sd)",
       title = "A. Residuals vs window position",
       subtitle = "saturation would bend this down late") + th
p1b <- ggplot(A3, aes(o2_abs_mid, mean_resid_sd)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = "absolute oxygen (mg/L)", y = NULL,
       title = "B. Residuals vs absolute oxygen",
       subtitle = "O2 limitation would bend this at low O2") + th
sv(p1a | p1b, "fig1_residual_structure", 10, 4.2)

# ---- Fig 2: why - the windows stop early -------------------------------------
p2 <- ggplot(A1, aes(O2_frac_remaining)) +
  geom_histogram(bins = 22, fill = "grey35") +
  geom_vline(xintercept = 0.2, colour = "#C0392B", linewidth = 0.9) +
  annotate("text", x = 0.2, y = Inf, label = "  20% - below here O2 limitation could bite",
           colour = "#C0392B", hjust = 0, vjust = 1.8, size = 3.1) +
  labs(x = "fraction of starting oxygen remaining at the window end",
       y = "curves",
       title = "The trimming stops each window long before oxygen runs out",
       subtitle = "median 44.5% remaining; only 7 of 75 end below 20%") + th
sv(p2, "fig2_windows_stop_early", 8.5, 4.2)

# ---- Fig 3: the window sensitivity is arithmetic -----------------------------
p3a <- ggplot(B1, aes(100 * r * 6, max_abs_dK_pct)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(alpha = 0.75, colour = "#2980B9") +
  labs(x = expression("arithmetic prediction  100 %*% r %*% 6 min  (%)"),
       y = "observed max |dK| (%)",
       title = "A. Window sensitivity of K is r x window shift",
       subtitle = "dashed = exact first-order prediction; R² of r alone = 0.977") + th
b3 <- B3 |> dplyr::mutate(predictor = forcats::fct_reorder(predictor, R2_dK))
p3b <- ggplot(b3, aes(R2_dK, predictor, fill = family)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c(saturation = "#C0392B", control = "grey45"), name = NULL) +
  labs(x = expression("R"^2*" on max |dK| (%)"), y = NULL,
       title = "B. What predicts the sensitivity") +
  th + theme(legend.position = "top")
sv(p3a | p3b, "fig3_arithmetic", 11, 4.6)

# ---- Fig 4: no model is preferred --------------------------------------------
c3 <- C3 |> dplyr::mutate(model = factor(model, levels = c("M1","M2o","M2g","M3o","M3g")))
p4a <- ggplot(c3, aes(model, median_dAIC)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_col(fill = "grey40", width = 0.6) +
  geom_text(aes(label = sprintf("%+.1f", median_dAIC)), vjust = -0.4, size = 3.3) +
  labs(x = NULL, y = "median dAIC vs M0",
       title = "A. Every extension is worse by AIC",
       subtitle = "positive = penalised for a parameter that buys no fit") + th
p4b <- ggplot(c3, aes(model, curves_AIC_better)) +
  geom_col(fill = "grey40", width = 0.6) +
  geom_text(aes(label = paste0(curves_AIC_better, "/75")), vjust = -0.4, size = 3.3) +
  coord_cartesian(ylim = c(0, 75)) +
  labs(x = NULL, y = "curves where dAIC < -2",
       title = "B. Curves actually preferring the extension") + th
sv(p4a | p4b, "fig4_model_comparison", 10, 4.2)

# ---- Fig 5: the power check --------------------------------------------------
d1 <- D1 |>
  tidyr::pivot_longer(c(pct_M1_preferred, pct_M2g_preferred),
                      names_to = "model", values_to = "pct") |>
  dplyr::mutate(model = ifelse(model == "pct_M1_preferred", "M1 (lag)", "M2g (saturation)"),
                arm = factor(arm, levels = c("none", "lag", "saturating")))
p5 <- ggplot(d1, aes(arm, pct, fill = model)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  geom_text(aes(label = sprintf("%.0f%%", pct)),
            position = position_dodge(width = 0.72), vjust = -0.35, size = 3.2) +
  scale_fill_manual(values = c("M1 (lag)" = "#E67E22", "M2g (saturation)" = "#C0392B"),
                    name = NULL) +
  coord_cartesian(ylim = c(0, 112)) +
  labs(x = "simulated truth", y = "% of curves where the extension wins on AIC",
       title = "The models CAN detect these effects when present",
       subtitle = "so the null on real data is 'no effect', not 'no power'") +
  th + theme(legend.position = "top")
sv(p5, "fig5_power_check", 8.5, 4.4)

# ---- Fig 6: the key test - sweep under each model ----------------------------
d3 <- D3 |>
  tidyr::pivot_longer(c(median_max_dK_pct, median_max_dR_pct),
                      names_to = "q", values_to = "v") |>
  dplyr::mutate(q = ifelse(q == "median_max_dK_pct", "K", "per-cell R"))
p6 <- ggplot(d3, aes(q, v, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_text(aes(label = sprintf("%.1f%%", v)),
            position = position_dodge(width = 0.7), vjust = -0.35, size = 3.3) +
  scale_fill_manual(values = c(M0 = "grey45", M2g = "#C0392B"), name = NULL) +
  labs(x = NULL, y = "median max |change| across the window sweep (%)",
       title = "The key test: window sensitivity does not fall under the extension",
       subtitle = "if M2g captured something real, these bars would drop") +
  th + theme(legend.position = "top")
sv(p6, "fig6_sweep", 8, 4.2)

message("D6 figures done.")
