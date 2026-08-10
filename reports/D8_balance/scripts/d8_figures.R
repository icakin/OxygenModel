# =============================================================================
# D8 figures - built ONLY from the artefact CSVs in runs/D8_analysis/
# =============================================================================
source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d8_common.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

FIGD <- file.path(D8_BASE, "reports", "D8_balance", "figures")
dir.create(FIGD, showWarnings = FALSE, recursive = TRUE)
sv <- function(p, n, w, h) {
  ggplot2::ggsave(file.path(FIGD, paste0(n, ".png")), p, width = w, height = h,
                  dpi = 200, bg = "white")
  message("  wrote figures/", n, ".png")
}
th <- ggplot2::theme_classic(base_size = 11)

A1 <- readr::read_csv(d8f("A1_pseudomonas_scale_free.csv"), show_col_types = FALSE)
A2 <- readr::read_csv(d8f("A2_activation_energies.csv"), show_col_types = FALSE)
A3 <- readr::read_csv(d8f("A3_argmin_and_foldchange.csv"), show_col_types = FALSE)
B2 <- readr::read_csv(d8f("B2_N0_expfactor_by_temperature.csv"), show_col_types = FALSE)
B4 <- readr::read_csv(d8f("B4_Topt_under_N0_treatments.csv"), show_col_types = FALSE)
B6 <- readr::read_csv(d8f("B6_Topt_gap_decomposition.csv"), show_col_types = FALSE)
C1 <- readr::read_csv(d8f("C1_massbalance_per_replicate.csv"), show_col_types = FALSE)
C5 <- readr::read_csv(d8f("C5_ratio_sensitivity_to_quota.csv"), show_col_types = FALSE)

# ---- Fig 1: E_r and E_K -----------------------------------------------------
e <- A2 |> dplyr::filter(grepl("^E_r|physical|E_K - E_r", quantity)) |>
  dplyr::mutate(quantity = factor(quantity, levels = rev(quantity)))
p1a <- ggplot(e, aes(E_eV, quantity)) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey55") +
  geom_errorbarh(aes(xmin = E_lo, xmax = E_hi), height = 0.12, linewidth = 0.8) +
  geom_point(size = 2.8, colour = "#C0392B") +
  labs(x = "activation energy (eV)", y = NULL,
       title = "A. Thermal sensitivity, scale-free",
       subtitle = "respiration is more sensitive than growth; the difference excludes 0") + th

arr <- A1 |> dplyr::group_by(Temperature) |>
  dplyr::summarise(y = stats::median(ratio_phys), .groups = "drop")
amin <- A3$value[A3$quantity == "argmin (K*O2_ref/r)  [solubility-corrected]"]
alo  <- A3$lo[A3$quantity == "argmin (K*O2_ref/r)  [solubility-corrected]"]
ahi  <- A3$hi[A3$quantity == "argmin (K*O2_ref/r)  [solubility-corrected]"]
Tpub <- A3$value[A3$quantity == "pipeline T_opt(CUE) from 11_temperature_cue.R"]
p1b <- ggplot(A1, aes(Temperature, ratio_phys)) +
  annotate("rect", xmin = alo, xmax = ahi, ymin = -Inf, ymax = Inf,
           fill = "#2980B9", alpha = 0.15) +
  geom_point(alpha = 0.55, size = 1.7) +
  geom_line(data = arr, aes(Temperature, y), linewidth = 0.9, colour = "grey30") +
  geom_vline(xintercept = amin, colour = "#2980B9", linewidth = 0.9) +
  geom_vline(xintercept = Tpub, colour = "#C0392B", linewidth = 0.9, linetype = 2) +
  annotate("text", x = amin, y = Inf, label = sprintf(" balance optimum %.1f°C", amin),
           colour = "#2980B9", hjust = 0, vjust = 1.6, size = 3.1) +
  annotate("text", x = Tpub, y = Inf, label = sprintf(" published T_opt(CUE) %.1f°C ", Tpub),
           colour = "#C0392B", hjust = 1, vjust = 3.2, size = 3.1) +
  scale_y_log10() +
  labs(x = "temperature (°C)", y = expression((K %.% O[2*ref])/r~"(log scale)"),
       title = "B. The balance ratio and its minimum",
       subtitle = "shaded = 95% bootstrap interval on the minimum") + th
sv(p1a | p1b, "fig1_scale_free", 11, 4.4)

# ---- Fig 2: the T_opt gap, decomposed ---------------------------------------
d <- B6 |> dplyr::slice(1:3) |>
  dplyr::mutate(step = factor(c("N0 temperature-independent\n(scale-free balance optimum)",
                                "+ N0 = exp(45r)\n(as the pipeline computes it)",
                                "+ Sharpe-Schoolfield fit\n(the published value)"),
                              levels = c("N0 temperature-independent\n(scale-free balance optimum)",
                                         "+ N0 = exp(45r)\n(as the pipeline computes it)",
                                         "+ Sharpe-Schoolfield fit\n(the published value)")))
p2 <- ggplot(d, aes(value_C, step)) +
  geom_segment(aes(x = min(d$value_C), xend = value_C, y = step, yend = step),
               colour = "grey80", linewidth = 1.2) +
  geom_point(size = 3.6, colour = "#C0392B") +
  geom_text(aes(label = sprintf("%.2f °C", value_C)), hjust = -0.25, size = 3.4) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(x = "T_opt(CUE) (°C)", y = NULL,
       title = "The 7.09 °C gap between the balance optimum and the published value",
       subtitle = "+2.51 °C from the 45-min back-projection, +4.58 °C from the functional form") +
  th
sv(p2, "fig2_topt_decomposition", 9.5, 3.6)

# ---- Fig 3: exp(r dt) by temperature, 45 vs 10 min --------------------------
b <- B2 |> tidyr::pivot_longer(c(median_expfac_45min, median_expfac_10min),
                               names_to = "dt", values_to = "f") |>
  dplyr::mutate(dt = ifelse(grepl("45", dt), "delta_t = 45 min (current)",
                            "delta_t = 10 min (pre-equilibrated)"))
p3 <- ggplot(b, aes(Temperature, f, colour = dt)) +
  geom_hline(yintercept = 1, linetype = 3) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
  scale_colour_manual(values = c("delta_t = 45 min (current)" = "#C0392B",
                                 "delta_t = 10 min (pre-equilibrated)" = "#2980B9"),
                      name = NULL) +
  labs(x = "temperature (°C)", y = expression("N"[0]*" inflation factor  exp(r"%.%Delta*"t)"),
       title = "The N₀ back-projection inflates biomass, and unevenly across temperature",
       subtitle = "it is the WARM/COLD DIFFERENTIAL, not the level, that contaminates a thermal comparison") +
  th + theme(legend.position = "top")
sv(p3, "fig3_N0_expfactor", 8.5, 4.4)

# ---- Fig 4: kinetic vs mass-balance CUE -------------------------------------
p4a <- ggplot(C1, aes(CUE_kinetic, CUE_massbalance)) +
  geom_abline(slope = 1, linetype = 2, colour = "grey45") +
  geom_point(aes(colour = frac_O2_post_window), size = 2.2) +
  scale_colour_viridis_c(name = "fraction of O₂\nconsumed after\nthe window") +
  coord_equal(xlim = c(0, 0.9), ylim = c(0, 0.9)) +
  labs(x = "kinetic CUE (fit window)", y = "mass-balance CUE (whole run)",
       title = "A. Two independent estimates",
       subtitle = "median ratio 1.077; no replicate outside [0,1]") + th
p4b <- ggplot(C5, aes(quota_scale, median_ratio)) +
  geom_hline(yintercept = 1, linetype = 3) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.6) +
  scale_x_log10(breaks = C5$quota_scale) +
  coord_cartesian(ylim = c(0.95, 1.2)) +
  labs(x = "carbon quota, relative to the value used", y = "median ratio (mass balance / kinetic)",
       title = "B. The ratio is only mildly quota-dependent",
       subtitle = "1.04–1.12 across a 16-fold range") + th
sv(p4a | p4b, "fig4_cue_comparison", 11, 4.6)

message("D8 figures done.")
