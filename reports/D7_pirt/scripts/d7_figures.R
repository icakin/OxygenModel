# =============================================================================
# D7 figures - built ONLY from the artefact CSVs in runs/D7_analysis/
# =============================================================================
source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d7_common.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

FIGD <- file.path(D7_BASE, "reports", "D7_pirt", "figures")
dir.create(FIGD, showWarnings = FALSE, recursive = TRUE)
sv <- function(p, n, w, h) {
  ggplot2::ggsave(file.path(FIGD, paste0(n, ".png")), p, width = w, height = h,
                  dpi = 200, bg = "white")
  message("  wrote figures/", n, ".png")
}
th <- ggplot2::theme_classic(base_size = 11)

A2 <- readr::read_csv(d7f("A2_gate_per_curve.csv"), show_col_types = FALSE)
A3 <- readr::read_csv(d7f("A3_gate_pirt_fits.csv"), show_col_types = FALSE)
A4 <- readr::read_csv(d7f("A4_gate_pseudomonas_per_curve.csv"), show_col_types = FALSE)
A5 <- readr::read_csv(d7f("A5_gate_pseudomonas_fits.csv"), show_col_types = FALSE)
D2 <- readr::read_csv(d7f("D2_identifiability_vs_spread.csv"), show_col_types = FALSE)

# ---- Fig 1: the gate - the fit under each N0 treatment ----------------------
long <- dplyr::bind_rows(
  tibble::tibble(G = A2$G_C_fg_cell_h, R = A2$R_default,   trt = "N0 default (r-dependent)"),
  tibble::tibble(G = A2$G_C_fg_cell_h, R = A2$R_fixed_tax, trt = "N0 fixed per taxon"),
  tibble::tibble(G = A2$G_C_fg_cell_h, R = A2$R_fixed_all, trt = "N0 fixed globally")) |>
  dplyr::mutate(trt = factor(trt, levels = unique(trt)))
lab <- A3 |> dplyr::transmute(trt = factor(dataset, levels = levels(long$trt)),
                              lab = sprintf("m = %.0f\nslope = %+.2f", m_intercept, slope_inv_Y))
p1 <- ggplot(long, aes(G, R)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey55") +
  geom_point(alpha = 0.65, size = 1.6, colour = "#2980B9") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "#C0392B", linewidth = 0.8) +
  geom_text(data = lab, aes(x = Inf, y = Inf, label = lab), hjust = 1.08, vjust = 1.3,
            size = 3.2, inherit.aes = FALSE) +
  facet_wrap(~ trt, scales = "free_y") +
  expand_limits(x = 0) +
  labs(x = "per-cell growth G (fg C / cell / h)", y = "per-cell respiration R",
       title = "The gate: the Pirt slope changes sign with how N0 is bookkept",
       subtitle = "same K and r throughout; only the r-dependence of N0 differs. x-axis extended to G = 0") +
  th
sv(p1, "fig1_gate", 11, 4.2)

# ---- Fig 2: Pseudomonas, the same test --------------------------------------
l2 <- dplyr::bind_rows(
  tibble::tibble(mu = A4$mu_per_hr, q = A4$q_default, trt = "N0 default (exp(45r))"),
  tibble::tibble(mu = A4$mu_per_hr, q = A4$q_fixed,   trt = "N0 fixed (no r)")) |>
  dplyr::mutate(trt = factor(trt, levels = unique(trt)))
lab2 <- A5 |> dplyr::transmute(trt = factor(dataset, levels = levels(l2$trt)),
                               lab = sprintf("m = %.2f\nslope = %+.2f", m_intercept, slope_inv_Y))
p2 <- ggplot(l2, aes(mu, q)) +
  geom_point(aes(colour = A4$Temperature[match(mu, A4$mu_per_hr)]), size = 2.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "#C0392B", linewidth = 0.8) +
  geom_text(data = lab2, aes(x = Inf, y = Inf, label = lab), hjust = 1.08, vjust = 1.3,
            size = 3.2, inherit.aes = FALSE) +
  scale_colour_viridis_c(name = "T (°C)") +
  facet_wrap(~ trt, scales = "free_y") + expand_limits(x = 0) +
  labs(x = expression("specific growth "*mu*" (per h)"),
       y = "specific respiration q",
       title = "Pseudomonas across temperature: the same sign flip",
       subtitle = "11_temperature_cue.R sets its own INOC_DELAY_MIN = 45, so its N0 carries r too") +
  th + theme(legend.position = "right")
sv(p2, "fig2_gate_pseudomonas", 11, 4.2)

# ---- Fig 3: identifiability -------------------------------------------------
p3a <- ggplot(D2, aes(G_spread_ratio, abs(corr_intercept_slope))) +
  geom_hline(yintercept = 1, linetype = 3) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.4) +
  scale_x_log10() + coord_cartesian(ylim = c(0.5, 1.02)) +
  labs(x = "spread in G (max/min, log scale)",
       y = "|corr(intercept, slope)|",
       title = "A. The intercept and slope stay anti-correlated",
       subtitle = "widening the growth range barely helps") + th
p3b <- ggplot(D2, aes(G_spread_ratio, pct_m_negative)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.4) +
  scale_x_log10() + coord_cartesian(ylim = c(0, 50)) +
  labs(x = "spread in G (max/min, log scale)",
       y = "% of simulations returning m < 0",
       title = "B. Even with a TRUE positive maintenance",
       subtitle = "~30% of datasets report a negative one") + th
sv(p3a | p3b, "fig3_identifiability", 10, 4.2)

message("D7 figures done.")
