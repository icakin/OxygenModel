# =============================================================================
# d2_figures.R - the three figures the D2 report generates for itself
# =============================================================================
# Plots artefact CSVs from runs/D2_analysis/ into reports/D2_n0_routes/figures/.
# No model is fitted here; every number plotted was computed by
# d2_partA_reconcile.R, d2_partB_routes.R or d2_partC_temperature.R.
#
# Run:  Rscript reports/D2_n0_routes/scripts/d2_figures.R
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(scales); library(patchwork) })

BASE <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", "..", ".."),
  mustWork = FALSE)
if (!dir.exists(file.path(BASE, "scripts"))) BASE <- getwd()
CMP <- file.path(BASE, "runs", "D2_analysis")
FIG <- file.path(BASE, "reports", "D2_n0_routes", "figures")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

rd <- function(f) suppressWarnings(suppressMessages(
  readr::read_csv(file.path(CMP, f), show_col_types = FALSE, progress = FALSE)))

thm <- ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(axis.title = ggplot2::element_text(size = 11),
                 legend.position = "top", legend.title = ggplot2::element_blank(),
                 strip.background = ggplot2::element_blank())

# ---------------------------------------------------------------------------
# Fig 1. The three statements of cell density, per taxon
# ---------------------------------------------------------------------------
a2 <- rd("A2_three_way_density_per_replicate.csv")

f1_d <- a2 %>%
  dplyr::select(Taxon, Replicate,
                `Stated OD600 = 0.0005` = stated_OD_mid,
                `Measured OD600` = measured_OD_mid,
                `Flow cytometry` = FC_cells_per_L,
                `N_inoc (used by the pipeline)` = N_inoc) %>%
  tidyr::pivot_longer(-c(Taxon, Replicate), names_to = "source", values_to = "cells_per_L") %>%
  dplyr::mutate(source = factor(source, levels = c(
    "Stated OD600 = 0.0005", "Measured OD600", "Flow cytometry",
    "N_inoc (used by the pipeline)")),
    Taxon = factor(Taxon, levels = sort(unique(Taxon), decreasing = TRUE)))

p1 <- ggplot2::ggplot(f1_d, ggplot2::aes(cells_per_L, Taxon, colour = source, shape = source)) +
  ggplot2::geom_point(size = 2, alpha = 0.85,
                      position = ggplot2::position_nudge(y = 0)) +
  ggplot2::scale_x_log10(labels = scales::label_log(),
                         breaks = 10^(6:11), limits = c(3e5, 3e11)) +
  ggplot2::scale_colour_manual(values = c("#999999", "#56B4E9", "#009E73", "#D55E00")) +
  ggplot2::scale_shape_manual(values = c(3, 16, 17, 15)) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
  ggplot2::labs(x = expression("Cell density at inoculation (cells L"^-1*", log scale)"), y = NULL) +
  thm
ggplot2::ggsave(file.path(FIG, "fig1_three_way_density.pdf"), p1, width = 8, height = 5.2)
ggplot2::ggsave(file.path(FIG, "fig1_three_way_density.png"), p1, width = 8, height = 5.2, dpi = 300)

# ---------------------------------------------------------------------------
# Fig 2. The exponential-growth check, and the N0 ratio against r
# ---------------------------------------------------------------------------
b1 <- rd("B1_exponential_growth_check.csv")
b2 <- rd("B2_N0_both_routes_per_replicate.csv")

p2a <- ggplot2::ggplot(b1, ggplot2::aes(fold_predicted, fold_observed)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  ggplot2::geom_point(ggplot2::aes(colour = !keep_backward), size = 2, alpha = 0.85) +
  ggplot2::scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "#D55E00"),
                               labels = c("kept", "excluded")) +
  ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
  ggplot2::labs(x = expression("Predicted fold change, exp("*italic(r)*" × duration)"),
                y = "Observed flow-cytometry fold change",
                title = "(A) The shared assumption") + thm

p2b <- ggplot2::ggplot(b2 %>% dplyr::filter(keep_backward),
                       ggplot2::aes(r_per_minute, ratio_back_over_forward)) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  ggplot2::geom_point(size = 2, alpha = 0.85, colour = "#0072B2") +
  ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                       colour = "black", linewidth = 0.6) +
  ggplot2::scale_y_log10() +
  ggplot2::labs(x = expression("Fitted growth rate "*italic(r)*" (min"^-1*")"),
                y = expression(N[0]^backward / N[0]^forward),
                title = "(B) The disagreement grows with r") + thm

p2 <- p2a + p2b
ggplot2::ggsave(file.path(FIG, "fig2_growth_check_and_ratio.pdf"), p2, width = 9.5, height = 4.2)
ggplot2::ggsave(file.path(FIG, "fig2_growth_check_and_ratio.png"), p2, width = 9.5, height = 4.2, dpi = 300)

# ---------------------------------------------------------------------------
# Fig 3. CUE against temperature under the three routes
# ---------------------------------------------------------------------------
c2 <- rd("C2_per_temperature_three_routes.csv") %>%
  dplyr::mutate(route = factor(route, levels = c("forward", "backward", "none"),
                               labels = c("forward (as published)",
                                          "backward",
                                          "none (lower bound)")))
c4 <- rd("C4_optima_three_routes.csv") %>%
  dplyr::mutate(route = factor(route, levels = c("forward", "backward", "none"),
                               labels = c("forward (as published)",
                                          "backward",
                                          "none (lower bound)")))

p3 <- ggplot2::ggplot(c2, ggplot2::aes(Temperature, CUE, colour = route, shape = route)) +
  ggplot2::geom_vline(data = c4, ggplot2::aes(xintercept = Topt_CUE_C, colour = route),
                      linetype = "dotted", linewidth = 0.5, show.legend = FALSE) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 2.4) +
  ggplot2::scale_colour_manual(values = c("#D55E00", "#0072B2", "#999999")) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::labs(x = "Temperature (°C)",
                y = "Carbon use efficiency (mean of 4 replicates)",
                caption = "Dotted lines: the fitted CUE optimum under each route (32.44, 32.37, 32.31 °C).") +
  thm
ggplot2::ggsave(file.path(FIG, "fig3_cue_three_routes.pdf"), p3, width = 7.5, height = 4.4)
ggplot2::ggsave(file.path(FIG, "fig3_cue_three_routes.png"), p3, width = 7.5, height = 4.4, dpi = 300)

message("d2_figures: wrote 3 figures to ", FIG)
