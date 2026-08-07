# =============================================================================
# config.R - Shared paths, constants, model, inputs, and helper functions
# =============================================================================
# Sourced by every non-app script in the pipeline. (The Shiny apps 03/04 are
# standalone and locate the project the same way.)
#
# Project: bacterial oxygen-respiration from dissolved-O2 time series.
# The grouping variable is TAXON (a bacterial strain/genus, e.g. Bacillus,
# Burkholderia, Arthrobacter, Yersinia). Each O2 series is identified by
# (Taxon, Replicate) - or, for the temperature-gradient script (10), by
# (Taxon, Temperature, Replicate).
#
# Organised like the companion Candida pipeline: central config, numbered
# scripts, interactive trim-selector + cell-size Shiny apps, and a manual
# fit-window / exclusion override system that the apps write and this config
# auto-loads. The pipeline produces FIGURES (plus the tables they are built
# from); the OD600 / flow-cytometry growth-rate comparison from the original
# OxygenModel.R is intentionally NOT included.
#
# Edit base_dir and the USER INPUT blocks below in this one place only.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ===== Project root and output directories ===================================
# base_dir is AUTO-DETECTED as the parent of the scripts/ folder. Set
# base_dir_manual only if detection fails.
base_dir_manual <- ""

.detect_base_dir <- function() {
  if (nzchar(base_dir_manual)) return(normalizePath(base_dir_manual, mustWork = FALSE))
  if (base::exists(".this_dir", inherits = TRUE)) {
    d <- base::get(".this_dir", inherits = TRUE)
    if (length(d) == 1 && !is.na(d) && nzchar(d) && dir.exists(d)) {
      return(dirname(normalizePath(d, mustWork = FALSE)))
    }
  }
  d <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  if ((is.null(d) || is.na(d) || !nzchar(d)) &&
      requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    d <- tryCatch(dirname(rstudioapi::getActiveDocumentContext()$path),
                  error = function(e) NA_character_)
  }
  if (is.null(d) || is.na(d) || !nzchar(d)) return(getwd())
  dirname(normalizePath(d, mustWork = FALSE))
}

base_dir    <- .detect_base_dir()
message("config.R: base_dir = ", base_dir)
data_dir    <- file.path(base_dir, "data")
tables_dir  <- file.path(base_dir, "results", "tables")
figures_dir <- file.path(base_dir, "results", "figures")
rds_dir     <- file.path(base_dir, "results", "rds")
for (d in c(data_dir, tables_dir, figures_dir, rds_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

tbl <- function(name) file.path(tables_dir, name)
fig <- function(name) file.path(figures_dir, name)

# ===== Input / intermediate / app files ======================================
OXYGEN_LONG_CSV <- file.path(data_dir, "Oxygen_Data_Long.csv")          # 01/02/05/07
CUE_CSV         <- file.path(data_dir, "Oxygen_Data_Filtered_CUE.csv")  # 10 (temperature)
# Optional per-curve inoculation table (05). If present it overrides the scalar
# N_inoculation_cells_per_L below, per (Taxon, Replicate). Either name is accepted.
NINOC_CSV       <- file.path(data_dir, "Ninoc.csv")
NINOC_CSV_LEGACY <- file.path(data_dir, "Ninoc_and_deltaTime_to_N0.csv")

LONG_CSV      <- file.path(tables_dir, "Oxygen_All_Long.csv")            # 01 -> 02
TRIMMED_CSV   <- file.path(tables_dir, "Oxygen_Data_Smoothed_Trimmed.csv") # 02
FILTERED_CSV  <- file.path(tables_dir, "Oxygen_Data_Filtered.csv")       # 02 -> 05/07/08
SKIPPED_CSV   <- file.path(tables_dir, "Skipped_Series_Log.csv")         # 02
TRIM_META_CSV <- file.path(tables_dir, "Oxygen_Trimmed_Series_Metadata.csv") # 02 -> 03/05
CURVE_KEY_CSV <- file.path(tables_dir, "Oxygen_Curve_Code_Key.csv")      # 02
FITCURVES_CSV <- file.path(tables_dir, "oxygen_fit_curves.csv")          # 05 -> 06
RESULTS_FIT_CSV   <- file.path(tables_dir, "oxygen_fit_results.csv")     # 05
RESULTS_FINAL_CSV <- file.path(tables_dir, "oxygen_results_with_R.csv")  # 05 -> 06/08

# Files written by the Shiny apps (03 trim-selector, 04 cell-sizes).
MANUAL_FIT_WINDOWS_CSV <- file.path(tables_dir, "manual_fit_windows.csv")
PLOT_EXCLUDE_CSV       <- file.path(tables_dir, "plot_exclude_points.csv")
# Per-taxon cell sizes are a user-provided INPUT (like Ninoc.csv), so it lives in data/.
TAXON_CELL_SIZES_CSV   <- file.path(data_dir, "taxon_cell_sizes.csv")

# ===== The normalised O2 respiration model ===================================
#   O2_norm(t) = O2_0 + (K / r) * (1 - exp(r * t))    (N0 fixed to 1 in the fit)
# Per-cell respiration R is reconstructed afterwards from N0 (see N0_METHOD below).
resp_model <- function(r, K, t, O2_0) {
  O2_0 + (K / r) * (1 - exp(r * t))
}

# ===== nlsLM fit bounds and QC thresholds ====================================
FIT_LOWER <- c(r = 1e-4, K = 1e-5, O2_0 = 0.8)
FIT_UPPER <- c(r = 0.1,  K = 0.5,  O2_0 = 1.2)

REL_SE_THRESHOLD <- 0.15
PVAL_THRESHOLD   <- 0.001
R2_THRESHOLD     <- 0.90
MAX_RESID_RANGE  <- 0.20
AIC_IMPROVEMENT  <- 10
MAPE_MAX         <- 0.15

# ===== Fit-window standardisation (used in 05_oxygen_fits.R) =================
# Ending each fit window at a per-curve point makes the window length vary and
# biases K. USE_DRAWDOWN_WINDOW = TRUE instead ends the window when O2 has fallen
# by a fixed FRACTION of its total drawdown, so every series is fit over the same
# slice of its curve.
#   FALSE = legacy behaviour (fit the whole trimmed/selected window).
#   TRUE  = fractional-drawdown end (standardised); FIT_DRAWDOWN_FRAC in (0,1).
USE_DRAWDOWN_WINDOW <- FALSE
FIT_DRAWDOWN_FRAC   <- 0.45

# ===== USER INPUT: inoculation density (Ninoc) - FALLBACK ONLY ===============
# Only used when N0_METHOD = "initial" (the default depletion-anchored method
# derives N0 from FC_Final instead - see below). Inoculation density (cells / L).
#
# TWO WAYS to supply it (05_oxygen_fits.R), for the "initial" fallback:
#   (a) Per-curve CSV: data/Ninoc.csv with columns
#         Taxon, Replicate, N_inoculation_cells_per_L, [delta_Ninoc_to_N0_min]
#       If delta_Ninoc_to_N0_min is present it overrides the trimming-derived
#       delay for that curve; otherwise the delay from 02 is used.
#   (b) The single scalar below, used for any curve not covered by the CSV.
N_inoculation_cells_per_L <- 1e9

# Loader for the optional per-curve Ninoc table (used by 05). Returns NULL if no
# file exists. Accepts either filename; maps the first 3-4 columns if the header
# names differ but the column order matches.
load_ninoc_table <- function() {
  p <- c(NINOC_CSV, NINOC_CSV_LEGACY)
  p <- p[file.exists(p)]
  if (!length(p)) return(NULL)
  tb <- tryCatch(readr::read_csv(p[1], show_col_types = FALSE), error = function(e) NULL)
  if (is.null(tb)) return(NULL)
  names(tb) <- trimws(names(tb))
  if (!all(c("Taxon", "Replicate", "N_inoculation_cells_per_L") %in% names(tb)) && ncol(tb) >= 3) {
    nm <- c("Taxon", "Replicate", "N_inoculation_cells_per_L", "delta_Ninoc_to_N0_min")
    names(tb)[seq_len(min(4, ncol(tb)))] <- nm[seq_len(min(4, ncol(tb)))]
  }
  if (!all(c("Taxon", "Replicate", "N_inoculation_cells_per_L") %in% names(tb))) return(NULL)
  tb$Taxon <- as.character(tb$Taxon); tb$Replicate <- as.character(tb$Replicate)
  tb$N_inoculation_cells_per_L <- as.numeric(tb$N_inoculation_cells_per_L)
  if (!"delta_Ninoc_to_N0_min" %in% names(tb)) tb$delta_Ninoc_to_N0_min <- NA_real_
  else tb$delta_Ninoc_to_N0_min <- as.numeric(tb$delta_Ninoc_to_N0_min)
  message("  Loaded Ninoc table from ", basename(p[1]), " (", nrow(tb), " rows).")
  tb[, c("Taxon", "Replicate", "N_inoculation_cells_per_L", "delta_Ninoc_to_N0_min")]
}

# N0 anchoring for the "initial" fallback ONLY (used when N0_METHOD = "initial";
# the default depletion-anchored method below does not use this flag).
#   TRUE  = N0 = N_inoc * exp(r * delta), delta = inoculation -> trim start.
#   FALSE = N0 = N_inoc (anchor at consumption onset).
N0_BACKPROJECT <- TRUE

# Minutes between inoculation and the FIRST O2 reading. Only used as a fallback
# additive offset when 02's per-curve delta is unavailable. Set to your protocol.
INOC_DELAY_MIN <- 0

# ===== Depletion-anchored N0 (default method) ================================
# N0 is derived from the FINAL flow-cytometry count (OD_r_FC_r.csv, FC_Final),
# projected back to the fit start over the interval during which growth actually
# occurred:
#     N0 = FC_Final * FC_TO_CELLS_PER_L * exp(-r * (t_depletion - fit_start))
# t_depletion = the time each vial reaches DEPLETION_FRAC of its starting O2, i.e.
# when O2 runs out and aerobic growth stops. This uses the true growth interval and
# removes the -r*delta geometric coupling of the initial-count back-projection
# (the artefact that made per-cell R anti-correlate with growth). r and the fit
# windows are unchanged; only N0 (and therefore R, CUE) changes.
# Set N0_METHOD <- "initial" to fall back to the old N_inoc * exp(r*delta) route.
N0_METHOD         <- "depletion"   # "depletion" (default) or "initial"
DEPLETION_FRAC    <- 0.10          # O2 fraction remaining that marks growth-stop (90% depleted)
FC_TO_CELLS_PER_L <- 909916        # FC events -> cells/L (dilution x sample-volume calibration);
                                   # sets ABSOLUTE scale only - slope and relative R are independent of it.
OD_FC_CSV         <- file.path(data_dir, "OD_r_FC_r.csv")

# Per-vial O2 depletion time from the full raw series (LONG_CSV). Returns
# Taxon, Replicate, t_depletion_min, O2_start, depleted (FALSE = threshold never
# reached within the recording; the recording end is then used as a lower bound).
load_depletion_table <- function() {
  if (!file.exists(LONG_CSV)) return(NULL)
  raw <- tryCatch(readr::read_csv(LONG_CSV, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(raw) || !all(c("Taxon", "Replicate", "Time", "Oxygen") %in% names(raw))) return(NULL)
  raw <- raw[order(raw$Taxon, raw$Replicate, raw$Time), ]
  out <- by(raw, list(raw$Taxon, raw$Replicate), function(g) {
    if (nrow(g) < 3) return(NULL)
    o0  <- mean(head(g$Oxygen, 3), na.rm = TRUE)
    hit <- which(g$Oxygen <= DEPLETION_FRAC * o0)
    reached <- length(hit) > 0
    data.frame(Taxon = as.character(g$Taxon[1]), Replicate = as.character(g$Replicate[1]),
               t_depletion_min = if (reached) g$Time[hit[1]] else max(g$Time, na.rm = TRUE),
               O2_start = o0, depleted = reached, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# Final flow-cytometry counts (events) per curve, from OD_r_FC_r.csv.
load_fc_final <- function() {
  if (!file.exists(OD_FC_CSV)) return(NULL)
  tb <- tryCatch(readr::read_csv(OD_FC_CSV, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(tb)) return(NULL)
  names(tb) <- trimws(names(tb))
  if (!all(c("Taxon", "Replicate", "FC_Final") %in% names(tb))) return(NULL)
  data.frame(Taxon = as.character(tb$Taxon), Replicate = as.character(tb$Replicate),
             FC_Final = as.numeric(tb$FC_Final), stringsAsFactors = FALSE)
}

# Depletion-anchored N0 (cells/L): back-project FC_Final to the fit start.
n0_depletion <- function(FC_Final, r_per_min, t_depletion_min, fit_start_min) {
  FC_Final * FC_TO_CELLS_PER_L * exp(-r_per_min * (t_depletion_min - fit_start_min))
}

# ===== USER INPUT: cell carbon ===============================================
# A bacterial cell modelled as a rod (prolate capsule): a cylinder of length
# (CELL_LENGTH_UM - CELL_WIDTH_UM) capped by two hemispheres of diameter
# CELL_WIDTH_UM, at CARBON_DENSITY_FG_PER_UM3 fg C / um^3. This sets carbon per
# cell, which scales growth (fg C / h). It does NOT change respiration or the
# temperature-response shapes.
#
# PER-TAXON sizes: run 04_cell_sizes.R to give each strain its own cell carbon.
# It writes tables/taxon_cell_sizes.csv, which 05_oxygen_fits.R uses per taxon;
# the values below are the GLOBAL fallback for any taxon not listed.
CELL_WIDTH_UM             <- 0.65
CELL_LENGTH_UM            <- 2.25
CARBON_DENSITY_FG_PER_UM3 <- 100
RESPIRATORY_QUOTIENT      <- 1

.cell_radius_um     <- CELL_WIDTH_UM / 2
.cell_cyl_length_um <- max(CELL_LENGTH_UM - CELL_WIDTH_UM, 0)
CELL_VOLUME_UM3 <-
  pi * (.cell_radius_um^2) * .cell_cyl_length_um +
  (4 / 3) * pi * (.cell_radius_um^3)
CELL_CARBON_FG_PER_CELL <- CELL_VOLUME_UM3 * CARBON_DENSITY_FG_PER_UM3

MG_TO_FG <- 1e12
MIN_TO_H <- 60
# O2 mass -> C mass conversion (mg C per mg O2, 1 mol O2 : 1 mol C respired).
M_C_G_PER_MOL  <- 12.011
M_O2_G_PER_MOL <- 31.998
O2_to_C_mass   <- (M_C_G_PER_MOL / M_O2_G_PER_MOL) * RESPIRATORY_QUOTIENT

# Per-taxon cell carbon: reads taxon_cell_sizes.csv (written by 04_cell_sizes.R)
# if present, else falls back to the global CELL_CARBON_FG_PER_CELL.
cell_carbon_of <- function(taxa) {
  taxa <- as.character(taxa)
  out  <- rep(CELL_CARBON_FG_PER_CELL, length(taxa))
  if (file.exists(TAXON_CELL_SIZES_CSV)) {
    cs <- tryCatch(readr::read_csv(TAXON_CELL_SIZES_CSV, show_col_types = FALSE),
                   error = function(e) NULL)
    if (!is.null(cs) && all(c("Taxon", "cell_carbon_fg") %in% names(cs))) {
      m <- stats::setNames(as.numeric(cs$cell_carbon_fg), as.character(cs$Taxon))
      hit <- m[taxa]
      out[is.finite(hit)] <- hit[is.finite(hit)]
    }
  }
  out
}

# ===== MANUAL fit-window overrides / excluded curves =========================
# Two tables keyed on (Taxon, Replicate):
#   PLOT_EXCLUDE_POINTS = curves discarded from the fits and every plot.
#   MANUAL_FIT_WINDOWS  = per-curve fit start/end (NA = automatic on that side).
# Start EMPTY: every curve is trimmed automatically by 02_trimming.R. Run
# 03_trim_selector.R to review curves and set these (it auto-saves the two CSVs
# below, which are loaded here when USE_APP_TRIM_FILES <- TRUE).
PLOT_EXCLUDE_POINTS <- data.frame(
  Taxon = character(0), Replicate = character(0), stringsAsFactors = FALSE
)
MANUAL_FIT_WINDOWS <- data.frame(
  Taxon = character(0), Replicate = character(0),
  fit_start = numeric(0), fit_end = numeric(0), stringsAsFactors = FALSE
)

USE_APP_TRIM_FILES <- TRUE
if (isTRUE(USE_APP_TRIM_FILES)) {
  if (file.exists(MANUAL_FIT_WINDOWS_CSV)) {
    .tmp <- tryCatch(read.csv(MANUAL_FIT_WINDOWS_CSV, stringsAsFactors = FALSE),
                     error = function(e) NULL)
    if (!is.null(.tmp) &&
        all(c("Taxon", "Replicate", "fit_start", "fit_end") %in% names(.tmp))) {
      MANUAL_FIT_WINDOWS <- .tmp
      message("  Loaded MANUAL_FIT_WINDOWS from app file (", nrow(.tmp), " curves).")
    }
  }
  if (file.exists(PLOT_EXCLUDE_CSV)) {
    .tmp <- tryCatch(read.csv(PLOT_EXCLUDE_CSV, stringsAsFactors = FALSE),
                     error = function(e) NULL)
    if (!is.null(.tmp) && all(c("Taxon", "Replicate") %in% names(.tmp))) {
      PLOT_EXCLUDE_POINTS <- .tmp
      message("  Loaded PLOT_EXCLUDE_POINTS from app file (", nrow(.tmp), " points).")
    }
  }
}

# ===== Data-scope filters ====================================================
# EXCLUDE_TAXA drops whole taxa from the analysis (05 onward). NULL = none.
EXCLUDE_TAXA <- NULL

drop_excluded_taxa <- function(df, what = "data") {
  if (is.null(EXCLUDE_TAXA) || !length(EXCLUDE_TAXA) || is.null(df) || nrow(df) == 0) return(df)
  if (!"Taxon" %in% names(df)) return(df)
  n0 <- nrow(df)
  df <- df[!(as.character(df$Taxon) %in% EXCLUDE_TAXA), , drop = FALSE]
  if (nrow(df) < n0)
    message(sprintf("EXCLUDE_TAXA: dropped %d of %d %s rows (%s).",
                    n0 - nrow(df), n0, what, paste(EXCLUDE_TAXA, collapse = ", ")))
  df
}

# Apply the app/exclusion list to a per-curve table (Taxon, Replicate keyed).
apply_excluded_curves <- function(df, what = "curves") {
  if (is.null(PLOT_EXCLUDE_POINTS) || nrow(PLOT_EXCLUDE_POINTS) == 0) return(df)
  if (!all(c("Taxon", "Replicate") %in% names(df))) return(df)
  ex <- paste(PLOT_EXCLUDE_POINTS$Taxon, toupper(PLOT_EXCLUDE_POINTS$Replicate))
  key <- paste(df$Taxon, toupper(df$Replicate))
  n0 <- nrow(df); df <- df[!(key %in% ex), , drop = FALSE]
  if (nrow(df) < n0)
    message(sprintf("PLOT_EXCLUDE_POINTS: dropped %d of %d %s.", n0 - nrow(df), n0, what))
  df
}

# ===== Taxon registry, palette and figure design constants ===================
TAXON_ORDER <- NULL   # e.g. c("Bacillus", "Burkholderia", "Arthrobacter", "Yersinia")

TAXON_PALETTE <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000",
  "#332288", "#88CCEE", "#117733", "#999933",
  "#CC6677", "#AA4499", "#661100", "#BBBBBB"
)

taxon_colors <- function(taxa) {
  taxa <- as.character(unique(taxa))
  if (!is.null(TAXON_ORDER)) taxa <- c(intersect(TAXON_ORDER, taxa), setdiff(taxa, TAXON_ORDER))
  pal <- rep(TAXON_PALETTE, length.out = length(taxa))
  stats::setNames(pal[seq_along(taxa)], taxa)
}

# The representative (Taxon, Replicate) series shown in the Fig 2 facet.
FIG2_SELECTED_COMBOS <- tibble::tribble(
  ~Taxon,          ~Replicate,
  "Bacillus",      "R4",
  "Burkholderia",  "R2",
  "Arthrobacter",  "R2",
  "Yersinia",      "R1"
)

# ===== Shared plot theme and figure-saving helpers ===========================
isme_theme <- function(base_size = 14) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = base_size + 2),
      axis.text  = ggplot2::element_text(size = base_size),
      axis.line  = ggplot2::element_line(linewidth = 0.8),
      axis.ticks = ggplot2::element_line(linewidth = 0.6),
      legend.position = "none"
    )
}

save_tiff <- function(plot, path, width, height, dpi = 600) {
  ggplot2::ggsave(filename = path, plot = plot, width = width, height = height,
                  dpi = dpi, device = "tiff", compression = "lzw")
  message("  wrote ", path)
}

save_png_pdf <- function(plot, path_noext, width, height, dpi = 600) {
  ggplot2::ggsave(paste0(path_noext, ".png"), plot, width = width, height = height,
                  dpi = dpi, bg = "white")
  tryCatch(
    ggplot2::ggsave(paste0(path_noext, ".pdf"), plot, width = width, height = height,
                    device = if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf,
                    bg = "white"),
    error = function(e) NULL
  )
  message("  wrote ", path_noext, ".png + .pdf")
}

message("config.R loaded: base_dir = ", base_dir)
message("  Project: bacterial O2 respiration (Taxon x Replicate) | figures-only pipeline")
message("  N_inoc = ", N_inoculation_cells_per_L, " cells/L | cell C = ",
        round(CELL_CARBON_FG_PER_CELL, 2), " fg | N0 method = ", N0_METHOD)
