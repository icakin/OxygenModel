# =============================================================================
# d2_partA_reconcile.R - PART A: the three statements of cell density
# =============================================================================
# Puts three independent claims about how many cells were in a vial onto one
# scale (cells per litre):
#
#   (i)   STATED   - the paper's OD600 = 0.0005 at aliquoting, and the OD600
#                    actually recorded per replicate (data/OD_r_FC_r.csv)
#   (ii)  MEASURED - the per-replicate flow-cytometry count FC_Initial, with the
#                    full dilution chain from the methods reconstructed
#   (iii) USED     - N_inoculation_cells_per_L in data/Ninoc.csv, as 05 uses it
#
# Also: hunts the provenance of N_inoc, confirms/refutes that Cell_Counts.csv
# FC_Initial is the rounded taxon mean of the per-replicate values, and
# recomputes r_FC with the per-replicate initials (PART A2).
#
# Read-only over data/ and results/. Writes to runs/D2_analysis/.
# Run:  Rscript reports/D2_n0_routes/scripts/d2_partA_reconcile.R
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(lme4); library(lmerTest) })

BASE <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", "..", ".."),
  mustWork = FALSE)
if (!dir.exists(file.path(BASE, "scripts"))) BASE <- getwd()
OUT <- file.path(BASE, "runs", "D2_analysis")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

rd <- function(p) suppressWarnings(suppressMessages(
  readr::read_csv(p, show_col_types = FALSE, progress = FALSE)))

# =============================================================================
# CONVERSION CONSTANTS - every one of these is an assumption, stated explicitly
# =============================================================================

# ---- (a) the flow-cytometry dilution chain, from the paper's methods ---------
#   step 1  fixation : 490 ul culture + 10 ul 25% glutaraldehyde -> 500 ul
#   step 2  PBS      :  50 ul of that + 450 ul PBS               -> 500 ul
#   step 3  stain    : 500 ul + 2.5 ul of 100x SYBR Green I      -> 502.5 ul
# The CytoFLEX S reports events per microlitre OF THE MEASURED SAMPLE, so the
# original culture density is the measured density times the product of these.
FC_DIL_FIX   <- 500 / 490      # 1.020408
FC_DIL_PBS   <- 500 / 50       # 10
FC_DIL_STAIN <- 502.5 / 500    # 1.005
FC_DILUTION  <- FC_DIL_FIX * FC_DIL_PBS * FC_DIL_STAIN   # 10.2551
UL_PER_L     <- 1e6

# ---- (b) OD600 -> cells, generic ---------------------------------------------
# No taxon-specific OD-to-cell factor is published for these 15 environmental
# isolates (see the report's QUESTIONS FOR THE AUTHORS). The generic reference
# is the E. coli value used by the standard bio-calculators:
#     OD600 = 1.0  <->  8e8 cells/mL     (Agilent / Sigma-Aldrich OD600 note)
# with a documented spread across strain, growth phase, cell size and instrument
# optics of roughly 5e8 to 2e9 cells/mL. A range, not a point estimate.
OD_CELLS_PER_ML_MID <- 8e8
OD_CELLS_PER_ML_LO  <- 5e8
OD_CELLS_PER_ML_HI  <- 2e9
ML_PER_L <- 1e3

# ---- (c) an optional size adjustment, clearly labelled as approximate --------
# Turbidity per cell rises with cell size, so a taxon with cells larger than
# E. coli gives fewer cells per OD unit. A first-order (LINEAR-IN-BIOVOLUME)
# adjustment uses this repo's own cell volumes. This is deliberately crude:
# Koch's light-scattering treatment gives a scaling steeper than linear in
# volume for cells of this size, so this is reported as a sensitivity only, not
# as a preferred conversion.
ECOLI_WIDTH_UM <- 0.5; ECOLI_LENGTH_UM <- 2.0
rod_volume <- function(w, l) {
  r <- w / 2; cyl <- pmax(l - w, 0)
  pi * r^2 * cyl + (4/3) * pi * r^3
}
ECOLI_VOL_UM3 <- rod_volume(ECOLI_WIDTH_UM, ECOLI_LENGTH_UM)

PAPER_OD_STATED <- 0.0005

# =============================================================================
# Load
# =============================================================================
cc <- rd(file.path(BASE, "data", "Cell_Counts.csv")); names(cc)[1] <- "Taxon"
od <- rd(file.path(BASE, "data", "OD_r_FC_r.csv"));   names(od)[1] <- "Taxon"
ni <- rd(file.path(BASE, "data", "Ninoc.csv"))
cs <- rd(file.path(BASE, "data", "taxon_cell_sizes.csv"))
fit <- rd(file.path(BASE, "runs", "D1_baseline", "tables", "oxygen_results_with_R.csv"))

# =============================================================================
# A1. Is Cell_Counts.csv FC_Initial the rounded taxon mean of OD_r_FC_r.csv?
# =============================================================================
mean_check <- od %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(n_reps = dplyr::n(),
                   perrep_values = paste(FC_Initial, collapse = ", "),
                   perrep_mean = mean(FC_Initial), .groups = "drop") %>%
  dplyr::left_join(
    cc %>% dplyr::group_by(Taxon) %>%
      dplyr::summarise(n_distinct_in_CellCounts = dplyr::n_distinct(FC_Initial),
                       CellCounts_value = dplyr::first(FC_Initial), .groups = "drop"),
    by = "Taxon") %>%
  dplyr::mutate(rounded_mean = round(perrep_mean),
                constant_across_replicates = n_distinct_in_CellCounts == 1,
                equals_rounded_mean = CellCounts_value == rounded_mean)
readr::write_csv(mean_check, file.path(OUT, "A1_CellCounts_is_taxon_mean.csv"))

# Do the OD columns agree between the two files?
od_agree <- tibble::tibble(
  quantity = c("OD_Initial", "OD_Final", "FC_Initial", "FC_Final"),
  n_rows = nrow(cc),
  n_identical = c(sum(cc$OD_Initial == od$OD_Initial),
                  sum(cc$OD_Final   == od$OD_Final),
                  sum(cc$FC_Initial == od$FC_Initial),
                  sum(cc$FC_Final   == od$FC_Final))
) %>% dplyr::mutate(n_differing = n_rows - n_identical)
readr::write_csv(od_agree, file.path(OUT, "A1_CellCounts_vs_ODFC_columns.csv"))

# =============================================================================
# A2. The three-way reconciliation table
# =============================================================================
tri <- od %>%
  dplyr::transmute(Taxon = as.character(Taxon), Replicate = as.character(Replicate),
                   OD_Initial, OD_Final,
                   FC_Initial_perrep = FC_Initial, FC_Final_perrep = FC_Final,
                   Duration_min = Time) %>%
  dplyr::left_join(cc %>% dplyr::transmute(Taxon = as.character(Taxon),
                                           Replicate = as.character(Replicate),
                                           FC_Initial_taxonmean = FC_Initial),
                   by = c("Taxon", "Replicate")) %>%
  dplyr::left_join(ni %>% dplyr::transmute(Taxon = as.character(Taxon),
                                           Replicate = as.character(Replicate),
                                           N_inoc_cells_per_L = N_inoculation_cells_per_L,
                                           delta_min = delta_Ninoc_to_N0_min),
                   by = c("Taxon", "Replicate")) %>%
  dplyr::left_join(cs %>% dplyr::transmute(Taxon = as.character(Taxon),
                                           cell_volume_um3, cell_carbon_fg),
                   by = "Taxon") %>%
  dplyr::mutate(
    # (i) STATED: the paper's nominal OD at aliquoting
    stated_OD_mid = PAPER_OD_STATED * OD_CELLS_PER_ML_MID * ML_PER_L,
    stated_OD_lo  = PAPER_OD_STATED * OD_CELLS_PER_ML_LO  * ML_PER_L,
    stated_OD_hi  = PAPER_OD_STATED * OD_CELLS_PER_ML_HI  * ML_PER_L,
    # (i-b) the OD actually RECORDED for this replicate
    measured_OD_mid = OD_Initial * OD_CELLS_PER_ML_MID * ML_PER_L,
    measured_OD_lo  = OD_Initial * OD_CELLS_PER_ML_LO  * ML_PER_L,
    measured_OD_hi  = OD_Initial * OD_CELLS_PER_ML_HI  * ML_PER_L,
    # size-adjusted variant (sensitivity only)
    size_factor       = ECOLI_VOL_UM3 / cell_volume_um3,
    measured_OD_sizeadj = measured_OD_mid * size_factor,
    # (ii) MEASURED: flow cytometry, full dilution chain
    FC_cells_per_L      = FC_Initial_perrep   * FC_DILUTION * UL_PER_L,
    FC_cells_per_L_taxonmean = FC_Initial_taxonmean * FC_DILUTION * UL_PER_L,
    FC_cells_per_L_nodil = FC_Initial_perrep * UL_PER_L,   # the naive x1e6
    FC_end_cells_per_L  = FC_Final_perrep     * FC_DILUTION * UL_PER_L,
    # (iii) USED
    N_inoc = N_inoc_cells_per_L,
    # pairwise ratios
    ratio_FC_over_Ninoc        = FC_cells_per_L / N_inoc,
    ratio_measuredOD_over_Ninoc = measured_OD_mid / N_inoc,
    ratio_statedOD_over_Ninoc   = stated_OD_mid / N_inoc,
    ratio_FC_over_measuredOD    = FC_cells_per_L / measured_OD_mid,
    ratio_FC_over_statedOD      = FC_cells_per_L / stated_OD_mid,
    ratio_measuredOD_over_statedOD = OD_Initial / PAPER_OD_STATED
  )
readr::write_csv(tri, file.path(OUT, "A2_three_way_density_per_replicate.csv"))

tri_taxon <- tri %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(
    n = dplyr::n(),
    OD_init_med = stats::median(OD_Initial),
    stated_cells_per_L = dplyr::first(stated_OD_mid),
    measuredOD_cells_per_L_med = stats::median(measured_OD_mid),
    FC_cells_per_L_med = stats::median(FC_cells_per_L),
    N_inoc_med = stats::median(N_inoc),
    r_FCoverNinoc = stats::median(ratio_FC_over_Ninoc),
    r_measODoverNinoc = stats::median(ratio_measuredOD_over_Ninoc),
    r_FCovermeasOD = stats::median(ratio_FC_over_measuredOD),
    r_measODoverStatedOD = stats::median(ratio_measuredOD_over_statedOD),
    .groups = "drop")
readr::write_csv(tri_taxon, file.path(OUT, "A2_three_way_density_by_taxon.csv"))

# overall spread of each pairwise ratio
ratio_summary <- tri %>%
  dplyr::select(dplyr::starts_with("ratio_")) %>%
  tidyr::pivot_longer(dplyr::everything(), names_to = "ratio", values_to = "v") %>%
  dplyr::group_by(ratio) %>%
  dplyr::summarise(n = sum(is.finite(v)),
                   min = min(v, na.rm = TRUE), q25 = stats::quantile(v, .25, na.rm = TRUE),
                   median = stats::median(v, na.rm = TRUE),
                   q75 = stats::quantile(v, .75, na.rm = TRUE), max = max(v, na.rm = TRUE),
                   fold_spread = max(v, na.rm = TRUE) / min(v, na.rm = TRUE),
                   .groups = "drop")
readr::write_csv(ratio_summary, file.path(OUT, "A2_pairwise_ratio_summary.csv"))

# the assumptions table, so every conversion is auditable
assumptions <- tibble::tribble(
  ~quantity, ~value, ~source,
  "FC fixation dilution (490 ul + 10 ul glutaraldehyde)", sprintf("%.6f", FC_DIL_FIX), "paper, Materials and methods",
  "FC PBS dilution (50 ul into 450 ul PBS)", sprintf("%.1f", FC_DIL_PBS), "paper, Materials and methods",
  "FC stain dilution (2.5 ul SYBR into 500 ul)", sprintf("%.4f", FC_DIL_STAIN), "paper, Materials and methods",
  "FC total dilution factor", sprintf("%.4f", FC_DILUTION), "product of the three above",
  "CytoFLEX S output unit", "events per uL of the MEASURED (diluted) sample", "ASSUMED - instrument default, not stated in the paper",
  "OD600 = 1 -> cells/mL, central", format(OD_CELLS_PER_ML_MID, scientific = TRUE), "E. coli reference used by the standard OD600 bio-calculators (Agilent; Sigma-Aldrich OD600 application note)",
  "OD600 = 1 -> cells/mL, low", format(OD_CELLS_PER_ML_LO, scientific = TRUE), "documented spread across strain / phase / optics",
  "OD600 = 1 -> cells/mL, high", format(OD_CELLS_PER_ML_HI, scientific = TRUE), "documented spread across strain / phase / optics",
  "Taxon-specific OD-to-cell factor", "NONE AVAILABLE", "no published factor exists for these 15 environmental isolates - see QUESTIONS FOR THE AUTHORS",
  "E. coli reference cell volume", sprintf("%.4f um^3", ECOLI_VOL_UM3), "0.5 x 2.0 um rod, using this repo's own prolate-capsule formula (config.R)",
  "Size adjustment", "linear in biovolume", "APPROXIMATE sensitivity only; true OD-per-cell scaling is steeper than linear (Koch light-scattering)",
  "Paper's stated OD at aliquoting", sprintf("%.4f", PAPER_OD_STATED), "paper, Materials and methods"
)
readr::write_csv(assumptions, file.path(OUT, "A2_conversion_assumptions.csv"))

# =============================================================================
# A3. Where does N_inoc come from? Test every candidate.
# =============================================================================
# For each candidate predictor X, compute the implied factor N_inoc / X per row.
# If N_inoc were built from X by a constant factor, that ratio would be constant
# (CV = 0). The CV of the implied factor is therefore the test statistic.
cand <- tri %>%
  dplyr::transmute(
    Taxon, Replicate, N_inoc,
    `OD_Initial`                    = OD_Initial,
    `OD_Final`                      = OD_Final,
    `FC_Initial (per replicate)`    = FC_Initial_perrep,
    `FC_Initial (taxon mean)`       = FC_Initial_taxonmean,
    `FC_Final (per replicate)`      = FC_Final_perrep,
    `OD_Initial / cell_volume_um3`  = OD_Initial / cell_volume_um3,
    `OD_Initial / cell_carbon_fg`   = OD_Initial / cell_carbon_fg,
    `OD_Initial * cell_volume_um3`  = OD_Initial * cell_volume_um3,
    `FC_Initial / cell_volume_um3`  = FC_Initial_perrep / cell_volume_um3,
    `OD_Initial * FC_Initial`       = OD_Initial * FC_Initial_perrep,
    `constant (1)`                  = 1
  )

prov <- purrr::map_dfr(setdiff(names(cand), c("Taxon", "Replicate", "N_inoc")), function(v) {
  x <- cand[[v]]
  f <- cand$N_inoc / x
  ok <- is.finite(f) & f > 0
  lf <- log(cand$N_inoc[ok]); lx <- log(x[ok])
  r <- if (stats::sd(lx) > 0) stats::cor(lf, lx) else NA_real_
  tibble::tibble(
    candidate = v, n = sum(ok),
    implied_factor_median = stats::median(f[ok]),
    implied_factor_min = min(f[ok]), implied_factor_max = max(f[ok]),
    implied_factor_CV = stats::sd(f[ok]) / mean(f[ok]),
    implied_factor_fold_spread = max(f[ok]) / min(f[ok]),
    pearson_r_log_log = r
  )
}) %>% dplyr::arrange(implied_factor_CV)
readr::write_csv(prov, file.path(OUT, "A3_Ninoc_provenance_candidates.csv"))

# Is N_inoc per-replicate or per-taxon?
prov_struct <- tri %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(n = dplyr::n(), n_distinct_Ninoc = dplyr::n_distinct(N_inoc),
                   Ninoc_min = min(N_inoc), Ninoc_max = max(N_inoc),
                   Ninoc_fold_within_taxon = max(N_inoc) / min(N_inoc), .groups = "drop")
readr::write_csv(prov_struct, file.path(OUT, "A3_Ninoc_structure_by_taxon.csv"))

# Within-taxon: does N_inoc track OD_Initial or FC_Initial at all?
prov_within <- tri %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(
    n = dplyr::n(),
    cor_logNinoc_logOD = suppressWarnings(stats::cor(log(N_inoc), log(OD_Initial))),
    cor_logNinoc_logFC = suppressWarnings(stats::cor(log(N_inoc), log(FC_Initial_perrep))),
    .groups = "drop")
readr::write_csv(prov_within, file.path(OUT, "A3_Ninoc_within_taxon_correlations.csv"))

# =============================================================================
# A4. r_FC with PER-REPLICATE initials vs the published taxon-mean initials
# =============================================================================
growth_base <- rd(file.path(BASE, "runs", "D1_baseline", "tables", "growth_rates_combined.csv"))

gr <- tri %>%
  dplyr::transmute(Taxon, Replicate, Duration_min,
                   r_OD600 = (log(OD_Final) - log(OD_Initial)) / Duration_min,
                   # AS PUBLISHED: 06_main_figures.R reads data/OD_r_FC_r.csv, whose
                   # FC_Initial IS per-replicate. Verified against D1's
                   # growth_rates_combined.csv to exactly 0.
                   r_FC_perrep     = (log(FC_Final_perrep) - log(FC_Initial_perrep)) / Duration_min,
                   # COUNTERFACTUAL: the taxon-mean initial in data/Cell_Counts.csv,
                   # a file no automated stage reads.
                   r_FC_taxonmean  = (log(FC_Final_perrep) - log(FC_Initial_taxonmean)) / Duration_min) %>%
  dplyr::left_join(growth_base %>% dplyr::select(Taxon, Replicate, r_O2), by = c("Taxon", "Replicate"))
readr::write_csv(gr, file.path(OUT, "A4_growth_rates_perrep_vs_published.csv"))

agreement <- function(d, xv, label) {
  dd <- d %>% dplyr::filter(is.finite(r_O2), is.finite(.data[[xv]]))
  f  <- suppressMessages(lme4::lmer(
    stats::as.formula(paste0("r_O2 ~ ", xv, " + (1 | Taxon)")), data = dd))
  fe <- lme4::fixef(f)
  ols <- stats::lm(stats::as.formula(paste0("r_O2 ~ ", xv)), data = dd)
  diffv <- dd$r_O2 - dd[[xv]]
  tibble::tibble(
    variant = label, n = nrow(dd),
    lmer_intercept = unname(fe[1]), lmer_slope = unname(fe[2]),
    R2_lmer_fitted = stats::cor(stats::fitted(f), dd$r_O2)^2,
    ols_intercept = unname(stats::coef(ols)[1]), ols_slope = unname(stats::coef(ols)[2]),
    R2_ols = summary(ols)$r.squared,
    BA_bias = mean(diffv), BA_sd = stats::sd(diffv),
    BA_loa_lo = mean(diffv) - 1.96 * stats::sd(diffv),
    BA_loa_hi = mean(diffv) + 1.96 * stats::sd(diffv),
    BA_pct_within = 100 * mean(diffv >= mean(diffv) - 1.96 * stats::sd(diffv) &
                               diffv <= mean(diffv) + 1.96 * stats::sd(diffv)))
}

a4 <- dplyr::bind_rows(
  agreement(gr, "r_OD600",        "r_O2 ~ r_OD600 (unchanged)"),
  agreement(gr, "r_FC_perrep",    "r_O2 ~ r_FC, per-replicate initial (AS PUBLISHED)"),
  agreement(gr, "r_FC_taxonmean", "r_O2 ~ r_FC, taxon-mean initial (COUNTERFACTUAL)")
)
readr::write_csv(a4, file.path(OUT, "A4_agreement_published_vs_perrep.csv"))

a4_taxon <- gr %>%
  dplyr::group_by(Taxon) %>%
  dplyr::summarise(n = dplyr::n(),
                   bias_as_published = mean(r_O2 - r_FC_perrep, na.rm = TRUE),
                   bias_taxonmean    = mean(r_O2 - r_FC_taxonmean, na.rm = TRUE),
                   .groups = "drop")
readr::write_csv(a4_taxon, file.path(OUT, "A4_bland_altman_by_taxon.csv"))

message("d2_partA: wrote ", length(list.files(OUT)), " files to ", OUT)
