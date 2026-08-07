# =============================================================================
# d1_compare.R - the three D1 comparisons, as machine-readable artefacts
# =============================================================================
# Reads three trees and writes comparison CSVs that the D1 Quarto report renders
# without recomputing anything:
#
#   COMMITTED : results/tables            (outputs of the ORIGINAL scripts, in git)
#   MODULAR   : runs/D1_baseline/tables   (this run of the restructured pipeline)
#   ORIGINAL  : runs/D1_originals/Tables  (this run of scripts/original_scripts/)
#
#   D1  MODULAR vs COMMITTED   - does the restructure reproduce what is in git?
#   D2  MODULAR vs ORIGINAL    - does the restructure preserve the published numbers?
#   D3  headline numbers       - the claims the paper actually makes.
#
# Writes to runs/D1_baseline/comparisons/. Touches nothing in results/ or data/.
# Run:  Rscript reports/D1_baseline/scripts/d1_compare.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(lme4); library(lmerTest)
})

BASE <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", "..", ".."),
  mustWork = FALSE)
if (!dir.exists(file.path(BASE, "scripts"))) BASE <- getwd()

COMMITTED <- file.path(BASE, "results", "tables")
MODULAR   <- file.path(BASE, "runs", "D1_baseline", "tables")
ORIGINAL  <- file.path(BASE, "runs", "D1_originals", "Tables")
OUT       <- file.path(BASE, "runs", "D1_baseline", "comparisons")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

IDENTICAL_TOL <- 1e-10
EQUIV_TOL     <- 1e-6

# Curated by the 03 app; INPUTS to any run, never regenerated. Excluded from D1.
APP_INPUT_TABLES <- c("manual_fit_windows.csv", "plot_exclude_points.csv")

rd <- function(p) suppressWarnings(suppressMessages(
  readr::read_csv(p, show_col_types = FALSE, progress = FALSE)))

# Candidate key columns, most specific first. Used to align rows across trees.
KEY_CANDIDATES <- c("Taxon", "Temperature", "Replicate", "curve_code", "series_id",
                    "Method", "response", "scenario", "draw_model", "N0_CV",
                    "noise_sd", "dt", "duration", "Time", "Time_min", "Time0_min",
                    "r_true", "R_true", "rep_id")

pick_keys <- function(a, b) {
  k <- intersect(KEY_CANDIDATES, intersect(names(a), names(b)))
  if (!length(k)) return(character(0))
  # keep only the smallest prefix of keys that makes rows unique in BOTH tables
  for (i in seq_along(k)) {
    kk <- k[seq_len(i)]
    if (!anyDuplicated(a[kk]) && !anyDuplicated(b[kk])) return(kk)
  }
  # No key set identifies rows uniquely (e.g. the synthetic-recovery grid, whose
  # replicate index is not exported). Fall back to row-order alignment rather
  # than a many-to-many join that would fabricate rows.
  character(0)
}

rel_diff <- function(x, y) {
  d <- abs(x - y)
  den <- pmax(abs(x), abs(y))
  out <- ifelse(den > 0, d / den, ifelse(d > 0, Inf, 0))
  out
}

classify <- function(max_rel, max_abs, ok_structure) {
  if (!ok_structure) return("DIVERGENT")
  # Nothing numeric to compare (all-character table, or an empty table) and the
  # structure/character check already passed: that is an exact match.
  if (is.na(max_rel) && is.na(max_abs)) return("IDENTICAL")
  if (!is.finite(max_rel) && !is.finite(max_abs)) return("DIVERGENT")
  m <- if (is.finite(max_rel)) max_rel else Inf
  if (m < IDENTICAL_TOL || (is.finite(max_abs) && max_abs == 0)) "IDENTICAL"
  else if (m < EQUIV_TOL) "NUMERICALLY EQUIVALENT"
  else "DIVERGENT"
}

# ---------------------------------------------------------------------------
# Compare one CSV pair -> (table row, per-column rows)
# ---------------------------------------------------------------------------
compare_csv <- function(name, path_ref, path_new, label_ref, label_new) {
  if (!file.exists(path_new)) {
    return(list(
      tbl = tibble(table = name, class = "NOT REGENERATED",
                   n_rows_ref = if (file.exists(path_ref)) nrow(rd(path_ref)) else NA_integer_,
                   n_rows_new = NA_integer_, n_cols_compared = 0L,
                   max_abs_diff = NA_real_, max_rel_diff = NA_real_,
                   worst_column = NA_character_,
                   note = paste0("no counterpart in ", label_new)),
      cols = tibble()))
  }
  if (!file.exists(path_ref)) {
    return(list(
      tbl = tibble(table = name, class = "NEW OUTPUT",
                   n_rows_ref = NA_integer_, n_rows_new = nrow(rd(path_new)),
                   n_cols_compared = 0L, max_abs_diff = NA_real_,
                   max_rel_diff = NA_real_, worst_column = NA_character_,
                   note = paste0("present in ", label_new, " only")),
      cols = tibble()))
  }

  a <- rd(path_ref); b <- rd(path_new)
  notes <- character(0)
  shared <- intersect(names(a), names(b))
  only_ref <- setdiff(names(a), names(b)); only_new <- setdiff(names(b), names(a))
  if (length(only_ref)) notes <- c(notes, paste0("cols only in ", label_ref, ": ",
                                                 paste(only_ref, collapse = "/")))
  if (length(only_new)) notes <- c(notes, paste0("cols only in ", label_new, ": ",
                                                 paste(only_new, collapse = "/")))
  if (nrow(a) != nrow(b)) notes <- c(notes, sprintf("row count %d vs %d", nrow(a), nrow(b)))

  keys <- pick_keys(a[shared], b[shared])
  if (length(keys)) {
    m <- dplyr::inner_join(a[shared], b[shared], by = keys, suffix = c(".__ref", ".__new"))
    n_matched <- nrow(m)
    if (n_matched < max(nrow(a), nrow(b)))
      notes <- c(notes, sprintf("%d of %d/%d rows matched on %s",
                                n_matched, nrow(a), nrow(b), paste(keys, collapse = "+")))
    value_cols <- setdiff(shared, keys)
  } else if (nrow(a) == nrow(b)) {
    m <- dplyr::bind_cols(
      rlang::set_names(a[shared], paste0(shared, ".__ref")),
      rlang::set_names(b[shared], paste0(shared, ".__new")))
    n_matched <- nrow(m)
    value_cols <- shared
    notes <- c(notes, "aligned by row order (no key columns)")
  } else {
    return(list(
      tbl = tibble(table = name, class = "DIVERGENT",
                   n_rows_ref = nrow(a), n_rows_new = nrow(b), n_cols_compared = 0L,
                   max_abs_diff = NA_real_, max_rel_diff = NA_real_,
                   worst_column = NA_character_,
                   note = paste(c(notes, "cannot align rows"), collapse = "; ")),
      cols = tibble()))
  }

  col_rows <- purrr::map_dfr(value_cols, function(cn) {
    x <- m[[paste0(cn, ".__ref")]]; y <- m[[paste0(cn, ".__new")]]
    if (is.numeric(x) && is.numeric(y)) {
      ok <- is.finite(x) & is.finite(y)
      na_mismatch <- sum(is.finite(x) != is.finite(y))
      ad <- if (any(ok)) max(abs(x[ok] - y[ok])) else NA_real_
      rdif <- if (any(ok)) max(rel_diff(x[ok], y[ok])) else NA_real_
      tibble(column = cn, type = "numeric", n = sum(ok),
             max_abs_diff = ad, max_rel_diff = rdif,
             n_mismatch = sum(ok & x != y) + na_mismatch)
    } else {
      xs <- as.character(x); ys <- as.character(y)
      nm <- sum(!(xs == ys | (is.na(xs) & is.na(ys))), na.rm = TRUE)
      tibble(column = cn, type = "non-numeric", n = length(xs),
             max_abs_diff = NA_real_, max_rel_diff = NA_real_, n_mismatch = nm)
    }
  })

  num <- col_rows %>% dplyr::filter(type == "numeric")
  max_abs <- suppressWarnings(max(num$max_abs_diff, na.rm = TRUE))
  max_rel <- suppressWarnings(max(num$max_rel_diff, na.rm = TRUE))
  if (!is.finite(max_abs)) max_abs <- NA_real_
  if (!is.finite(max_rel)) max_rel <- if (nrow(num)) Inf else NA_real_
  worst <- if (nrow(num) && any(is.finite(num$max_rel_diff)))
    num$column[which.max(replace(num$max_rel_diff, !is.finite(num$max_rel_diff), -1))] else NA_character_

  chr_mismatch <- sum(col_rows$n_mismatch[col_rows$type == "non-numeric"], na.rm = TRUE)
  ok_structure <- length(only_ref) == 0 && length(only_new) == 0 &&
    nrow(a) == nrow(b) && n_matched == nrow(a) && chr_mismatch == 0
  if (chr_mismatch > 0) notes <- c(notes, paste0(chr_mismatch, " non-numeric cell mismatches"))

  cls <- classify(max_rel, max_abs, ok_structure)
  # A table that matches numerically but differs in shape is reported as such.
  if (cls %in% c("IDENTICAL", "NUMERICALLY EQUIVALENT") && !ok_structure)
    cls <- paste0(cls, " (SHAPE DIFFERS)")

  list(
    tbl = tibble(table = name, class = cls, n_rows_ref = nrow(a), n_rows_new = nrow(b),
                 n_cols_compared = nrow(col_rows),
                 max_abs_diff = max_abs, max_rel_diff = max_rel, worst_column = worst,
                 note = if (length(notes)) paste(notes, collapse = "; ") else ""),
    cols = col_rows %>% dplyr::mutate(table = name, .before = 1))
}

compare_txt <- function(name, path_ref, path_new, label_new) {
  if (!file.exists(path_new))
    return(tibble(file = name, class = "NOT REGENERATED", n_lines_ref = NA_integer_,
                  n_lines_new = NA_integer_, n_lines_differing = NA_integer_,
                  note = paste0("no counterpart in ", label_new)))
  if (!file.exists(path_ref))
    return(tibble(file = name, class = "NEW OUTPUT", n_lines_ref = NA_integer_,
                  n_lines_new = length(readLines(path_new, warn = FALSE)),
                  n_lines_differing = NA_integer_, note = ""))
  a <- readLines(path_ref, warn = FALSE); b <- readLines(path_new, warn = FALSE)
  n <- max(length(a), length(b))
  a2 <- c(a, rep("", n - length(a))); b2 <- c(b, rep("", n - length(b)))
  d <- sum(trimws(a2) != trimws(b2))

  # Line-by-line text diffing is defeated by print-width line wraps, which shift
  # every subsequent line. Compare the NUMERIC CONTENT independently: pull every
  # number out of each file, in order, and check they agree.
  nums <- function(v) {
    m <- regmatches(v, gregexpr("-?[0-9]+\\.?[0-9]*(e[-+]?[0-9]+)?", v, perl = TRUE))
    suppressWarnings(as.numeric(unlist(m)))
  }
  na <- nums(a); nb <- nums(b)
  num_note <- if (length(na) != length(nb)) {
    sprintf("numeric tokens: %d vs %d", length(na), length(nb))
  } else {
    ok <- is.finite(na) & is.finite(nb)
    mx <- if (any(ok)) max(rel_diff(na[ok], nb[ok])) else 0
    sprintf("%d numeric tokens, max relative difference %.3g", length(na), mx)
  }
  numbers_match <- length(na) == length(nb) &&
    { ok <- is.finite(na) & is.finite(nb)
      !any(ok) || max(rel_diff(na[ok], nb[ok])) < EQUIV_TOL }

  cls <- if (d == 0) "IDENTICAL"
         else if (numbers_match) "NUMERICALLY EQUIVALENT (formatting differs)"
         else "DIVERGENT"
  tibble(file = name, class = cls,
         n_lines_ref = length(a), n_lines_new = length(b), n_lines_differing = d,
         note = num_note)
}

# ---------------------------------------------------------------------------
# D1 : MODULAR vs COMMITTED
# ---------------------------------------------------------------------------
message("== D1: modular vs committed ==")
committed_csv <- sort(setdiff(basename(Sys.glob(file.path(COMMITTED, "*.csv"))), APP_INPUT_TABLES))
committed_txt <- sort(basename(Sys.glob(file.path(COMMITTED, "*.txt"))))

d1 <- purrr::map(committed_csv, ~compare_csv(.x, file.path(COMMITTED, .x),
                                             file.path(MODULAR, .x), "committed", "modular"))
d1_tables <- purrr::map_dfr(d1, "tbl")
d1_cols   <- purrr::map_dfr(d1, "cols")
d1_txt    <- purrr::map_dfr(committed_txt, ~compare_txt(.x, file.path(COMMITTED, .x),
                                                        file.path(MODULAR, .x), "modular"))
# tables the modular run produces that are not in results/
d1_new <- setdiff(basename(Sys.glob(file.path(MODULAR, "*.csv"))),
                  basename(Sys.glob(file.path(COMMITTED, "*.csv"))))

readr::write_csv(d1_tables, file.path(OUT, "D1_modular_vs_committed_tables.csv"))
readr::write_csv(d1_cols,   file.path(OUT, "D1_modular_vs_committed_columns.csv"))
readr::write_csv(d1_txt,    file.path(OUT, "D1_modular_vs_committed_textfiles.csv"))

# figures: names only (TIFF/PDF will never be byte-identical)
fig_committed <- sort(basename(Sys.glob(file.path(BASE, "results", "figures", "*"))))
fig_committed <- fig_committed[fig_committed != ".gitkeep"]
fig_modular   <- sort(basename(Sys.glob(file.path(BASE, "runs", "D1_baseline", "figures", "*"))))
fig_original  <- sort(basename(Sys.glob(file.path(BASE, "runs", "D1_originals", "plots", "*"))))
d1_figs <- tibble(figure = sort(unique(c(fig_committed, fig_modular)))) %>%
  dplyr::mutate(
    in_committed = figure %in% fig_committed,
    in_modular   = figure %in% fig_modular,
    in_original_run = figure %in% fig_original,
    status = dplyr::case_when(
      in_committed &  in_modular ~ "regenerated",
      in_committed & !in_modular &  in_original_run ~ "NOT regenerated by modular (original run makes it)",
      in_committed & !in_modular ~ "NOT regenerated by modular",
      TRUE ~ "new in modular run"))
readr::write_csv(d1_figs, file.path(OUT, "D1_figure_inventory.csv"))

# ---------------------------------------------------------------------------
# D2 : MODULAR vs ORIGINAL (both re-run now, same inputs)
# ---------------------------------------------------------------------------
message("== D2: modular vs original ==")
orig_csv <- sort(basename(Sys.glob(file.path(ORIGINAL, "*.csv"))))
orig_txt <- sort(basename(Sys.glob(file.path(ORIGINAL, "*.txt"))))
d2 <- purrr::map(orig_csv, ~compare_csv(.x, file.path(ORIGINAL, .x),
                                        file.path(MODULAR, .x), "original", "modular"))
d2_tables <- purrr::map_dfr(d2, "tbl")
d2_cols   <- purrr::map_dfr(d2, "cols")
d2_txt    <- purrr::map_dfr(orig_txt, ~compare_txt(.x, file.path(ORIGINAL, .x),
                                                   file.path(MODULAR, .x), "modular"))
readr::write_csv(d2_tables, file.path(OUT, "D2_modular_vs_original_tables.csv"))
readr::write_csv(d2_cols,   file.path(OUT, "D2_modular_vs_original_columns.csv"))
readr::write_csv(d2_txt,    file.path(OUT, "D2_modular_vs_original_textfiles.csv"))

# also: original vs committed, i.e. do the originals still reproduce what is in git?
message("== D2b: original vs committed ==")
d2b <- purrr::map(intersect(committed_csv, orig_csv),
                  ~compare_csv(.x, file.path(COMMITTED, .x), file.path(ORIGINAL, .x),
                               "committed", "original-rerun"))
readr::write_csv(purrr::map_dfr(d2b, "tbl"), file.path(OUT, "D2b_original_vs_committed_tables.csv"))

# ---- D2 per-QUANTITY view (not just per file) -------------------------------
QUANTITIES <- tribble(
  ~quantity,                        ~file,                                  ~column,
  "growth rate r (min^-1)",         "oxygen_fit_results.csv",               "r_per_minute",
  "scaling K",                      "oxygen_fit_results.csv",               "K",
  "fitted intercept O2_0",          "oxygen_fit_results.csv",               "O2_0",
  "reference O2 (O2_ref)",          "oxygen_fit_results.csv",               "O2_ref",
  "fit window length T_end (min)",  "oxygen_fit_results.csv",               "T_end_min",
  "total O2 consumed C_tot",        "oxygen_fit_results.csv",               "C_tot_mg_per_L",
  "pseudo R^2",                     "oxygen_fit_results.csv",               "pseudo_R2",
  "reconstructed N0 (cells/L)",     "oxygen_results_with_R.csv",            "N0_cells_per_L",
  "per-cell respiration R",         "oxygen_results_with_R.csv",            "R",
  "respiration R (fg C/cell/h)",    "oxygen_results_with_R.csv",            "R_C_fg_cell_h",
  "growth G (fg C/cell/h)",         "oxygen_results_with_R.csv",            "G_C_fg_cell_h",
  "r vs OD600 growth rate",         "growth_rates_combined.csv",            "r_OD600",
  "r vs flow cytometry",            "growth_rates_combined.csv",            "r_FC",
  "r from O2",                      "growth_rates_combined.csv",            "r_O2",
  "cutoff sensitivity r (main)",    "oxygen_model_results_main.csv",        "r_per_minute",
  "cutoff sensitivity r (O2>=0.5)", "oxygen_model_results_O2_ge_0.5.csv",   "r_per_minute",
  "MC relative SD of R (median)",   "N0_MC_ourmodel_overall_summary_ALL.csv","R_rel_sd_med",
  "CUE at each temperature",        "oxygen_model_results_good_only_NEWformula.csv", "carbon_use_efficiency",
  "growth C flux vs temperature",   "oxygen_model_results_good_only_NEWformula.csv", "growth_C_fg_per_hr",
  "respiration C flux vs temp",     "oxygen_model_results_good_only_NEWformula.csv", "resp_C_fg_per_hr",
  "thermal Topt (deg C)",           "SharpeSchoolfield_Temperature_Params_NEWformula.csv", "Topt_C",
  "activation energy E (eV)",       "SharpeSchoolfield_Temperature_Params_NEWformula.csv", "E_eV",
  "deactivation energy Eh (eV)",    "SharpeSchoolfield_Temperature_Params_NEWformula.csv", "Eh_eV",
  "synthetic recovery r bias",      "Table_S2_synthetic_parameter_recovery_summary.csv", "r_bias",
  "synthetic recovery R bias",      "Table_S2_synthetic_parameter_recovery_summary.csv", "R_bias"
)

quantity_row <- function(cols_tbl, q) {
  hit <- cols_tbl %>% dplyr::filter(table == q$file, column == q$column)
  if (!nrow(hit)) return(tibble(quantity = q$quantity, file = q$file, column = q$column,
                                n = NA_integer_, max_abs_diff = NA_real_,
                                max_rel_diff = NA_real_, verdict = "NOT COMPARABLE"))
  tibble(quantity = q$quantity, file = q$file, column = q$column, n = hit$n[1],
         max_abs_diff = hit$max_abs_diff[1], max_rel_diff = hit$max_rel_diff[1],
         verdict = classify(hit$max_rel_diff[1], hit$max_abs_diff[1], TRUE))
}
d2_quant <- purrr::map_dfr(seq_len(nrow(QUANTITIES)), ~quantity_row(d2_cols, QUANTITIES[.x, ]))
d1_quant <- purrr::map_dfr(seq_len(nrow(QUANTITIES)), ~quantity_row(d1_cols, QUANTITIES[.x, ]))
readr::write_csv(d2_quant, file.path(OUT, "D2_per_quantity.csv"))
readr::write_csv(d1_quant, file.path(OUT, "D1_per_quantity.csv"))

# ---------------------------------------------------------------------------
# D3 : the paper's headline numbers
# ---------------------------------------------------------------------------
message("== D3: headline numbers ==")

# (a) growth-rate agreement vs OD600 and FC, exactly as 06_main_figures.R does it
agreement_for <- function(growth_csv, tree) {
  if (!file.exists(growth_csv)) return(tibble())
  g <- rd(growth_csv)
  out <- purrr::map_dfr(c("r_OD600", "r_FC"), function(xv) {
    d <- g %>% dplyr::filter(is.finite(r_O2), is.finite(.data[[xv]]))
    if (nrow(d) < 5) return(tibble())
    f <- suppressMessages(lme4::lmer(
      stats::as.formula(paste0("r_O2 ~ ", xv, " + (1 | Taxon)")), data = d))
    fe <- lme4::fixef(f)
    # marginal (fixed-effects-only) R2 alongside the script's fitted-vs-observed R2
    r2_fitted <- stats::cor(stats::fitted(f), d$r_O2, use = "complete.obs")^2
    r2_ols    <- summary(stats::lm(stats::as.formula(paste0("r_O2 ~ ", xv)), data = d))$r.squared
    ols       <- stats::coef(stats::lm(stats::as.formula(paste0("r_O2 ~ ", xv)), data = d))
    tibble(tree = tree, comparison = paste0("r_O2 ~ ", xv), n = nrow(d),
           lmer_intercept = unname(fe[1]), lmer_slope = unname(fe[2]),
           R2_lmer_fitted = r2_fitted,
           ols_intercept = unname(ols[1]), ols_slope = unname(ols[2]), R2_ols = r2_ols)
  })
  out
}

# (b) Bland-Altman bias per method, exactly as Fig 5 computes it
ba_for <- function(growth_csv, tree) {
  if (!file.exists(growth_csv)) return(tibble())
  g <- rd(growth_csv)
  purrr::map_dfr(c("r_OD600", "r_FC"), function(xv) {
    d <- g %>% dplyr::transmute(Taxon, m1 = r_O2, m2 = .data[[xv]]) %>%
      dplyr::mutate(avg = (m1 + m2) / 2, diff = m1 - m2) %>%
      dplyr::filter(is.finite(avg), is.finite(diff))
    if (!nrow(d)) return(tibble())
    bias <- mean(d$diff); sdd <- stats::sd(d$diff)
    tibble(tree = tree, comparison = paste0("O2 - ", sub("^r_", "", xv)), n = nrow(d),
           bias = bias, sd_diff = sdd,
           loa_lower = bias - 1.96 * sdd, loa_upper = bias + 1.96 * sdd,
           pct_within_loa = 100 * mean(d$diff >= bias - 1.96 * sdd &
                                       d$diff <= bias + 1.96 * sdd),
           bias_rel_to_mean_r = bias / mean(d$m1))
  })
}

# (c) per-taxon Bland-Altman bias - "no taxon-specific bias" is a claim about this
ba_by_taxon <- function(growth_csv, tree) {
  if (!file.exists(growth_csv)) return(tibble())
  rd(growth_csv) %>%
    dplyr::select(Taxon, r_O2, r_OD600, r_FC) %>%
    tidyr::pivot_longer(c(r_OD600, r_FC), names_to = "method", values_to = "r_other") %>%
    dplyr::filter(is.finite(r_O2), is.finite(r_other)) %>%
    dplyr::group_by(tree = tree, method, Taxon) %>%
    dplyr::summarise(n = dplyr::n(), bias = mean(r_O2 - r_other),
                     .groups = "drop")
}

trees <- c(committed = COMMITTED, modular = MODULAR, original = ORIGINAL)
d3_agreement <- purrr::imap_dfr(trees, ~agreement_for(file.path(.x, "growth_rates_combined.csv"), .y))
d3_ba        <- purrr::imap_dfr(trees, ~ba_for(file.path(.x, "growth_rates_combined.csv"), .y))
d3_ba_taxon  <- purrr::imap_dfr(trees, ~ba_by_taxon(file.path(.x, "growth_rates_combined.csv"), .y))
readr::write_csv(d3_agreement, file.path(OUT, "D3_growth_rate_agreement.csv"))
readr::write_csv(d3_ba,        file.path(OUT, "D3_bland_altman.csv"))
readr::write_csv(d3_ba_taxon,  file.path(OUT, "D3_bland_altman_by_taxon.csv"))

# (d) Pseudomonas thermal parameters
d3_tpc <- purrr::imap_dfr(trees, function(dir, nm) {
  p <- file.path(dir, "SharpeSchoolfield_Temperature_Params_NEWformula.csv")
  if (!file.exists(p)) return(tibble())
  rd(p) %>% dplyr::mutate(tree = nm, .before = 1)
})
readr::write_csv(d3_tpc, file.path(OUT, "D3_thermal_parameters.csv"))

# (e) where the CUE optimum falls relative to the growth optimum
d3_opt <- d3_tpc %>%
  dplyr::select(tree, response, best_model, Topt_C) %>%
  tidyr::pivot_wider(names_from = response, values_from = c(best_model, Topt_C))
readr::write_csv(d3_opt, file.path(OUT, "D3_optima.csv"))

# observed (model-free) optima, straight off the per-temperature CUE table
d3_obs_opt <- purrr::imap_dfr(trees, function(dir, nm) {
  p <- file.path(dir, "oxygen_model_results_good_only_NEWformula.csv")
  if (!file.exists(p)) return(tibble())
  x <- rd(p)
  x %>% dplyr::group_by(Temperature) %>%
    dplyr::summarise(n = dplyr::n(),
                     mean_growth = mean(growth_C_fg_per_hr, na.rm = TRUE),
                     mean_resp   = mean(resp_C_fg_per_hr, na.rm = TRUE),
                     mean_r_per_hour = mean(r_per_hour, na.rm = TRUE),
                     mean_CUE    = mean(carbon_use_efficiency, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::mutate(tree = nm, .before = 1)
})
readr::write_csv(d3_obs_opt, file.path(OUT, "D3_cue_by_temperature.csv"))

message("d1_compare: wrote ", length(list.files(OUT)), " files to ", OUT)
