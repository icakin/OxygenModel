# =============================================================================
# d2_partC_temperature.R - PART C: the published temperature result, three ways
# =============================================================================
# Compares the three N0 routes already run into runs/D2_{forward,backward,none}:
#
#   forward   N0 = 550 * exp(r * 45) cells/uL           AS PUBLISHED
#   backward  N0 = 550 * exp(r * 45) * scale            scale transferred from the
#             bacterial dataset (Pseudomonas is in both); the temperature
#             experiment has NO endpoint biomass, so a true backward route
#             cannot be computed from it.
#   none      N0 = 550 cells/uL                         LOWER BOUND, not a candidate
#
# Reports all four optima per route plus the Sharpe-Schoolfield / Arrhenius
# parameters, and adds the Topt of the quantity Fig 7 panel A actually PLOTS
# (growth_C_fg_per_hr), which the published parameter table does not contain.
#
# Run:  Rscript reports/D2_n0_routes/scripts/d2_partC_temperature.R
# =============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(minpack.lm) })

BASE <- normalizePath(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", "..", ".."),
  mustWork = FALSE)
if (!dir.exists(file.path(BASE, "scripts"))) BASE <- getwd()
OUT <- file.path(BASE, "runs", "D2_analysis")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

rd <- function(p) suppressWarnings(suppressMessages(
  readr::read_csv(p, show_col_types = FALSE, progress = FALSE)))

ROUTES <- c(forward = "FORWARD (as published)",
            backward = "BACKWARD (scale transferred from the bacterial data)",
            none = "NONE (N0 = N_inoc; LOWER BOUND, not a candidate)")

# =============================================================================
# C1. The parameter tables the pipeline itself writes
# =============================================================================
tpc <- purrr::imap_dfr(ROUTES, function(lab, r) {
  p <- file.path(BASE, "runs", paste0("D2_", r), "tables",
                 "SharpeSchoolfield_Temperature_Params_NEWformula.csv")
  if (!file.exists(p)) return(tibble::tibble())
  rd(p) %>% dplyr::mutate(route = r, route_label = lab, .before = 1)
})
readr::write_csv(tpc, file.path(OUT, "C1_thermal_parameters_three_routes.csv"))

# =============================================================================
# C2. Per-temperature means under each route
# =============================================================================
per_T <- purrr::imap_dfr(ROUTES, function(lab, r) {
  p <- file.path(BASE, "runs", paste0("D2_", r), "tables",
                 "oxygen_model_results_good_only_NEWformula.csv")
  if (!file.exists(p)) return(tibble::tibble())
  rd(p) %>%
    dplyr::group_by(Temperature) %>%
    dplyr::summarise(n = dplyr::n(),
                     N0_cells_per_mL = mean(N0_cells_per_mL),
                     r_per_hour = mean(r_per_hour),
                     growth_C_fg_per_hr = mean(growth_C_fg_per_hr),
                     resp_C_fg_per_hr = mean(resp_C_fg_per_hr),
                     resp_to_growth_C = mean(resp_to_growth_C),
                     CUE = mean(carbon_use_efficiency), .groups = "drop") %>%
    dplyr::mutate(route = r, .before = 1)
})
readr::write_csv(per_T, file.path(OUT, "C2_per_temperature_three_routes.csv"))

# =============================================================================
# C3. Topt of the quantity Fig 7 panel A PLOTS (growth_C_fg_per_hr)
# =============================================================================
# The published parameter table reports Topt for `growth_r_per_hr`, i.e. the raw
# growth rate, which contains no N0 at all - so it cannot move between routes by
# construction. Fig 7 panel A plots growth in CARBON units, which does contain
# N0. Its Topt is computed here so the comparison is complete.
k_B <- 8.617e-5
safe_exp <- function(z) exp(pmin(700, z))
ln_SS_two <- function(T_C, lnR0, E, El, Tl, Eh, Th) {
  T <- T_C + 273.15; TlK <- Tl + 273.15; ThK <- Th + 273.15
  lnR0 - E/(k_B*T) - log1p(safe_exp(El/k_B*(1/T - 1/TlK))) - log1p(safe_exp(Eh/k_B*(1/ThK - 1/T)))
}
ln_SS_one <- function(T_C, lnR0, E, Eh, Th) {
  T <- T_C + 273.15; ThK <- Th + 273.15
  lnR0 - E/(k_B*T) - log1p(safe_exp(Eh/k_B*(1/ThK - 1/T)))
}
ln_boltz <- function(T_C, lnR0, E) { T <- T_C + 273.15; lnR0 - E/(k_B*T) }

topt_of <- function(df, yvar) {
  d <- df %>% dplyr::filter(is.finite(.data[[yvar]]), .data[[yvar]] > 0)
  if (nrow(d) < 4) return(tibble::tibble(best_model = NA_character_, Topt_C = NA_real_,
                                         E_eV = NA_real_))
  T_min <- min(d$Temperature); T_max <- max(d$Temperature)
  T_opt_obs <- d$Temperature[which.max(d[[yvar]])]
  f2 <- try(nlsLM(as.formula(paste0("log(", yvar, ") ~ ln_SS_two(Temperature, lnR0, E, El, Tl, Eh, Th)")),
    data = d, start = list(lnR0 = log(max(d[[yvar]])), E = .6, El = .4, Tl = T_min + 2, Eh = 1.5, Th = T_opt_obs + 3),
    lower = c(lnR0=-Inf,E=.1,El=.1,Tl=0,Eh=.1,Th=0), upper = c(lnR0=Inf,E=2.5,El=2.5,Tl=60,Eh=5,Th=60),
    control = nls.lm.control(maxiter = 600)), silent = TRUE)
  f1 <- try(nlsLM(as.formula(paste0("log(", yvar, ") ~ ln_SS_one(Temperature, lnR0, E, Eh, Th)")),
    data = d, start = list(lnR0 = log(max(d[[yvar]])), E = .6, Eh = 1.5, Th = T_opt_obs + 5),
    lower = c(lnR0=-Inf,E=.1,Eh=.5,Th=0), upper = c(lnR0=Inf,E=2.5,Eh=6,Th=60),
    control = nls.lm.control(maxiter = 500)), silent = TRUE)
  f0 <- try(nlsLM(as.formula(paste0("log(", yvar, ") ~ ln_boltz(Temperature, lnR0, E)")),
    data = d, start = list(lnR0 = log(max(d[[yvar]])), E = .6),
    lower = c(lnR0=-Inf,E=.1), upper = c(lnR0=Inf,E=2.5),
    control = nls.lm.control(maxiter = 400)), silent = TRUE)
  aic <- vapply(list(two=f2, one=f1, bol=f0),
                function(f) if (inherits(f, "try-error")) Inf else AIC(f), numeric(1))
  if (all(is.infinite(aic))) return(tibble::tibble(best_model = NA_character_, Topt_C = NA_real_, E_eV = NA_real_))
  best <- names(which.min(aic)); fit <- list(two=f2, one=f1, bol=f0)[[best]]
  g <- tibble::tibble(Temperature = seq(T_min, T_max, length.out = 300))
  g$Pred <- safe_exp(predict(fit, newdata = g))
  g <- g %>% dplyr::filter(is.finite(Pred))
  pars <- coef(summary(fit))
  tibble::tibble(best_model = best,
                 Topt_C = if (nrow(g) > 1) g$Temperature[which.max(g$Pred)] else NA_real_,
                 E_eV = if ("E" %in% rownames(pars)) pars["E", "Estimate"] else NA_real_)
}

c3 <- purrr::imap_dfr(ROUTES, function(lab, r) {
  p <- file.path(BASE, "runs", paste0("D2_", r), "tables",
                 "oxygen_model_results_good_only_NEWformula.csv")
  if (!file.exists(p)) return(tibble::tibble())
  d <- rd(p) %>% dplyr::filter(fit_ok, growth_C_fg_per_hr > 0, resp_C_fg_per_hr > 0)
  topt_of(d, "growth_C_fg_per_hr") %>%
    dplyr::mutate(route = r, response = "growth_C_fg_per_hr (Fig 7 panel A)", .before = 1)
})
readr::write_csv(c3, file.path(OUT, "C3_growth_carbon_Topt_three_routes.csv"))

# =============================================================================
# C4. The comparison table the report quotes
# =============================================================================
D1 <- c(growth = 33.97993, resp = 37.52508, cue = 32.44147)
opt <- tpc %>%
  dplyr::select(route, response, Topt_C) %>%
  tidyr::pivot_wider(names_from = response, values_from = Topt_C) %>%
  dplyr::left_join(c3 %>% dplyr::select(route, Topt_growthC = Topt_C), by = "route") %>%
  dplyr::transmute(
    route,
    Topt_growth_rate_C = growth_r_per_hr,
    Topt_growth_carbon_C = Topt_growthC,
    Topt_respiration_C = resp_C_fg_per_hr,
    Topt_CUE_C = carbon_use_efficiency,
    gap_CUE_minus_growth_C = carbon_use_efficiency - growth_r_per_hr,
    shift_vs_D1_growth = growth_r_per_hr - D1[["growth"]],
    shift_vs_D1_resp   = resp_C_fg_per_hr - D1[["resp"]],
    shift_vs_D1_CUE    = carbon_use_efficiency - D1[["cue"]],
    CUE_in_30_35_band  = carbon_use_efficiency >= 30 & carbon_use_efficiency <= 35)
readr::write_csv(opt, file.path(OUT, "C4_optima_three_routes.csv"))

# Does the paper's stated conclusion still hold under each route?
verdict <- purrr::imap_dfr(ROUTES, function(lab, r) {
  d <- per_T %>% dplyr::filter(route == r) %>% dplyr::arrange(Temperature)
  o <- opt %>% dplyr::filter(route == r)
  tibble::tibble(
    route = r, route_label = lab,
    CUE_at_lowest_T = d$CUE[1], CUE_at_peak = max(d$CUE),
    T_of_observed_CUE_peak = d$Temperature[which.max(d$CUE)],
    CUE_at_highest_T = d$CUE[nrow(d)],
    increases_with_T = d$CUE[which.max(d$CUE)] > d$CUE[1],
    declines_at_highest_T = d$CUE[nrow(d)] < max(d$CUE),
    fitted_CUE_Topt = o$Topt_CUE_C,
    peaks_near_growth_optimum_30_35 = o$Topt_CUE_C >= 30 & o$Topt_CUE_C <= 35,
    conclusion_holds = (d$CUE[which.max(d$CUE)] > d$CUE[1]) &&
      (d$CUE[nrow(d)] < max(d$CUE)) &&
      (o$Topt_CUE_C >= 30 && o$Topt_CUE_C <= 35))
})
readr::write_csv(verdict, file.path(OUT, "C4_conclusion_holds_three_routes.csv"))

message("d2_partC: wrote C1-C4 artefacts to ", OUT)
