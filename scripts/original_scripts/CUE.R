################################################################################
# README – O₂-based growth–respiration analysis across a temperature gradient
#
# This script:
#   1) Fits the normalised O₂ growth–respiration model to each
#      Temperature × Replicate time series (N₀ fixed to 1 in the model).
#   2) Reconstructs N₀ at the start of O₂ measurements from a fixed
#      inoculation density (N_inoc) and the fitted growth rate r.
#   3) Converts growth and respiration to carbon units
#      (fg C mL⁻¹ h⁻¹ and fg C mL⁻¹ min⁻¹), and computes
#      respiration:growth ratios and carbon use efficiency (CUE).
#   4) Fits temperature–response curves (Sharpe–Schoolfield and Arrhenius),
#      selects the best model by AIC, and extracts Topt and activation
#      energies (E, Eh) with approximate 95% confidence intervals.
#
# Inputs (in data/):
#   - Oxygen_Data_Filtered_CUE.csv
#       Trimmed O₂ time series with columns:
#       Taxon, Temperature, Replicate, Time, Oxygen
#
# Tables (written to Tables/):
#   - oxygen_model_results_good_only_NEWformula.csv
#   - SharpeSchoolfield_Temperature_Params_NEWformula.csv
#
# Plots (written to plots/):
#   - oxygen_dynamics_all_models_NEWformula.pdf
#   - oxygen_dynamics_fullsize_per_page_NEWformula.pdf
#   - Fig_7_SharpeSchoolfield_Temperature_Fits_NEWformula.pdf
################################################################################

## ───────────────────────── 1.  Libraries ──────────────────────────────── ##
suppressPackageStartupMessages({
  library(tidyverse)
  library(minpack.lm)
  library(patchwork)
  library(zoo)
})

## ───────────────────────── 2.  QC / Fit thresholds ───────────────────── ##
REL_SE_THRESHOLD <- 0.15
PVAL_THRESHOLD   <- 0.001
R2_THRESHOLD     <- 0.90
MAX_RESID_RANGE  <- 0.20
AIC_IMPROVEMENT  <- 10
MAPE_MAX         <- 0.15

## ───────────────────────── 3.  Project paths & inoculation N_inoc ────── ##
data_dir   <- "data"
plots_dir  <- "plots"
tables_dir <- "Tables"

dir.create(data_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# Fixed inoculation density at t = 0 (before O2 measurements)
INOC_CELLS_PER_uL <- 550       # cells µL⁻¹
INOC_DELAY_MIN    <- 45        # minutes between inoculation and first O2 reading

## ───────────────────────── 4.  Cell-carbon constants ─────────────────── ##
# Rod geometry
cell_width  <- 0.65           # µm
cell_length <- 2.25           # µm
cell_radius <- cell_width / 2

cell_vol_um3 <- pi * cell_radius^2 * (cell_length - cell_width) +    # cylinder
  (4/3) * pi * cell_radius^3                                         # caps

C_density_fg_per_um3 <- 100
cell_C_fg            <- cell_vol_um3 * C_density_fg_per_um3          # fg C cell⁻¹

# Helper: mg O₂ L⁻¹ → mol O₂ mL⁻¹  (MW O₂ = 32 g)
mgL_to_mol_per_mL <- function(mg_per_L) (mg_per_L * 1e-3) / 32 / 1000

## ───────────────────────── 5.  Load data (from data/ folder) ─────────── ##
OXYGEN_CSV <- file.path(data_dir, "Oxygen_Data_Filtered_CUE.csv")

oxygen_data <- read_csv(OXYGEN_CSV, show_col_types = FALSE) %>%
  mutate(
    Taxon       = as.character(Taxon),
    Temperature = as.numeric(Temperature),
    Replicate   = as.character(Replicate)
  )

## ───────────────────────── 6.  Model function (N0 fixed to 1) ───────── ##
# Our model:
#   O2_norm(t) = O2_0 + (K / r) * (1 - exp(r * t))
resp_model <- function(r, K, t, O2_0) {
  O2_0 + (K / r) * (1 - exp(r * t))
}

## ───────────────────────── 7.  Plot theme ────────────────────────────── ##
isme_theme <- function() {
  theme_classic(base_size = 14) +
    theme(
      axis.title = element_text(size = 16),
      axis.text  = element_text(size = 14),
      axis.line  = element_line(linewidth = 0.8),
      axis.ticks = element_line(linewidth = 0.6),
      legend.position = "none"
    )
}

## ───────────────────────── 8.  Prepare outputs ───────────────────────── ##
results <- tibble(
  Taxon                 = character(),
  Temperature           = numeric(),
  Replicate             = character(),
  N0_cells_per_mL       = numeric(),   # N0 at start of O₂ (from N_inoc & r)
  r_per_minute          = numeric(),
  r_per_hour            = numeric(),
  K                     = numeric(),   # model parameter (1/min, total scaling)
  resp_rate             = numeric(),   # per-cell respiration (mg O₂ cell⁻¹ min⁻¹)
  O2_0                  = numeric(),
  O2_ref                = numeric(),   # mean initial O₂ (mg L⁻¹) used for normalisation
  AICc                  = numeric(),
  lnO2_change_per_min   = numeric(),
  pseudo_R2             = numeric(),
  fit_ok                = logical(),
  growth_C_fg_per_hr    = numeric(),
  growth_C_fg_per_min   = numeric(),
  resp_C_fg_per_hr      = numeric(),
  resp_C_fg_per_min     = numeric(),
  resp_to_growth_C      = numeric(),
  carbon_use_efficiency = numeric()
)

plots_list <- list()

## ───────────────────────── 9.  Fitting loop (new model) ─────────────── ##
grouped <- oxygen_data %>% dplyr::group_by(Taxon, Temperature, Replicate)
combos  <- dplyr::group_keys(grouped)

for (i in seq_len(nrow(combos))) {
  Tax  <- combos$Taxon[i]
  Temp <- combos$Temperature[i]
  Rep  <- combos$Replicate[i]
  
  df <- grouped %>%
    dplyr::filter(Taxon == Tax, Temperature == Temp, Replicate == Rep) %>%
    dplyr::arrange(Time)
  if (nrow(df) < 5) next
  
  # Derivative-based window detection
  df <- df %>%
    dplyr::mutate(
      dO2   = c(NA_real_, diff(Oxygen)),
      dt    = c(NA_real_, diff(Time)),
      dO2dt = dO2 / dt,
      sm    = zoo::rollmean(dO2dt, 3, fill = NA_real_, align = "right")
    )
  idx <- which(df$sm < -1e-7)[1]
  if (!is.na(idx)) idx <- idx + 15
  if (is.na(idx) || idx > nrow(df) - 2) idx <- 1
  
  df <- df[idx:nrow(df), ] %>%
    dplyr::mutate(Time0 = Time - min(Time, na.rm = TRUE))
  
  # Normalise O₂ for fitting
  O0 <- mean(head(df$Oxygen, 3), na.rm = TRUE)  # mg O₂ L⁻¹
  df <- df %>% dplyr::mutate(Oxygen_norm = Oxygen / O0)
  
  # Starting values for r and K
  r_start <- {
    seg    <- head(df, max(3, floor(0.3 * nrow(df))))
    slopes <- abs(diff(log(pmax(seg$Oxygen_norm, 1e-6))) / diff(seg$Time0))
    pmin(pmax(max(slopes, na.rm = TRUE), 1e-4), 5e-2)
  }
  K_start <- {
    slope <- suppressWarnings(min(diff(df$Oxygen_norm) / diff(df$Time0), na.rm = TRUE))
    if (!is.finite(slope)) slope <- -1e-3
    k_guess <- abs(slope)
    k_guess <- pmin(pmax(k_guess, 1e-5), 0.1)
    k_guess
  }
  
  fit <- tryCatch(
    nlsLM(
      Oxygen_norm ~ resp_model(r, K, Time0, O2_0),
      data    = df,
      start   = list(r = r_start, K = K_start, O2_0 = 1),
      lower   = c(r = 1e-4,  K = 1e-5, O2_0 = 0.8),
      upper   = c(r = 0.1,   K = 0.5,  O2_0 = 1.2),
      control = nls.lm.control(maxiter = 300)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) next
  
  # Outlier pruning and re-fit
  pred <- predict(fit, df)
  df_kept <- df[abs(df$Oxygen_norm - pred) < 2 * sd(df$Oxygen_norm - pred, na.rm = TRUE), ]
  fit <- tryCatch(update(fit, data = df_kept), error = function(e) fit)
  
  # Extract parameters
  pars  <- coef(summary(fit))
  r_est <- pars["r", "Estimate"]           # min⁻¹
  K_est <- pars["K", "Estimate"]           # 1/min
  
  pseudo_R2 <- 1 - sum(residuals(fit)^2) /
    sum((df_kept$Oxygen_norm - mean(df_kept$Oxygen_norm))^2)
  
  # QC
  fit_ok <- all(
    abs(pars[, "Std. Error"] / pars[, "Estimate"]) < REL_SE_THRESHOLD,
    pars[, "Pr(>|t|)"] < PVAL_THRESHOLD,
    pseudo_R2 >= R2_THRESHOLD,
    diff(range(residuals(fit), na.rm = TRUE)) < MAX_RESID_RANGE,
    is.finite(K_est),
    K_est > 0, K_est < 0.5,
    dplyr::between(r_est, 1e-4, 0.1),
    AIC(lm(Oxygen_norm ~ 1, data = df_kept)) - AIC(fit) >= AIC_IMPROVEMENT,
    mean(abs(residuals(fit) / df_kept$Oxygen_norm)) < MAPE_MAX
  )
  if (!fit_ok) next
  
  # Log-change (for diagnostics)
  lnchg <- (log(dplyr::last(df_kept$Oxygen_norm)) - log(dplyr::first(df_kept$Oxygen_norm))) /
    (dplyr::last(df_kept$Time0) - dplyr::first(df_kept$Time0))
  
  ## ── Reconstruct N0 at O₂ start from N_inoc and r ─────────────────────
  N0_O2start_cells_per_uL <- INOC_CELLS_PER_uL * exp(r_est * INOC_DELAY_MIN)
  N0_cells_per_mL         <- N0_O2start_cells_per_uL * 1e3
  N0_cells_per_L          <- N0_O2start_cells_per_uL * 1e6
  
  ## ── Carbon-unit conversions ───────────────────────────────────────────
  biomass_fgC_per_mL   <- N0_cells_per_mL * cell_C_fg
  growth_C_fg_per_hr   <- r_est * 60 * biomass_fgC_per_mL
  growth_C_fg_per_min  <- growth_C_fg_per_hr / 60
  
  O0_mol_per_mL            <- mgL_to_mol_per_mL(O0)
  resp_molO2_per_mL_per_hr <- K_est * O0_mol_per_mL * 60
  resp_C_fg_per_hr         <- resp_molO2_per_mL_per_hr * 12e15
  resp_C_fg_per_min        <- resp_C_fg_per_hr / 60
  
  resp_rate_per_cell <- if (is.finite(N0_cells_per_L) && N0_cells_per_L > 0) {
    K_est * O0 / N0_cells_per_L
  } else {
    NA_real_
  }
  
  carbon_use_efficiency <- growth_C_fg_per_hr /
    (growth_C_fg_per_hr + resp_C_fg_per_hr)
  
  resp_to_growth_C <- resp_C_fg_per_hr / growth_C_fg_per_hr
  
  results <- results %>%
    dplyr::add_row(
      Taxon                 = Tax,
      Temperature           = Temp,
      Replicate             = Rep,
      N0_cells_per_mL       = N0_cells_per_mL,
      r_per_minute          = r_est,
      r_per_hour            = r_est * 60,
      K                     = K_est,
      resp_rate             = resp_rate_per_cell,
      O2_0                  = coef(fit)["O2_0"],
      O2_ref                = O0,
      AICc                  = AIC(fit),
      lnO2_change_per_min   = lnchg,
      pseudo_R2             = pseudo_R2,
      fit_ok                = TRUE,
      growth_C_fg_per_hr    = growth_C_fg_per_hr,
      growth_C_fg_per_min   = growth_C_fg_per_min,
      resp_C_fg_per_hr      = resp_C_fg_per_hr,
      resp_C_fg_per_min     = resp_C_fg_per_min,
      resp_to_growth_C      = resp_to_growth_C,
      carbon_use_efficiency = carbon_use_efficiency
    )
  
  df_kept$Pred <- predict(fit, df_kept)
  label <- paste(Tax, Temp, Rep, sep = "_")
  plots_list[[label]] <-
    ggplot(df_kept, aes(Time0, Oxygen_norm)) +
    geom_point(size = 2.5) +
    geom_line(aes(y = Pred), colour = "red", linewidth = 0.8) +
    labs(
      title = label,
      x     = expression(Time ~ "(min)"),
      y     = expression(Normalised ~ O[2])
    ) +
    isme_theme()
}

## ───────────────────────── 10.  Save results & plots ─────────────────── ##
write_csv(
  results,
  file.path(tables_dir, "oxygen_model_results_good_only_NEWformula.csv")
)

if (length(plots_list) > 0) {
  pdf(
    file.path(plots_dir, "oxygen_dynamics_all_models_NEWformula.pdf"),
    width = 14, height = 10
  )
  print(patchwork::wrap_plots(plots_list))
  dev.off()
  
  pdf(
    file.path(plots_dir, "oxygen_dynamics_fullsize_per_page_NEWformula.pdf"),
    width = 8, height = 6
  )
  purrr::walk(plots_list, print)
  dev.off()
}

################################################################################
# 11.  Temperature-response fits (3 panels + Topt & activation energies)
#   • Two-sided Sharpe–Schoolfield, one-sided SS, and Arrhenius
#   • Lowest AIC wins (per response variable)
#   • Growth activation energy is based on r_per_hour (intrinsic growth)
#   • CUE is fitted on logit scale (bounded 0–1)
#
# IMPORTANT: ONLY CHANGE = tighter CUE bounds for Eh and Th so the curve
# must crash within the observed range (i.e. at ~40°C).
################################################################################

k_B <- 8.617e-5

safe_exp <- function(z) exp(pmin(700, z))
predict_safe <- function(fit, newdata) {
  out <- try(predict(fit, newdata = newdata), silent = TRUE)
  if (inherits(out, "try-error")) rep(NA_real_, nrow(newdata)) else out
}

clamp01  <- function(p, eps = 1e-4) pmin(pmax(p, eps), 1 - eps)
logit    <- function(p) log(p / (1 - p))
invlogit <- function(x) 1 / (1 + exp(-x))

ln_SS_two <- function(T_C, lnR0, E, El, Tl, Eh, Th) {
  T   <- T_C + 273.15
  TlK <- Tl  + 273.15
  ThK <- Th  + 273.15
  lnR0 - E /(k_B *   T) -
    log1p( safe_exp(El / k_B * (1/T    - 1/TlK) ) ) -
    log1p( safe_exp(Eh / k_B * (1/ThK  - 1/T   ) ) )
}
ln_SS_one <- function(T_C, lnR0, E, Eh, Th) {
  T   <- T_C + 273.15
  ThK <- Th  + 273.15
  lnR0 - E /(k_B * T) -
    log1p( safe_exp(Eh / k_B * (1/ThK - 1/T)) )
}
ln_boltz <- function(T_C, lnR0, E) {
  T <- T_C + 273.15
  lnR0 - E /(k_B * T)
}

eta_SS_two <- function(T_C, eta0, E, El, Tl, Eh, Th) {
  T   <- T_C + 273.15
  TlK <- Tl  + 273.15
  ThK <- Th  + 273.15
  eta0 - E /(k_B *   T) -
    log1p( safe_exp(El / k_B * (1/T    - 1/TlK) ) ) -
    log1p( safe_exp(Eh / k_B * (1/ThK  - 1/T   ) ) )
}
eta_SS_one <- function(T_C, eta0, E, Eh, Th) {
  T   <- T_C + 273.15
  ThK <- Th  + 273.15
  eta0 - E /(k_B * T) -
    log1p( safe_exp(Eh / k_B * (1/ThK - 1/T)) )
}
eta_boltz <- function(T_C, eta0, E) {
  T <- T_C + 273.15
  eta0 - E /(k_B * T)
}

fit_T_plot <- function(df, yvar, ylab, show_x = TRUE) {
  
  if (yvar == "carbon_use_efficiency") {
    
    df <- df %>%
      dplyr::filter(is.finite(.data[[yvar]])) %>%
      dplyr::mutate(
        cue_adj   = clamp01(.data[[yvar]]),
        cue_logit = logit(cue_adj)
      )
    
    base <- ggplot(df, aes(Temperature, cue_adj)) +
      geom_point(size = 2.5) +
      labs(
        x = if (show_x) expression(Temperature~"(°C)") else NULL,
        y = ylab
      ) +
      isme_theme()
    
    if (nrow(df) < 4) return(base)
    
    T_min     <- min(df$Temperature)
    T_max     <- max(df$Temperature)
    T_opt_obs <- df$Temperature[which.max(df$cue_adj)]
    
    # ------------------ CUE bounds: FORCE crash at Tmax ------------------
    Th_low  <- T_max - 3
    Th_high <- T_max + 0.2
    Eh_low  <- 6
    Eh_high <- 300
    # --------------------------------------------------------------------
    
    start_two <- list(
      eta0 = logit(max(df$cue_adj)),
      E  = 0.6,
      El = 0.4,  Tl = T_min + 2,
      Eh = 10,   Th = T_max - 1
    )
    fit_two <- try(
      nlsLM(
        cue_logit ~ eta_SS_two(Temperature, eta0, E, El, Tl, Eh, Th),
        data   = df,
        start  = start_two,
        lower  = c(eta0=-Inf, E=0.1, El=0.1, Tl=0,  Eh=Eh_low,  Th=Th_low),
        upper  = c(eta0= Inf, E=2.5, El=2.5, Tl=60, Eh=Eh_high, Th=Th_high),
        control = nls.lm.control(maxiter = 1200)
      ),
      silent = TRUE
    )
    
    start_one <- list(
      eta0 = logit(max(df$cue_adj)),
      E  = 0.6,
      Eh = 10,
      Th = T_max - 1
    )
    fit_one <- try(
      nlsLM(
        cue_logit ~ eta_SS_one(Temperature, eta0, E, Eh, Th),
        data   = df,
        start  = start_one,
        lower  = c(eta0=-Inf, E=0.1, Eh=Eh_low,  Th=Th_low),
        upper  = c(eta0= Inf, E=2.5, Eh=Eh_high, Th=Th_high),
        control = nls.lm.control(maxiter = 1200)
      ),
      silent = TRUE
    )
    
    fit_bol <- try(
      nlsLM(
        cue_logit ~ eta_boltz(Temperature, eta0, E),
        data  = df,
        start = list(eta0 = logit(max(df$cue_adj)), E = 0.6),
        lower = c(eta0=-Inf, E=0.1),
        upper = c(eta0= Inf, E=2.5),
        control = nls.lm.control(maxiter = 800)
      ),
      silent = TRUE
    )
    
    get_AIC <- function(f) if (inherits(f, "try-error")) Inf else AIC(f)
    models <- list(two = fit_two, one = fit_one, bol = fit_bol)
    aics   <- purrr::map_dbl(models, get_AIC)
    aics[!is.finite(aics)] <- Inf
    
    if (all(is.infinite(aics))) {
      message("⚠ ", yvar, ": no model converged.")
      return(base)
    }
    
    best_name <- names(which.min(aics))
    fit       <- models[[best_name]]
    
    grid <- tibble::tibble(Temperature = seq(T_min, T_max, length.out = 300))
    eta_pred  <- predict_safe(fit, grid)
    grid$Pred <- invlogit(eta_pred)
    grid <- grid %>% dplyr::filter(is.finite(Pred))
    
    if (nrow(grid) > 1) base <- base + geom_line(data = grid, aes(Temperature, Pred), linewidth = 1)
    
    message("✅ ", yvar, ": using ", best_name, " model (logit-bounded).")
    return(base)
  }
  
  df <- df %>%
    dplyr::filter(is.finite(.data[[yvar]]), .data[[yvar]] > 0)
  
  base <- ggplot(df, aes(Temperature, .data[[yvar]])) +
    geom_point(size = 2.5) +
    labs(
      x = if (show_x) expression(Temperature~"(°C)") else NULL,
      y = ylab
    ) +
    isme_theme()
  
  if (nrow(df) < 4) return(base)
  
  T_min     <- min(df$Temperature)
  T_max     <- max(df$Temperature)
  T_opt_obs <- df$Temperature[which.max(df[[yvar]])]
  
  start_two <- list(
    lnR0 = log(max(df[[yvar]])),
    E  = 0.6,
    El = 0.4,  Tl = T_min + 2,
    Eh = 1.5,  Th = T_opt_obs + 3
  )
  fit_two <- try(
    nlsLM(
      as.formula(
        paste0("log(", yvar,
               ") ~ ln_SS_two(Temperature, lnR0, E, El, Tl, Eh, Th)")
      ),
      data   = df,
      start  = start_two,
      lower  = c(lnR0=-Inf, E=0.1, El=0.1, Tl= 0,  Eh=0.1, Th= 0),
      upper  = c(lnR0= Inf, E=2.5, El=2.5, Tl=60, Eh=5.0, Th=60),
      control = nls.lm.control(maxiter = 600)
    ),
    silent = TRUE
  )
  
  start_one <- list(
    lnR0 = log(max(df[[yvar]])),
    E  = 0.6,
    Eh = if (yvar == "growth_C_fg_per_hr") 3.0 else 1.5,
    Th = if (yvar == "growth_C_fg_per_hr") T_opt_obs + 2 else T_opt_obs + 5
  )
  fit_one <- try(
    nlsLM(
      as.formula(
        paste0("log(", yvar,
               ") ~ ln_SS_one(Temperature, lnR0, E, Eh, Th)")
      ),
      data   = df,
      start  = start_one,
      lower  = c(lnR0=-Inf, E=0.1, Eh=0.5, Th=0),
      upper  = c(lnR0= Inf, E=2.5, Eh=6.0, Th=60),
      control = nls.lm.control(maxiter = 500)
    ),
    silent = TRUE
  )
  
  fit_bol <- try(
    nlsLM(
      as.formula(
        paste0("log(", yvar,
               ") ~ ln_boltz(Temperature, lnR0, E)")
      ),
      data  = df,
      start = list(lnR0 = log(max(df[[yvar]])), E = 0.6),
      lower = c(lnR0=-Inf, E=0.1),
      upper = c(lnR0= Inf,  E=2.5),
      control = nls.lm.control(maxiter = 400)
    ),
    silent = TRUE
  )
  
  get_AIC <- function(f) if (inherits(f, "try-error")) Inf else AIC(f)
  
  models <- list(two = fit_two, one = fit_one, bol = fit_bol)
  aics   <- purrr::map_dbl(models, get_AIC)
  aics[!is.finite(aics)] <- Inf
  
  if (all(is.infinite(aics))) {
    message("⚠ ", yvar, ": no model converged.")
    return(base)
  }
  
  best_name <- names(which.min(aics))
  fit       <- models[[best_name]]
  
  grid <- tibble::tibble(Temperature = seq(T_min, T_max, length.out = 300))
  lnpred    <- predict_safe(fit, grid)
  grid$Pred <- safe_exp(lnpred)
  grid <- grid %>% dplyr::filter(is.finite(Pred))
  
  if (nrow(grid) > 1) base <- base + geom_line(data = grid, aes(Temperature, Pred), linewidth = 1)
  
  message("✅ ", yvar, ": using ", best_name, " model (plot only).")
  base
}

fit_T_params_only <- function(df, yvar) {
  
  stats <- tibble::tibble(
    response   = yvar,
    best_model = NA_character_,
    Topt_C     = NA_real_,
    E_eV       = NA_real_,
    E_lwr_eV   = NA_real_,
    E_upr_eV   = NA_real_,
    Eh_eV      = NA_real_,
    Eh_lwr_eV  = NA_real_,
    Eh_upr_eV  = NA_real_
  )
  
  if (yvar == "carbon_use_efficiency") {
    
    df <- df %>%
      dplyr::filter(is.finite(.data[[yvar]])) %>%
      dplyr::mutate(
        cue_adj   = clamp01(.data[[yvar]]),
        cue_logit = logit(cue_adj)
      )
    
    if (nrow(df) < 4) {
      message("⚠ ", yvar, ": fewer than 4 points, skipping TPC fit.")
      return(stats)
    }
    
    T_min     <- min(df$Temperature)
    T_max     <- max(df$Temperature)
    T_opt_obs <- df$Temperature[which.max(df$cue_adj)]
    
    # ------------------ CUE bounds: FORCE crash at Tmax ------------------
    Th_low  <- T_max - 3
    Th_high <- T_max + 0.2
    Eh_low  <- 6
    Eh_high <- 300
    # --------------------------------------------------------------------
    
    start_two <- list(
      eta0 = logit(max(df$cue_adj)),
      E  = 0.6,
      El = 0.4,  Tl = T_min + 2,
      Eh = 10,   Th = T_max - 1
    )
    fit_two <- try(
      nlsLM(
        cue_logit ~ eta_SS_two(Temperature, eta0, E, El, Tl, Eh, Th),
        data   = df,
        start  = start_two,
        lower  = c(eta0=-Inf, E=0.1, El=0.1, Tl=0,  Eh=Eh_low,  Th=Th_low),
        upper  = c(eta0= Inf, E=2.5, El=2.5, Tl=60, Eh=Eh_high, Th=Th_high),
        control = nls.lm.control(maxiter = 1200)
      ),
      silent = TRUE
    )
    
    start_one <- list(
      eta0 = logit(max(df$cue_adj)),
      E  = 0.6,
      Eh = 10,
      Th = T_max - 1
    )
    fit_one <- try(
      nlsLM(
        cue_logit ~ eta_SS_one(Temperature, eta0, E, Eh, Th),
        data   = df,
        start  = start_one,
        lower  = c(eta0=-Inf, E=0.1, Eh=Eh_low,  Th=Th_low),
        upper  = c(eta0= Inf, E=2.5, Eh=Eh_high, Th=Th_high),
        control = nls.lm.control(maxiter = 1200)
      ),
      silent = TRUE
    )
    
    fit_bol <- try(
      nlsLM(
        cue_logit ~ eta_boltz(Temperature, eta0, E),
        data  = df,
        start = list(eta0 = logit(max(df$cue_adj)), E = 0.6),
        lower = c(eta0=-Inf, E=0.1),
        upper = c(eta0= Inf, E=2.5),
        control = nls.lm.control(maxiter = 800)
      ),
      silent = TRUE
    )
    
    get_AIC <- function(f) if (inherits(f, "try-error")) Inf else AIC(f)
    models  <- list(two = fit_two, one = fit_one, bol = fit_bol)
    aics    <- purrr::map_dbl(models, get_AIC)
    aics[!is.finite(aics)] <- Inf
    
    if (all(is.infinite(aics))) {
      message("⚠ ", yvar, ": no TPC model converged.")
      return(stats)
    }
    
    best_name        <- names(which.min(aics))
    fit              <- models[[best_name]]
    stats$best_model <- best_name
    
    grid <- tibble::tibble(Temperature = seq(T_min, T_max, length.out = 300))
    eta_pred  <- predict_safe(fit, grid)
    grid$Pred <- invlogit(eta_pred)
    grid <- grid %>% dplyr::filter(is.finite(Pred))
    
    if (nrow(grid) > 1) stats$Topt_C <- grid$Temperature[which.max(grid$Pred)]
    
    pars <- coef(summary(fit))
    
    if ("E" %in% rownames(pars)) {
      E_est <- pars["E", "Estimate"]
      E_se  <- pars["E", "Std. Error"]
      stats$E_eV     <- E_est
      stats$E_lwr_eV <- E_est - 1.96 * E_se
      stats$E_upr_eV <- E_est + 1.96 * E_se
    }
    if ("Eh" %in% rownames(pars)) {
      Eh_est <- pars["Eh", "Estimate"]
      Eh_se  <- pars["Eh", "Std. Error"]
      stats$Eh_eV     <- Eh_est
      stats$Eh_lwr_eV <- Eh_est - 1.96 * Eh_se
      stats$Eh_upr_eV <- Eh_est + 1.96 * Eh_se
    }
    
    message("✅ ", yvar, ": best model = ", best_name,
            ", Topt ≈ ", round(stats$Topt_C, 1), " °C (logit-bounded)")
    return(stats)
  }
  
  df <- df %>% dplyr::filter(is.finite(.data[[yvar]]), .data[[yvar]] > 0)
  
  if (nrow(df) < 4) {
    message("⚠ ", yvar, ": fewer than 4 points, skipping TPC fit.")
    return(stats)
  }
  
  T_min     <- min(df$Temperature)
  T_max     <- max(df$Temperature)
  T_opt_obs <- df$Temperature[which.max(df[[yvar]])]
  
  start_two <- list(
    lnR0 = log(max(df[[yvar]])),
    E  = 0.6,
    El = 0.4,  Tl = T_min + 2,
    Eh = 1.5,  Th = T_opt_obs + 3
  )
  fit_two <- try(
    nlsLM(
      as.formula(
        paste0("log(", yvar,
               ") ~ ln_SS_two(Temperature, lnR0, E, El, Tl, Eh, Th)")
      ),
      data   = df,
      start  = start_two,
      lower  = c(lnR0 = -Inf, E = 0.1, El = 0.1, Tl = 0,  Eh = 0.1, Th = 0),
      upper  = c(lnR0 =  Inf, E = 2.5, El = 2.5, Tl = 60, Eh = 5.0, Th = 60),
      control = nls.lm.control(maxiter = 600)
    ),
    silent = TRUE
  )
  
  start_one <- list(
    lnR0 = log(max(df[[yvar]])),
    E  = 0.6,
    Eh = if (yvar == "growth_r_per_hr") 3.0 else 1.5,
    Th = if (yvar == "growth_r_per_hr") T_opt_obs + 2 else T_opt_obs + 5
  )
  fit_one <- try(
    nlsLM(
      as.formula(
        paste0("log(", yvar,
               ") ~ ln_SS_one(Temperature, lnR0, E, Eh, Th)")
      ),
      data   = df,
      start  = start_one,
      lower  = c(lnR0 = -Inf, E = 0.1, Eh = 0.5, Th = 0),
      upper  = c(lnR0 =  Inf, E = 2.5, Eh = 6.0, Th = 60),
      control = nls.lm.control(maxiter = 500)
    ),
    silent = TRUE
  )
  
  fit_bol <- try(
    nlsLM(
      as.formula(
        paste0("log(", yvar,
               ") ~ ln_boltz(Temperature, lnR0, E)")
      ),
      data  = df,
      start = list(lnR0 = log(max(df[[yvar]])), E = 0.6),
      lower = c(lnR0 = -Inf, E = 0.1),
      upper = c(lnR0 =  Inf, E = 2.5),
      control = nls.lm.control(maxiter = 400)
    ),
    silent = TRUE
  )
  
  get_AIC <- function(f) if (inherits(f, "try-error")) Inf else AIC(f)
  models  <- list(two = fit_two, one = fit_one, bol = fit_bol)
  aics    <- purrr::map_dbl(models, get_AIC)
  aics[!is.finite(aics)] <- Inf
  
  if (all(is.infinite(aics))) {
    message("⚠ ", yvar, ": no TPC model converged.")
    return(stats)
  }
  
  best_name        <- names(which.min(aics))
  fit              <- models[[best_name]]
  stats$best_model <- best_name
  
  grid <- tibble::tibble(Temperature = seq(T_min, T_max, length.out = 300))
  lnpred    <- predict_safe(fit, grid)
  grid$Pred <- safe_exp(lnpred)
  grid      <- grid %>% dplyr::filter(is.finite(Pred))
  
  if (nrow(grid) > 1) stats$Topt_C <- grid$Temperature[which.max(grid$Pred)]
  
  pars <- coef(summary(fit))
  
  if ("E" %in% rownames(pars)) {
    E_est <- pars["E", "Estimate"]
    E_se  <- pars["E", "Std. Error"]
    stats$E_eV     <- E_est
    stats$E_lwr_eV <- E_est - 1.96 * E_se
    stats$E_upr_eV <- E_est + 1.96 * E_se
  }
  
  if ("Eh" %in% rownames(pars)) {
    Eh_est <- pars["Eh", "Estimate"]
    Eh_se  <- pars["Eh", "Std. Error"]
    stats$Eh_eV     <- Eh_est
    stats$Eh_lwr_eV <- Eh_est - 1.96 * Eh_se
    stats$Eh_upr_eV <- Eh_est + 1.96 * Eh_se
  }
  
  message("✅ ", yvar, ": best model = ", best_name,
          ", Topt ≈ ", round(stats$Topt_C, 1), " °C")
  
  stats
}

results_filtered <- results %>%
  dplyr::filter(fit_ok,
                growth_C_fg_per_hr > 0,
                resp_C_fg_per_hr   > 0)

print(paste("Number of rows in results_filtered:", nrow(results_filtered)))

p_growth <- fit_T_plot(results_filtered, "growth_C_fg_per_hr",
                       expression(Growth~(fg~C~h^{-1})), TRUE)
p_resp   <- fit_T_plot(results_filtered, "resp_C_fg_per_hr",
                       expression(Respiration~(fg~C~h^{-1})), TRUE)
p_cue    <- fit_T_plot(results_filtered, "carbon_use_efficiency",
                       "Carbon Use Efficiency", TRUE)

tpc_params <- dplyr::bind_rows(
  fit_T_params_only(
    results_filtered %>%
      dplyr::filter(r_per_hour > 0) %>%
      dplyr::mutate(growth_r_per_hr = r_per_hour),
    "growth_r_per_hr"
  ) %>% dplyr::mutate(response = "growth_r_per_hr"),
  
  fit_T_params_only(results_filtered, "resp_C_fg_per_hr"),
  fit_T_params_only(results_filtered, "carbon_use_efficiency")
)

readr::write_csv(
  tpc_params,
  file.path(tables_dir, "SharpeSchoolfield_Temperature_Params_NEWformula.csv")
)

combo_plot <- patchwork::wrap_plots(p_growth, p_resp, p_cue, ncol = 3) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    tag_levels = "A",
    tag_prefix = "(", tag_suffix = ")"
  ) &
  theme(plot.margin = margin(5, 10, 5, 10))

ggsave(
  file.path(plots_dir, "Fig_7_SharpeSchoolfield_Temperature_Fits_NEWformula.pdf"),
  combo_plot, width = 12, height = 4, dpi = 600
)

##################################################################################
#                                   End of script                                #
##################################################################################
