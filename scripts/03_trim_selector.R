# =============================================================================
# 03_trim_selector.R
# Click on each O2 curve to set its fit-window START and END (manual trimming).
# -----------------------------------------------------------------------------
# - Shows one curve at a time (full O2 vs Time trace).
# - Blue curve = the respiration model fitted to the current window (solid inside
#   the window, dashed where it is projected beyond it) - this is the model 03
#   will fit. Orange dashed line = suggested end (a guide, computed live from your
#   start). Tick "Show model guide curve" to toggle the blue line.
# - Choose "Click sets: Start" or "End", then click on the plot at the time you
#   want. Green line = your manual start, red line = your manual end.
# - Navigate curves with the dropdown or Prev/Next.
# - The set of manual windows is written to tables/manual_fit_windows.csv and
#   printed as a ready-to-paste MANUAL_FIT_WINDOWS block.
# - SELECT / DESELECT SAMPLES: tick "Don't include this sample" to discard a
#   curve (or untick to keep it). The discard list is written to
#   tables/plot_exclude_points.csv and printed as a ready-to-paste PLOT_EXCLUDE_POINTS
#   block for config.R. The currently-discarded samples (in that CSV) are
#   RESTORED on startup. config.R loads it when USE_APP_TRIM_FILES <- TRUE.
# - Curves are keyed by Taxon + Replicate (OTU = Candida isolate code 1..18).
# - Persists: reloads tables/manual_fit_windows.csv + plot_exclude_points.csv.
#
# RUN:
#   RStudio: open this file -> "Run App".
#   Terminal: Rscript scripts/03_trim_selector.R
#
# Needs 02_longdata.R + 03_trimming.R to have run (for the data + metadata).
#
# PATHS: this app anchors strictly to THIS project. base_dir is the PARENT of
# the scripts/ folder (identical rule to config.R), so every file it reads
# or writes stays inside this project's tables/ folder. There are no hard-coded
# external paths, so nothing is ever generated outside this project.
# =============================================================================

# ---- Locate data ------------------------------------------------------------
# Robustly find the folder THIS script lives in — whether launched via RStudio
# "Run App", source(), or Rscript — then derive the project root from it.
.script_dir <- local({
  d <- tryCatch({
    if (requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable() &&
        nzchar(rstudioapi::getActiveDocumentContext()$path)) {
      dirname(rstudioapi::getActiveDocumentContext()$path)
    } else NA_character_
  }, error = function(e) NA_character_)
  if (length(d) == 0 || is.na(d) || !nzchar(d)) {
    d <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  }
  if (length(d) == 0 || is.na(d) || !nzchar(d)) d <- getwd()
  normalizePath(d, mustWork = FALSE)
})

# base_dir = parent of scripts/ (same rule as config.R); tables/ lives there.
base_dir   <- dirname(.script_dir)
tables_dir <- file.path(base_dir, "results", "tables")

# Fallback ONLY to tables/ dirs relative to the working directory — never to any
# hard-coded external location — so files can only ever land inside this project.
if (!dir.exists(tables_dir)) {
  cand <- c(file.path(base_dir, "results", "tables"), file.path(getwd(), "results", "tables"), file.path(getwd(), "..", "results", "tables"))
  hit  <- cand[dir.exists(cand)]
  if (length(hit)) tables_dir <- normalizePath(hit[1], mustWork = FALSE)
}

find_file <- function(fname) {
  p <- file.path(tables_dir, fname)
  if (file.exists(p)) normalizePath(p, mustWork = FALSE) else NA_character_
}

long_path <- find_file("Oxygen_All_Long.csv")
meta_path <- find_file("Oxygen_Trimmed_Series_Metadata.csv")
if (is.na(long_path)) stop("Could not find Oxygen_All_Long.csv in ", tables_dir,
                           ". Run 02_longdata.R first.")
if (is.na(meta_path)) stop("Could not find Oxygen_Trimmed_Series_Metadata.csv in ",
                           tables_dir, ". Run 03_trimming.R first.")
message("Trim selector reading/writing tables in: ", tables_dir)
out_csv <- file.path(tables_dir, "manual_fit_windows.csv")

# ---- Packages ---------------------------------------------------------------
need <- c("shiny", "ggplot2", "readr", "dplyr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0) {
  stop("Install missing package(s) first:\n  install.packages(c(",
       paste(sprintf('"%s"', miss), collapse = ", "), "))")
}
library(shiny); library(ggplot2); suppressPackageStartupMessages(library(dplyr))

# Model used for the guide curve (identical to the pipeline's resp_model).
# IMPORTANT: the guide fit here (start values + bounds in fit_resp / winfit below)
# is kept IDENTICAL to 05_oxygen_fits.R's fit_one, so the blue guide line is the
# exact model 03 will fit. If you change bounds here, change them in 03 too.
resp_model  <- function(r, K, t, O2_0) O2_0 + (K / r) * (1 - exp(r * t))
has_minpack <- requireNamespace("minpack.lm", quietly = TRUE)

# Suggested exponential-phase bounds from the RAW trace (reference markers only —
# they do NOT constrain your selection). start = time of the O2 peak (where the
# draw-down begins); end = where O2 has fallen to within 5% of its post-peak
# minimum (draw-down essentially complete). No smoothing is used.
exp_phase_bounds <- function(tt, yy) {
  ok <- is.finite(tt) & is.finite(yy)
  tt <- tt[ok]; yy <- yy[ok]
  if (length(tt) < 5) return(c(NA_real_, NA_real_))
  ord <- order(tt); tt <- tt[ord]; yy <- yy[ord]
  win <- which(tt <= min(tt) + 400)          # peak is always early; avoid late noise
  if (!length(win)) win <- seq_along(tt)
  pk      <- win[which.max(yy[win])]
  peak_t  <- tt[pk]; o2pk <- yy[pk]
  seg     <- which(tt >= peak_t)
  o2min   <- min(yy[seg])
  thr     <- o2min + 0.05 * (o2pk - o2min)
  hit     <- seg[yy[seg] <= thr]
  end_t   <- if (length(hit)) tt[hit[1]] else tt[length(tt)]
  c(peak_t, end_t)
}

# 95%-down endpoint measured from a GIVEN start time: take O2 at that start,
# find the floor after it, and return the first time O2 falls to within 5% of
# that floor. Lets the end suggestion follow wherever you click the start.
suggest_end_95 <- function(tt, yy, start_t, frac = 0.95) {
  ok <- is.finite(tt) & is.finite(yy)
  tt <- tt[ok]; yy <- yy[ok]
  ord <- order(tt); tt <- tt[ord]; yy <- yy[ord]
  seg <- which(tt >= start_t)
  if (length(seg) < 3) return(NA_real_)
  o2_start <- yy[seg[1]]
  o2min    <- min(yy[seg])
  thr      <- o2min + (1 - frac) * (o2_start - o2min)   # frac of the way down
  hit      <- seg[yy[seg] <= thr]
  if (length(hit)) tt[hit[1]] else tt[length(tt)]
}

# Fit resp_model to a (time, oxygen) vector; return coefs + R2 + RMSE.
fit_resp <- function(tt, yy) {
  out <- list(ok = FALSE, co = NULL, r2 = NA_real_, rmse = NA_real_, n = length(tt))
  if (!has_minpack || length(tt) < 6) return(out)
  t0     <- tt - min(tt)
  slope0 <- suppressWarnings(stats::median(diff(yy) / diff(t0), na.rm = TRUE))
  K0     <- if (is.finite(slope0)) max(abs(slope0), 1e-6) else 1e-3
  ft <- try(minpack.lm::nlsLM(
    Oxygen ~ resp_model(r, K, Time0, O2_0),
    data  = data.frame(Time0 = t0, Oxygen = yy),
    start = list(r = 1e-3, K = K0, O2_0 = yy[1]),
    lower = c(r = 1e-6, K = 1e-10, O2_0 = min(yy) - 1),
    upper = c(r = 0.15, K = 1,     O2_0 = max(yy) + 1),
    control = minpack.lm::nls.lm.control(maxiter = 200)
  ), silent = TRUE)
  if (inherits(ft, "try-error")) return(out)
  co    <- coef(ft)
  preds <- resp_model(co[["r"]], co[["K"]], t0, co[["O2_0"]])
  resid <- yy - preds
  sstot <- sum((yy - mean(yy))^2)
  out$ok <- TRUE; out$co <- co
  out$rmse <- sqrt(mean(resid^2))
  out$r2   <- if (sstot > 1e-12) 1 - sum(resid^2) / sstot else NA_real_
  out
}

# =============================================================================
# STABLE-r WINDOW OPTIMISER  (curve-intrinsic)
# =============================================================================
# The plain "Auto" detector below is a GEOMETRY rule: it asks "where does the
# curve bend down and level off?". It never asks whether the r it produces is
# trustworthy - and a clean fit to the WRONG window gives a precise but wrong r.
#
# This optimiser instead asks: "which window gives an r I can BELIEVE?"
#   - it tries many candidate windows (several starts x several end fractions),
#   - fits the exponential model in each and keeps the ones that fit well,
#   - and picks the window whose r sits at the CENTRE of the good windows.
# If r barely moves as the window is nudged, that r is real (stability_cv small).
# If r swings wildly, this curve's r is not a measurement - and it tells you so.
#
# It judges each curve on its OWN trace only. It never looks at the temperature
# curve, so it cannot manufacture a nice Sharpe-Schoolfield shape.
# Start candidates: how far AFTER the O2 peak the fit window may begin.
#
# KEEP THIS NARROW. The exponential ACCELERATION - the entire signal for r -
# sits in the shallow stretch immediately after the O2 peak. Allow the window to
# start much later and it can skip that region, leaving a short, nearly straight
# window where r cannot be measured; the curve is then wrongly flagged ZERO
# GROWTH (this happened to Clade1_2069 T=32 R3 when the range was widened).
#
# I briefly widened this to ~120 points to chase glab's small K. Widening the
# START RANGE was the wrong fix and was reverted. But the underlying suspicion
# was right after all: glab's K WAS a windowing problem - just at the END, not
# the start. See STAB_RT_TARGETS below.
STAB_START_OFFSETS <- c(0, 3, 6, 9, 12)
STAB_END_FRACS     <- c(0.35, 0.45, 0.55, 0.65, 0.75)  # of the total drawdown
STAB_R2_MIN        <- 0.90
STAB_MIN_PTS       <- 8

# Window ends offered by TARGET rt, as well as by fraction-of-drawdown above.
# Without these, a slow-respiring / fast-growing organism can only ever be given
# long windows (every drawdown target lies far out in its trace), and K ends up
# extrapolated back through e^4 or more.
#
# ONLY 3.0. Measured on glab (n = 175 curves), sd(log K) by target:
#     none (rt 3.95) 0.582 | 3.0 (rt 2.16) 0.365 | 2.5 (rt 1.89) 0.447
#                          | 2.0 (rt 1.68) 0.598 | 1.5 (rt 1.50) 0.723
# A target BELOW 3 cuts the window before enough curvature has developed, r goes
# soft, and the fit is WORSE than doing nothing at all. The two-pass construction
# overshoots (the provisional r is biased high), so a target of 3.0 lands at an
# achieved rt of ~2.2 - which is the optimum. Do not "improve" this by lowering it.
# Set to numeric(0) to restore the old behaviour exactly.
STAB_RT_TARGETS <- c(3.0)

# ---- Cap on how far the exponential may extrapolate --------------------------
# rt = r * window = the number of E-FOLDINGS of growth the fit spans.
#   rt ~ 2  -> biomass grows ~7-fold across the window (healthy exponential phase)
#   rt ~ 4.5-> biomass grows ~90-fold  (the model is extrapolating hard)
# The respiration model assumes UNLIMITED exponential growth. Past ~3 e-foldings
# that is untenable - the cells are running out of oxygen and substrate - and the
# fit inflates the biomass integral, which deflates per-cell respiration
# (respiration = K / N0 = total O2 / biomass integral).
#
# In this dataset the clades and para sit at rt = 1.8-2.7 with essentially NO
# windows above 4, so this cap does nothing to them. glab was running at rt = 4.5
# (108 of 179 windows above 4), which is exactly the over-extrapolation you can
# see as a "huge range" on its curves. A uniform rule that only bites where the
# assumption is actually being stretched.
#
# NOTE: this filters the window END only. It never moves the START, so it cannot
# skip the curvature that determines r (the mistake the plateau guard made).
# Set to Inf to disable.
# 3.0, not 3.5: measured on this dataset, glab's K noise is minimised at rt ~ 2.2,
# and 3.0 is the cap that lets the rt-targeted candidates (2.0/2.5) actually win.
# Only 7% of non-glab curves are above 3.0, and trimming them changes nothing
# measurable (Clade4: sd(log K) 0.218 -> 0.219, median respiration unchanged).
STAB_MAX_RT <- 6.0   # raised from the fungal 3.0 for fast bacterial growers
# (e.g. Yersinia ~0.046 min^-1) that fit cleanly (R2>=0.999) out to ~75% drawdown
# but at rt 4-5; a cap of 3 clipped them to ~15% of the curve. STAB_END_FRACS
# (0.75) is the real end-limiter; 6 admits fast growers while still guarding
# against runaway exponential extrapolation.

# Plateau guard: require the O2 slope at the window start to be at least this
# fraction of the curve's steepest slope.
#
# >>> DEFAULT 0 (OFF), ON PURPOSE. <<<
# It looks sensible - "don't start the fit on a flat shoulder" - but it BACKFIRES.
# The shallow stretch just after the O2 peak is precisely where the exponential
# ACCELERATION lives, and that acceleration is the entire signal for r. Skipping
# it leaves a short, nearly straight window where r cannot be measured, and the
# curve then gets wrongly flagged ZERO GROWTH even though it grew perfectly well
# (seen on Clade1_2069 T=32 R3: window pushed to 157-221 min, missing the curved
# 85-160 min region entirely).
#
# The two goals genuinely conflict:
#   - r  needs the curved region right after the peak (shallow slope).
#   - K  is the consumption rate AT the window start, so it is ~0 if you start there.
# You cannot fix both by MOVING the window start.
#
# CORRECTION (this comment used to blame the N0 back-projection and recommend
# N0_BACKPROJECT <- FALSE; that was wrong). glab's bad K is a windowing problem
# after all - but it is the window END, not the start. Its windows spanned rt ~ 4,
# so K was extrapolated backwards through e^4 (~55x). Capping rt fixes it:
# sd(log K) 0.58 -> 0.37, and glab's respiration moves by 2.5x. See
# STAB_RT_TARGETS / STAB_MAX_RT above. Do NOT set N0_BACKPROJECT <- FALSE.
#
# Set to ~0.05 only if you have a specific curve with a genuinely dead flat start.
STAB_MIN_START_SLOPE_FRAC <- 0

stable_window <- function(tt, yy) {
  ok <- is.finite(tt) & is.finite(yy)
  tt <- tt[ok]; yy <- yy[ok]
  ord <- order(tt); tt <- tt[ord]; yy <- yy[ord]
  n <- length(tt)
  if (n < STAB_MIN_PTS) return(NULL)

  win <- which(tt <= min(tt) + 400)
  if (!length(win)) win <- seq_along(tt)
  pk <- win[which.max(yy[win])]

  # Steepest decline anywhere after the peak - used to judge whether a candidate
  # start is on the flat plateau or in the real consumption phase.
  slope_of <- function(i) {
    j <- min(i + 5L, n)
    if (j <= i) return(0)
    (yy[j] - yy[i]) / (tt[j] - tt[i])          # negative while O2 falls
  }
  max_slope <- min(vapply(seq(pk, max(pk, n - 6L)), slope_of, numeric(1)))  # most negative

  cand <- list()

  # Fit one window [s, e] and push it onto the candidate list if it is usable.
  try_window <- function(s, e) {
    if (e - s + 1L < STAB_MIN_PTS) return(NULL)
    ft <- fit_resp(tt[s:e], yy[s:e])
    if (!isTRUE(ft$ok)) return(NULL)
    rr <- tryCatch(unname(ft$co[["r"]]), error = function(e) NA_real_)
    if (!is.finite(ft$r2) || ft$r2 < STAB_R2_MIN) return(NULL)
    if (!is.finite(rr) || rr <= 0) return(NULL)
    data.frame(start = tt[s], end = tt[e], r = rr, r2 = ft$r2,
               rt = rr * (tt[e] - tt[s]))   # e-foldings this window spans
  }

  for (so in STAB_START_OFFSETS) {
    s <- pk + so
    if (s > n - STAB_MIN_PTS) next
    o0 <- yy[s]; omin <- min(yy[s:n]); dd <- o0 - omin
    if (!is.finite(dd) || dd <= 0) next
    # Reject starts sitting on the PLATEAU: at the start of the window O2 must
    # already be falling at a decent fraction of its steepest rate. Opening the
    # window on the flat shoulder drives K -> 0 and destroys respiration.
    if (is.finite(max_slope) && max_slope < 0) {
      if (slope_of(s) > STAB_MIN_START_SLOPE_FRAC * max_slope) next
    }

    # ---- (i) ends chosen by FRACTION OF DRAWDOWN (the original rule) --------
    r_prov <- NA_real_
    for (f in STAB_END_FRACS) {
      thr <- o0 - f * dd
      hit <- which(yy[s:n] <= thr)
      e <- if (length(hit)) (s + hit[1] - 1L) else n
      cc <- try_window(s, e)
      if (is.null(cc)) next
      cand[[length(cand) + 1]] <- cc
      r_prov <- cc$r                       # keep the last (widest) r as a seed
    }

    # ---- (ii) ends chosen by TARGET rt  --- THE FIX -------------------------
    # The drawdown-fraction rule above was the ONLY way a window end could ever
    # be chosen. That is fine for a brisk respirer, but for an organism that
    # grows fast and consumes O2 slowly (C. glabrata here) every drawdown target
    # lands far out in the trace: 67% of its windows ran past rt = 3, median
    # rt = 3.6, versus 1.8-2.7 for every other taxon.
    #
    # Why that matters: K is the O2 consumption rate AT THE WINDOW START, and the
    # model infers it by extrapolating BACKWARDS across the window. rt is exactly
    # how far back that is, in e-foldings. At rt = 4 the extrapolation spans
    # e^4 ~ 55x, and K comes out both noisy AND BIASED - on this dataset,
    # capping glab at rt ~ 2.2 cut sd(log K) by 37% and moved its respiration by
    # 2.5x. We would have published that bias.
    #
    # So: also offer ends defined directly by rt. Two-pass, because rt depends on
    # r which depends on the window: take the provisional r from the widest
    # drawdown window, cut to the length that hits the target, refit.
    #
    # This is a UNIFORM rule, not a glab special case. It touches 7% of non-glab
    # curves and changes them by nothing measurable (Clade4, the most affected at
    # 16%: sd(log K) 0.218 -> 0.219; median respiration unchanged to 4 s.f.).
    if (length(STAB_RT_TARGETS) && is.finite(r_prov) && r_prov > 0) {
      for (tgt in STAB_RT_TARGETS) {
        wlen <- tgt / r_prov
        ok <- which(tt[s:n] <= tt[s] + wlen)
        if (!length(ok)) next
        e <- s + max(ok) - 1L
        cc <- try_window(s, e)
        if (!is.null(cc)) cand[[length(cand) + 1]] <- cc
      }
    }
  }
  if (!length(cand)) return(NULL)
  cd <- do.call(rbind, cand)
  cd <- cd[!duplicated(cd[, c("start", "end")]), , drop = FALSE]

  # Drop windows that extrapolate the exponential too far (STAB_MAX_RT).
  #
  # The old code had a fallback here: "if EVERY window is over the cap, keep the
  # 3 shortest anyway." That fallback is why glab sailed straight past the cap -
  # 67% of its windows were above rt = 3 and it simply kept them. With the
  # rt-targeted ends added above, a short window now EXISTS for these curves, so
  # the cap can finally be enforced. The fallback is kept only as a last resort
  # (a curve so short that even rt = 2 cannot be reached), and it now takes the
  # single shortest window rather than three.
  if (is.finite(STAB_MAX_RT)) {
    keep <- is.finite(cd$rt) & cd$rt <= STAB_MAX_RT
    if (any(keep)) {
      cd <- cd[keep, , drop = FALSE]
    } else {
      cd <- cd[order(cd$rt), , drop = FALSE]
      cd <- cd[1L, , drop = FALSE]          # the single shortest extrapolation
    }
  }

  medr <- stats::median(cd$r)
  madr <- 1.4826 * stats::median(abs(cd$r - medr))
  cv   <- if (is.finite(medr) && medr > 0) madr / medr else NA_real_
  best <- cd[which.min(abs(cd$r - medr)), , drop = FALSE]

  # ---- ZERO GROWTH: two independent ways a curve can have no usable r --------
  # (a) NO CURVATURE. rt = r * window is the dimensionless curvature. If it is
  #     small the trace is effectively a STRAIGHT decreasing line: constant
  #     respiration, no growth. (An R2/RMSE check cannot see this - a straight
  #     line fits the model beautifully.)
  # (b) NO SIGNAL. The culture barely consumed any oxygen, so there is nothing to
  #     fit and the "curvature" is just noise. This is the dead-culture case seen
  #     at 42-44 degC: only ~0.3-1.7 mg/L of O2 drawn down (the typical curve uses
  #     ~4.7), yet the model reports near-maximal growth. Impossible - you cannot
  #     build biomass without paying for it in respiration. Test (a) MISSES these
  #     (their rt is 0.7-2.0), which is why we need this second test.
  s_i <- which.min(abs(tt - best$start[1]))
  e_i <- which.min(abs(tt - best$end[1]))

  # Oxygen used INSIDE the chosen window - reported for information only.
  drawdown_win <- yy[s_i] - min(yy[s_i:e_i], na.rm = TRUE)

  # Oxygen used over the WHOLE trace from the peak onwards. THIS is the one the
  # dead-culture test must use.
  #
  # BUG THIS FIXES: no_signal used to be measured inside the fit window. That was
  # harmless only because the windows used to be long. Once STAB_RT_TARGETS made
  # windows shorter, the slowest respirer in the set (glab) stopped consuming
  # 2 mg/L *within its window* and was flagged NO GROWTH - a healthy culture
  # condemned because we shortened its window. "Dead culture" is a property of the
  # TRACE, not of whichever window we happened to pick, so it is measured on the
  # trace.
  drawdown <- yy[s_i] - min(yy[s_i:length(yy)], na.rm = TRUE)

  rt <- best$r[1] * (best$end[1] - best$start[1])
  no_curv   <- is.finite(rt) && rt <= CURV_MIN_RT
  no_signal <- is.finite(drawdown) && drawdown < MIN_O2_DRAWDOWN

  list(start = best$start[1], end = best$end[1], r = best$r[1], r2 = best$r2[1],
       n_good = nrow(cd), stability_cv = cv, rt = rt,
       drawdown = drawdown, drawdown_win = drawdown_win,
       no_curv = no_curv, no_signal = no_signal,
       no_growth = no_curv || no_signal)
}

# ZERO-GROWTH thresholds.
#   CURV_MIN_RT     : below this, r * window cannot be told from a straight line.
#   MIN_O2_DRAWDOWN : below this (mg/L of O2 consumed over the window) there is
#                     no metabolic signal at all - the culture is dead and any r
#                     is fitted noise. Median curve in this dataset uses ~4.7.
CURV_MIN_RT     <- 0.5
MIN_O2_DRAWDOWN <- 2.0

# Auto-detect a window. Search START times from the O2 peak forward (up to
# peak + start_search minutes, never before the peak), and for each start take
# the LONGEST window that keeps R2 >= r2_target and RMSE <= rmse_max (>= min_pts
# points). Pick the overall longest such window (tie-break: higher R2). Falls
# back to a peak-anchored best-R2 window if none qualify. Returns c(start, end).
auto_detect_window <- function(tt, yy, r2_target, rmse_max, min_pts,
                               start_search = 60, start_mode = "peak",
                               start_drawdown_frac = 0.05,
                               t_peak_min = 10, t_peak_max = 200,
                               start_step = 15, end_step = 6) {
  ok <- is.finite(tt) & is.finite(yy)
  tt <- tt[ok]; yy <- yy[ok]
  ord <- order(tt); tt <- tt[ord]; yy <- yy[ord]
  if (length(tt) < min_pts + 1) return(c(NA_real_, NA_real_))
  win <- which(tt >= t_peak_min & tt <= t_peak_max)
  if (!length(win)) win <- seq_along(tt)
  pk_idx <- win[which.max(yy[win])]
  peak   <- tt[pk_idx]
  o2_pk  <- yy[pk_idx]

  if (identical(start_mode, "drawdown")) {
    # Fixed rule: start where O2 has fallen start_drawdown_frac of its peak->min
    # drop. Same point for every curve (no per-curve floating).
    seg    <- which(tt >= peak)
    o2_min <- min(yy[seg])
    target <- o2_pk - start_drawdown_frac * (o2_pk - o2_min)
    hit    <- seg[yy[seg] <= target]
    starts <- if (length(hit)) tt[hit[1]] else peak
  } else {
    # Peak mode: candidate starts from peak .. peak + start_search.
    s_pool <- tt[tt >= peak & tt <= peak + max(0, start_search)]
    if (!length(s_pool)) s_pool <- peak
    starts <- peak; last <- peak
    for (s in s_pool) if (s - last >= start_step) { starts <- c(starts, s); last <- s }
  }

  best <- list(start = NA_real_, end = NA_real_, len = -Inf, r2 = -Inf)
  fb   <- list(start = NA_real_, end = NA_real_, r2 = -Inf)   # best-R2 fallback
  for (s in starts) {
    after <- which(tt > s)
    if (length(after) < min_pts) next
    cand <- after[seq(min_pts, length(after), by = end_step)]
    for (ei in cand) {
      idx <- which(tt >= s & tt <= tt[ei])
      if (length(idx) < min_pts) next
      fr <- fit_resp(tt[idx], yy[idx])
      if (!isTRUE(fr$ok) || !is.finite(fr$r2)) next
      pass <- fr$r2 >= r2_target &&
        (!is.finite(rmse_max) || (is.finite(fr$rmse) && fr$rmse <= rmse_max))
      len <- tt[ei] - s
      if (pass && (len > best$len || (len == best$len && fr$r2 > best$r2))) {
        best <- list(start = s, end = tt[ei], len = len, r2 = fr$r2)
      }
      if (fr$r2 > fb$r2) { fb$r2 <- fr$r2; fb$start <- s; fb$end <- tt[ei] }
    }
  }
  if (is.finite(best$len) && best$len > 0) return(c(best$start, best$end))
  c(fb$start, fb$end)
}

# ---- Load -------------------------------------------------------------------
# NOTE: in THIS project the second experimental dimension is OTU (the unique
# Candida isolate code, integer 1..18). Curves are keyed by Taxon + Replicate,
# matching Oxygen_All_Long.csv and PLOT_EXCLUDE_POINTS in config.R.
long <- readr::read_csv(long_path, show_col_types = FALSE) %>%
  dplyr::mutate(Taxon = as.character(Taxon),
                Replicate = toupper(as.character(Replicate)),
                key = paste(Taxon, Replicate, sep = "_"))

meta <- readr::read_csv(meta_path, show_col_types = FALSE) %>%
  dplyr::mutate(Taxon = as.character(Taxon),
                Replicate = toupper(as.character(Replicate)),
                key = paste(Taxon, Replicate, sep = "_"),
                main_run_start_time = as.numeric(main_run_start_time),
                steepest_drop_time  = as.numeric(steepest_drop_time)) %>%
  dplyr::select(key, auto_start = main_run_start_time, auto_end = steepest_drop_time)

# Taxon is the display name directly (no OTU-code indirection in this project).
curves <- long %>%
  dplyr::distinct(Taxon, Replicate, key) %>%
  dplyr::arrange(Taxon, Replicate) %>%
  dplyr::mutate(label = sprintf("%s  |  %s", Taxon, Replicate))

# preload existing manual windows
init <- data.frame(key = character(0), fit_start = numeric(0), fit_end = numeric(0),
                   stringsAsFactors = FALSE)
if (file.exists(out_csv)) {
  prev <- tryCatch(readr::read_csv(out_csv, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(prev) && all(c("Taxon", "Replicate") %in% names(prev))) {
    init <- prev %>%
      dplyr::mutate(key = paste(as.character(Taxon), toupper(Replicate), sep = "_"),
                    fit_start = suppressWarnings(as.numeric(fit_start)),
                    fit_end   = suppressWarnings(as.numeric(fit_end))) %>%
      dplyr::select(key, fit_start, fit_end)
  }
}

# preload existing DISCARDED samples so they are restored on startup. Read from
# tables/plot_exclude_points.csv (seeded from PLOT_EXCLUDE_POINTS in 05_oxygen_fits.R).
excl_csv <- file.path(tables_dir, "plot_exclude_points.csv")
excl_init <- character(0)
if (file.exists(excl_csv)) {
  pe <- tryCatch(readr::read_csv(excl_csv, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(pe) && all(c("Taxon", "Replicate") %in% names(pe))) {
    excl_init <- paste(as.character(pe$Taxon), toupper(pe$Replicate), sep = "_")
  }
}

# ---- UI ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Trim selector — click to set each curve's fit window"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      # Filter the curve list down first.
      selectInput("flt_iso", "Taxon:",
                  choices = c("All", sort(unique(curves$Taxon))),
                  selected = "All"),
      actionButton("scan_all", "1) Scan all curves (flags only, changes nothing)",
                   class = "btn-primary", width = "100%"),
      helpText("Works out which curves are zero-growth or unreliable. It does NOT ",
               "change any window and does NOT exclude anything - it just computes ",
               "the flags so you can find them below."),
      selectInput("flt_status", "2) Show only:",
                  choices = c("All curves",
                              "ZERO GROWTH (flagged)",
                              "UNRELIABLE r (flagged)",
                              "OK",
                              "EXCLUDED - all",
                              "EXCLUDED - but NOT zero-growth",
                              "EXCLUDED - zero-growth only"),
                  selected = "All curves"),
      actionButton("refresh_filter", "Refresh list", class = "btn-default"),
      helpText("The EXCLUDED views let you review what you have binned. ",
               "'EXCLUDED - but NOT zero-growth' shows the ones you removed by ",
               "hand, so you can double-check them. Press 'Refresh list' after ",
               "ticking / unticking to update the view."),
      verbatimTextOutput("scan_summary"),
      selectInput("curve", "Curve:", choices = setNames(curves$key, curves$label)),
      fluidRow(column(6, actionButton("prev", "◀ Prev", width = "100%")),
               column(6, actionButton("nxt",  "Next ▶", width = "100%"))),
      br(),
      radioButtons("mode", "Click sets:", c("Start", "End"), selected = "Start", inline = TRUE),
      checkboxInput("show_guide", "Show model guide curve (blue)", value = TRUE),
      checkboxInput("show_exp", "Show suggested exponential phase (orange dotted)", value = TRUE),
      numericInput("exp_frac", "Suggested end at % of draw-down (from your start)",
                   value = 85, min = 50, max = 99.9, step = 1),
      actionButton("clear_one", "Clear this curve (back to auto)"),
      actionButton("clear_all", "Clear ALL windows (deactivate optimiser)"),
      helpText("Click sets Start first (resets to 'Start' on each new curve), then switch to End. Green = your start, Red = your end, orange dashed = suggested end (computed from your start, at the % below). It's a guide only — you can still select any window you like. The long low-O2 plateau is cropped from view."),
      hr(),
      strong("Auto-detect window"),
      fluidRow(
        column(3, numericInput("r2_target", "R2 >=", value = 0.99,
                               min = 0.5, max = 1, step = 0.005)),
        column(3, numericInput("rmse_max", "RMSE <=", value = 0.08,
                               min = 0, step = 0.01)),
        column(3, numericInput("min_pts", "Min pts", value = 50, min = 6, step = 5)),
        column(3, numericInput("start_search", "Start +min", value = 60,
                               min = 0, step = 15))
      ),
      radioButtons("start_mode", "Start rule:",
                   c("Peak (+ optional forward search)" = "peak",
                     "Fixed % drawdown from peak (consistent)" = "drawdown"),
                   selected = "peak"),
      conditionalPanel(
        condition = "input.start_mode == 'drawdown'",
        numericInput("start_drawdown", "Start drawdown %", value = 5,
                     min = 0, max = 50, step = 1)
      ),
      fluidRow(
        column(6, actionButton("auto_one", "Auto THIS curve", width = "100%")),
        column(6, actionButton("auto_all", "Auto ALL curves", width = "100%"))
      ),
      tags$br(),
      strong("Stable-r optimiser (recommended):"),
      helpText("Tries many windows and picks the one whose growth rate r is ",
               "STABLE (barely changes when the window is nudged) - i.e. the r ",
               "you can actually believe. Judged from this curve alone; it never ",
               "looks at the temperature curve."),
      fluidRow(
        column(6, actionButton("stab_one", "Stable-r THIS curve",
                               width = "100%", class = "btn-success")),
        column(6, actionButton("stab_all", "Stable-r ALL curves",
                               width = "100%", class = "btn-success"))
      ),
      verbatimTextOutput("stab_info"),
      fluidRow(
        column(12, actionButton("stab_iso", "Stable-r THIS ISOLATE (all its curves)",
                                width = "100%", class = "btn-success"))
      ),
      helpText("'Stable-r THIS ISOLATE' runs stable-r on every temperature / replicate ",
               "curve of the isolate chosen in the 'Isolate' filter above (or, if that ",
               "is 'All', the isolate of the curve currently on screen)."),
      helpText("Start = O2 peak; end = longest window keeping R2 >= target. You can still nudge afterwards."),
      hr(),
      strong("Trim ALL by time interval"),
      fluidRow(
        column(6, numericInput("all_start", "All start (min)", value = NA)),
        column(6, numericInput("all_end",   "All end (min)",   value = NA))
      ),
      actionButton("apply_all_time", "Apply interval to ALL curves"),
      helpText("Sets the same start/end (min) on every curve. Leave a box blank to keep that side unchanged."),
      checkboxInput("exclude_this",
                    "Don't include this sample (exclude from plots)",
                    value = FALSE),
      actionButton("include_all", "Include ALL points (clear every exclusion)"),
      helpText("Un-excludes every sample at once, so nothing is dropped from the ",
               "plots."),
      tags$br(),
      actionButton("restore_nogrowth", "Restore zero-growth curves (KEEP my other exclusions)",
                   class = "btn-info"),
      helpText("Un-excludes ONLY the straight-line / zero-growth curves, leaving ",
               "every other exclusion you made untouched. Use this if you binned ",
               "the zero-growth ones and want them back: their respiration is real, ",
               "and their CUE = 0 is a genuine result."),
      tags$br(),
      actionButton("excl_unreliable", "Exclude ALL unreliable-r curves",
                   class = "btn-warning"),
      helpText("Excludes the curves flagged UNRELIABLE r (no window gives a stable ",
               "growth rate). Needs a 'Stable-r ALL curves' scan first. Optional - ",
               "nothing is excluded unless you press this."),
      tags$br(),
      actionButton("excl_nogrowth", "Exclude ALL zero-growth curves",
                   class = "btn-warning"),
      helpText("Finds and excludes every curve with ZERO GROWTH - either no ",
               "curvature (a straight decreasing line) or no signal at all (the ",
               "culture consumed almost no O2 and is dead). Their r is unmeasurable, ",
               "so they are discarded from growth AND respiration. Nothing is ",
               "discarded automatically anywhere else - this button is the decision. ",
               "Use 'Include ALL points' to undo."),
      hr(),
      strong(textOutput("count")),
      div(style = "max-height:200px; overflow-y:auto;", tableOutput("tbl")),
      downloadButton("dl", "Download manual_fit_windows.csv"),
      hr(),
      strong(textOutput("excl_count")),
      strong("Windows + exclusions auto-save to tables/ (loaded by config.R when USE_APP_TRIM_FILES <- TRUE):"),
      verbatimTextOutput("combo"),
      downloadButton("excl_dl", "Download plot_exclude_points.csv")
    ),
    mainPanel(
      width = 8,
      plotOutput("plot", height = "640px", click = "plot_click"),
      textOutput("curinfo")
    )
  )
)

# ---- Server -----------------------------------------------------------------
server <- function(input, output, session) {
  wins <- reactiveVal(init)
  excl <- reactiveVal(excl_init)

  # Results of the last scan (status per curve). ADVICE ONLY - it never excludes
  # anything and never changes a window by itself.
  scan_res <- reactiveVal(NULL)

  # Compute the status of every curve. Changes NOTHING - no windows, no exclusions.
  observeEvent(input$scan_all, {
    n <- nrow(curves); st <- vector("list", n)
    withProgress(message = "Scanning all curves (nothing is being changed)...", value = 0, {
      for (i in seq_len(n)) {
        k <- curves$key[i]
        d <- long %>% dplyr::filter(key == k, is.finite(Time), is.finite(Oxygen)) %>%
          dplyr::arrange(Time)
        sw <- stable_window(d$Time, d$Oxygen)
        st[[i]] <- if (is.null(sw)) {
          data.frame(key = k, label = curves$label[i], status = "UNRELIABLE r",
                     r = NA_real_, rt = NA_real_, O2_used = NA_real_,
                     stability_cv = NA_real_, stringsAsFactors = FALSE)
        } else {
          data.frame(key = k, label = curves$label[i],
                     status = if (isTRUE(sw$no_growth)) "ZERO GROWTH" else "OK",
                     r = round(sw$r, 5), rt = round(sw$rt, 2),
                     O2_used = round(sw$drawdown, 2),
                     stability_cv = round(sw$stability_cv, 3), stringsAsFactors = FALSE)
        }
        incProgress(1 / n)
      }
    })
    scan_res(do.call(rbind, st))
    sr <- scan_res()
    showNotification(
      sprintf(paste0("Scan done (nothing changed). %d OK | %d ZERO GROWTH | %d UNRELIABLE r.\n",
                     "Use 'Show only:' to step through them."),
              sum(sr$status == "OK"), sum(sr$status == "ZERO GROWTH"),
              sum(sr$status == "UNRELIABLE r")),
      type = "message", duration = 12)
  })

  output$scan_summary <- renderText({
    sr <- scan_res()
    if (is.null(sr)) return("No scan yet - press 'Scan all curves' to compute the flags.")
    sprintf("Scan: %d OK | %d ZERO GROWTH | %d UNRELIABLE r   (nothing excluded)",
            sum(sr$status == "OK"), sum(sr$status == "ZERO GROWTH"),
            sum(sr$status == "UNRELIABLE r"))
  })

  # The curves currently VISIBLE in the dropdown (after all filters). Prev/Next
  # navigate this, not the full 1080.
  visible_keys <- reactiveVal(curves$key)

  # ---- Filter the curve list (1080 curves is unusable unfiltered) ------------
  observeEvent(list(input$flt_iso, input$flt_status, input$refresh_filter), {
    cc <- curves
    if (!is.null(input$flt_iso) && input$flt_iso != "All") {
      cc <- cc[cc$Taxon == input$flt_iso, , drop = FALSE]
    }
    fs <- input$flt_status
    if (!is.null(fs) && fs != "All curves") {
      ex <- excl()                       # what you have currently binned
      sr <- scan_res()                   # flags (NULL until you run a scan)

      needs_scan <- fs %in% c("ZERO GROWTH (flagged)", "UNRELIABLE r (flagged)", "OK",
                              "EXCLUDED - but NOT zero-growth",
                              "EXCLUDED - zero-growth only")
      if (needs_scan && is.null(sr)) {
        showNotification("Press 'Scan all curves' first - that computes the flags.",
                         type = "warning")
        return()
      }
      zg <- if (is.null(sr)) character(0) else sr$key[sr$status == "ZERO GROWTH"]

      keys <- switch(
        fs,
        "ZERO GROWTH (flagged)"          = zg,
        "UNRELIABLE r (flagged)"         = sr$key[sr$status == "UNRELIABLE r"],
        "OK"                             = sr$key[sr$status == "OK"],
        "EXCLUDED - all"                 = ex,
        "EXCLUDED - but NOT zero-growth" = setdiff(ex, zg),
        "EXCLUDED - zero-growth only"    = intersect(ex, zg),
        cc$key)

      if (length(keys) == 0) {
        showNotification(sprintf("Nothing matches '%s'.", fs), type = "warning")
        return()
      }
      cc <- cc[cc$key %in% keys, , drop = FALSE]
    }
    if (nrow(cc) == 0) {
      showNotification("No curves match that filter.", type = "warning")
      return()
    }
    # keep the current curve selected if it survives the filter
    sel <- if (!is.null(input$curve) && input$curve %in% cc$key) input$curve else cc$key[1]
    visible_keys(cc$key)          # Prev/Next now walk exactly this list
    updateSelectInput(session, "curve",
                      choices = setNames(cc$key, cc$label), selected = sel)
    showNotification(sprintf("%d curve(s) shown. Prev/Next step through these only.",
                             nrow(cc)), duration = 4)
  }, ignoreInit = TRUE)

  # On each new curve: reset the click mode to "Start" and sync the exclude box.
  observeEvent(input$curve, {
    updateRadioButtons(session, "mode", selected = "Start")
    updateCheckboxInput(session, "exclude_this", value = input$curve %in% excl())
  }, ignoreInit = FALSE)

  observeEvent(input$exclude_this, {
    k <- input$curve; cur <- excl()
    excl(if (isTRUE(input$exclude_this)) union(cur, k) else setdiff(cur, k))
  }, ignoreInit = TRUE)

  # Include ALL points: clear every exclusion at once (empties the exclusion set,
  # which auto-saves to plot_exclude_points.csv and unticks the current box).
  observeEvent(input$include_all, {
    excl(character(0))
    updateCheckboxInput(session, "exclude_this", value = FALSE)
  })

  get_row <- function(k) {
    w <- wins(); r <- w[w$key == k, , drop = FALSE]
    if (nrow(r) == 0) list(fit_start = NA_real_, fit_end = NA_real_)
    else list(fit_start = r$fit_start[1], fit_end = r$fit_end[1])
  }
  set_val <- function(k, side, val) {
    w <- wins()
    if (!(k %in% w$key)) w <- rbind(w, data.frame(key = k, fit_start = NA_real_, fit_end = NA_real_))
    w[w$key == k, side] <- val
    # drop fully-empty rows
    w <- w[!(is.na(w$fit_start) & is.na(w$fit_end)), , drop = FALSE]
    wins(w)
  }

  # Prev/Next must walk the list you can actually SEE (i.e. after the Isolate /
  # Temperature / Show-only filters), not all 1080 curves.
  observeEvent(input$prev, {
    ks <- visible_keys()
    i <- match(input$curve, ks)
    if (!is.na(i) && i > 1) updateSelectInput(session, "curve", selected = ks[i - 1])
  })
  observeEvent(input$nxt, {
    ks <- visible_keys()
    i <- match(input$curve, ks)
    if (!is.na(i) && i < length(ks)) updateSelectInput(session, "curve", selected = ks[i + 1])
  })

  observeEvent(input$plot_click, {
    k <- input$curve; x <- input$plot_click$x
    if (is.null(x)) return()
    set_val(k, if (input$mode == "Start") "fit_start" else "fit_end", round(x, 1))
  })
  observeEvent(input$clear_one, {
    w <- wins(); wins(w[w$key != input$curve, , drop = FALSE])
  })

  observeEvent(input$apply_all_time, {
    s <- suppressWarnings(as.numeric(input$all_start))
    e <- suppressWarnings(as.numeric(input$all_end))
    if (!is.finite(s) && !is.finite(e)) {
      showNotification("Enter a start and/or end time first.", type = "warning")
      return()
    }
    w <- wins()
    for (k in curves$key) {
      if (!(k %in% w$key))
        w <- rbind(w, data.frame(key = k, fit_start = NA_real_, fit_end = NA_real_))
      if (is.finite(s)) w[w$key == k, "fit_start"] <- s
      if (is.finite(e)) w[w$key == k, "fit_end"]   <- e
    }
    w <- w[!(is.na(w$fit_start) & is.na(w$fit_end)), , drop = FALSE]
    wins(w)
    showNotification("Applied time interval to all curves.", type = "message")
  })

  observeEvent(input$auto_one, {
    k <- input$curve
    d <- long %>% dplyr::filter(key == k, is.finite(Time), is.finite(Oxygen)) %>%
      dplyr::arrange(Time)
    sd_frac <- (if (is.null(input$start_drawdown)) 5 else input$start_drawdown) / 100
    se <- auto_detect_window(d$Time, d$Oxygen, input$r2_target, input$rmse_max,
                             input$min_pts, input$start_search,
                             start_mode = input$start_mode, start_drawdown_frac = sd_frac)
    if (all(is.finite(se))) {
      set_val(k, "fit_start", round(se[1], 1))
      set_val(k, "fit_end",   round(se[2], 1))
    } else {
      showNotification("Auto-detect could not find a valid window for this curve.",
                       type = "warning")
    }
  })

  observeEvent(input$auto_all, {
    n <- nrow(curves)
    w <- wins()
    sd_frac <- (if (is.null(input$start_drawdown)) 5 else input$start_drawdown) / 100
    withProgress(message = "Auto-detecting all curves...", value = 0, {
      for (i in seq_len(n)) {
        k <- curves$key[i]
        d <- long %>% dplyr::filter(key == k, is.finite(Time), is.finite(Oxygen)) %>%
          dplyr::arrange(Time)
        se <- auto_detect_window(d$Time, d$Oxygen, input$r2_target, input$rmse_max,
                                 input$min_pts, input$start_search,
                                 start_mode = input$start_mode, start_drawdown_frac = sd_frac)
        if (all(is.finite(se))) {
          if (!(k %in% w$key))
            w <- rbind(w, data.frame(key = k, fit_start = NA_real_, fit_end = NA_real_))
          w[w$key == k, "fit_start"] <- round(se[1], 1)
          w[w$key == k, "fit_end"]   <- round(se[2], 1)
        }
        incProgress(1 / n)
      }
    })
    w <- w[!(is.na(w$fit_start) & is.na(w$fit_end)), , drop = FALSE]
    wins(w)
    showNotification("Auto-detect ALL done. Review and tweak as needed.", type = "message")
  })

  # ---- STABLE-r optimiser (opt-in: nothing happens unless you click) ---------

  observeEvent(input$stab_one, {
    k <- input$curve
    d <- long %>% dplyr::filter(key == k, is.finite(Time), is.finite(Oxygen)) %>%
      dplyr::arrange(Time)
    sw <- stable_window(d$Time, d$Oxygen)
    if (is.null(sw)) {
      showNotification(paste0("No window gives a usable r for this curve - ",
                              "its growth rate is not a reliable measurement."),
                       type = "warning", duration = 8)
      return()
    }
    set_val(k, "fit_start", round(sw$start, 1))
    set_val(k, "fit_end",   round(sw$end, 1))
    if (isTRUE(sw$no_growth)) {
      showNotification(
        sprintf(paste0("NO GROWTH: this trace is a straight decreasing line ",
                       "(r x window = %.2f). Respiration is real, but r is not ",
                       "measurable. It will be dropped from the GROWTH fit ",
                       "automatically - do NOT tick 'don't include', that would ",
                       "also throw away its good respiration."), sw$rt),
        type = "warning", duration = 12)
    } else {
      showNotification(sprintf("Stable-r window set (r = %.5f, %d good windows, stability CV = %.2f).",
                               sw$r, sw$n_good, sw$stability_cv), type = "message")
    }
  })

  observeEvent(input$stab_all, {
    n <- nrow(curves); w <- wins(); n_ok <- 0L; n_bad <- 0L; n_ng <- 0L
    st <- vector("list", n)
    withProgress(message = "Stable-r: optimising all curves...", value = 0, {
      for (i in seq_len(n)) {
        k <- curves$key[i]
        d <- long %>% dplyr::filter(key == k, is.finite(Time), is.finite(Oxygen)) %>%
          dplyr::arrange(Time)
        sw <- stable_window(d$Time, d$Oxygen)
        if (!is.null(sw)) {
          if (!(k %in% w$key))
            w <- rbind(w, data.frame(key = k, fit_start = NA_real_, fit_end = NA_real_))
          w[w$key == k, "fit_start"] <- round(sw$start, 1)
          w[w$key == k, "fit_end"]   <- round(sw$end, 1)
          n_ok <- n_ok + 1L
          if (isTRUE(sw$no_growth)) n_ng <- n_ng + 1L
          st[[i]] <- data.frame(
            key = k, label = curves$label[i],
            status = if (isTRUE(sw$no_growth)) "ZERO GROWTH" else "OK",
            r = round(sw$r, 5), rt = round(sw$rt, 2),
            O2_used = round(sw$drawdown, 2),
            stability_cv = round(sw$stability_cv, 3),
            stringsAsFactors = FALSE)
        } else {
          n_bad <- n_bad + 1L
          st[[i]] <- data.frame(
            key = k, label = curves$label[i], status = "UNRELIABLE r",
            r = NA_real_, rt = NA_real_, O2_used = NA_real_, stability_cv = NA_real_,
            stringsAsFactors = FALSE)
        }
        incProgress(1 / n)
      }
    })
    w <- w[!(is.na(w$fit_start) & is.na(w$fit_end)), , drop = FALSE]
    wins(w)
    scan_res(do.call(rbind, st))
    showNotification(
      sprintf(paste0("Stable-r done. Windows set on %d curves.\n",
                     "FLAGGED (advice only - NOTHING has been excluded):\n",
                     "  %d zero growth   |   %d unreliable r\n",
                     "Use 'Show only:' to step through them and tick each yourself."),
              n_ok, n_ng, n_bad),
      type = "message", duration = 15)
  })

  # Stable-r for one isolate GROUP: every T/replicate curve of the isolate chosen
  # in the 'Isolate' filter (or the current curve's isolate when the filter = All).
  observeEvent(input$stab_iso, {
    target <- if (!is.null(input$flt_iso) && input$flt_iso != "All") {
      input$flt_iso
    } else {
      cur <- curves[curves$key == input$curve, , drop = FALSE]
      if (nrow(cur) == 0) NA_character_ else cur$Taxon[1]
    }
    if (is.na(target) || !nzchar(target)) {
      showNotification("Choose an isolate in the 'Isolate' filter first (or select a curve).",
                       type = "warning", duration = 8); return()
    }
    grp <- curves[curves$Taxon == target, , drop = FALSE]
    if (nrow(grp) == 0) {
      showNotification(sprintf("No curves found for isolate '%s'.", target),
                       type = "warning", duration = 8); return()
    }
    ng <- nrow(grp); w <- wins(); n_ok <- 0L; n_bad <- 0L; n_ng <- 0L
    st <- vector("list", ng)
    withProgress(message = sprintf("Stable-r: isolate %s (%d curves)...", target, ng), value = 0, {
      for (i in seq_len(ng)) {
        k <- grp$key[i]
        d <- long %>% dplyr::filter(key == k, is.finite(Time), is.finite(Oxygen)) %>%
          dplyr::arrange(Time)
        sw <- stable_window(d$Time, d$Oxygen)
        if (!is.null(sw)) {
          if (!(k %in% w$key))
            w <- rbind(w, data.frame(key = k, fit_start = NA_real_, fit_end = NA_real_))
          w[w$key == k, "fit_start"] <- round(sw$start, 1)
          w[w$key == k, "fit_end"]   <- round(sw$end, 1)
          n_ok <- n_ok + 1L
          if (isTRUE(sw$no_growth)) n_ng <- n_ng + 1L
          st[[i]] <- data.frame(
            key = k, label = grp$label[i],
            status = if (isTRUE(sw$no_growth)) "ZERO GROWTH" else "OK",
            r = round(sw$r, 5), rt = round(sw$rt, 2),
            O2_used = round(sw$drawdown, 2),
            stability_cv = round(sw$stability_cv, 3),
            stringsAsFactors = FALSE)
        } else {
          n_bad <- n_bad + 1L
          st[[i]] <- data.frame(
            key = k, label = grp$label[i], status = "UNRELIABLE r",
            r = NA_real_, rt = NA_real_, O2_used = NA_real_, stability_cv = NA_real_,
            stringsAsFactors = FALSE)
        }
        incProgress(1 / ng)
      }
    })
    w <- w[!(is.na(w$fit_start) & is.na(w$fit_end)), , drop = FALSE]
    wins(w)
    scan_res(do.call(rbind, st))
    showNotification(
      sprintf(paste0("Stable-r done for isolate %s: windows set on %d of %d curves.\n",
                     "FLAGGED (advice only - NOTHING excluded): %d zero growth | %d unreliable r."),
              target, n_ok, ng, n_ng, n_bad),
      type = "message", duration = 15)
  })

  # Live stability readout for the curve on screen (does not change anything).
  output$stab_info <- renderText({
    k <- input$curve
    req(k)
    d <- long %>% dplyr::filter(key == k, is.finite(Time), is.finite(Oxygen)) %>%
      dplyr::arrange(Time)
    sw <- stable_window(d$Time, d$Oxygen)
    if (is.null(sw)) return("Stable-r: NO reliable window for this curve (r not measurable).")
    verdict <- if (is.finite(sw$stability_cv) && sw$stability_cv < 0.10) {
      "r is STABLE -> trustworthy"
    } else if (is.finite(sw$stability_cv) && sw$stability_cv < 0.25) {
      "r is moderately stable"
    } else {
      "r SWINGS with the window -> do not trust this r"
    }
    curv <- if (isTRUE(sw$no_growth)) {
      why <- if (isTRUE(sw$no_signal) && isTRUE(sw$no_curv)) {
        sprintf("no curvature (r x window = %.2f) AND no signal (only %.2f mg/L O2 used)",
                sw$rt, sw$drawdown)
      } else if (isTRUE(sw$no_signal)) {
        sprintf("NO SIGNAL - the culture consumed only %.2f mg/L O2 (typical: ~4.7).\n    The cells are dead; this 'growth' is fitted noise",
                sw$drawdown)
      } else {
        sprintf("NO CURVATURE - straight decreasing line (r x window = %.2f <= %.2f).\n    Constant respiration, no growth",
                sw$rt, CURV_MIN_RT)
      }
      paste0("\n*** ZERO GROWTH: ", why, ".\n",
             "    r is NOT measurable here. Use 'Exclude ALL zero-growth curves'\n",
             "    to discard this and every other one like it.")
    } else {
      sprintf("\ncurvature: r x window = %.2f | O2 used = %.2f mg/L  -> real growth signal",
              sw$rt, sw$drawdown)
    }
    paste0(sprintf("Stable-r suggestion: %.1f - %.1f min | r = %.5f | R2 = %.3f\n%d good windows | stability CV = %.2f  (%s)",
                   sw$start, sw$end, sw$r, sw$r2, sw$n_good, sw$stability_cv, verdict),
           curv)
  })

  # Exclude every zero-growth (straight-line) curve in one go. NOTE: this removes
  # them from RESPIRATION too - see the warning next to the button.
  observeEvent(input$excl_nogrowth, {
    n <- nrow(curves); hits <- character(0)
    withProgress(message = "Finding zero-growth (straight-line) curves...", value = 0, {
      for (i in seq_len(n)) {
        k <- curves$key[i]
        d <- long %>% dplyr::filter(key == k, is.finite(Time), is.finite(Oxygen)) %>%
          dplyr::arrange(Time)
        sw <- stable_window(d$Time, d$Oxygen)
        if (!is.null(sw) && isTRUE(sw$no_growth)) hits <- c(hits, k)
        incProgress(1 / n)
      }
    })
    if (length(hits) == 0) {
      showNotification("No zero-growth curves found - nothing excluded.", type = "message")
      return()
    }
    excl(union(excl(), hits))
    updateCheckboxInput(session, "exclude_this", value = input$curve %in% excl())
    showNotification(
      sprintf(paste0("Excluded %d zero-growth curve(s) - they are now dropped from ",
                     "EVERYTHING, respiration included. Use 'Include ALL points' to undo."),
              length(hits)),
      type = "warning", duration = 12)
  })

  # Exclude every curve flagged UNRELIABLE r. Opt-in: only runs when you click.
  observeEvent(input$excl_unreliable, {
    sr <- scan_res()
    if (is.null(sr)) {
      showNotification("Run 'Stable-r ALL curves' first to compute the flags.",
                       type = "warning"); return()
    }
    hits <- sr$key[sr$status == "UNRELIABLE r"]
    if (!length(hits)) {
      showNotification("No curves are flagged UNRELIABLE r.", type = "message"); return()
    }
    excl(union(excl(), hits))
    updateCheckboxInput(session, "exclude_this", value = input$curve %in% excl())
    showNotification(sprintf("Excluded %d unreliable-r curve(s). 'Include ALL points' undoes it.",
                             length(hits)), type = "warning", duration = 10)
  })

  # Surgically un-exclude ONLY the zero-growth curves, keeping every other
  # exclusion the user made by hand.
  observeEvent(input$restore_nogrowth, {
    cur <- excl()
    if (!length(cur)) {
      showNotification("Nothing is excluded.", type = "message"); return()
    }
    hits <- character(0)
    withProgress(message = "Checking excluded curves...", value = 0, {
      for (i in seq_along(cur)) {
        k <- cur[i]
        d <- long %>% dplyr::filter(key == k, is.finite(Time), is.finite(Oxygen)) %>%
          dplyr::arrange(Time)
        sw <- stable_window(d$Time, d$Oxygen)
        if (!is.null(sw) && isTRUE(sw$no_growth)) hits <- c(hits, k)
        incProgress(1 / length(cur))
      }
    })
    if (!length(hits)) {
      showNotification("None of your exclusions are zero-growth curves - nothing changed.",
                       type = "message"); return()
    }
    excl(setdiff(cur, hits))
    updateCheckboxInput(session, "exclude_this", value = input$curve %in% excl())
    showNotification(
      sprintf(paste0("Restored %d zero-growth curve(s). Your other %d exclusion(s) ",
                     "were kept."), length(hits), length(setdiff(cur, hits))),
      type = "message", duration = 10)
  })

  # Deactivate / revert: drop every manual window, back to automatic trimming.
  observeEvent(input$clear_all, {
    wins(data.frame(key = character(0), fit_start = numeric(0), fit_end = numeric(0),
                    stringsAsFactors = FALSE))
    showNotification("Cleared ALL manual windows - back to automatic trimming.",
                     type = "message")
  })

  # Fit resp_model to the CURRENT window and report how well the raw data matches
  # it (R2 and RMSE). Shared by the plot guide and the readout below the plot.
  winfit <- reactive({
    k <- input$curve
    d <- long %>% dplyr::filter(key == k, is.finite(Time), is.finite(Oxygen)) %>%
      dplyr::arrange(Time)
    a <- meta[meta$key == k, , drop = FALSE]
    r <- get_row(k)
    auto_start <- if (nrow(a)) a$auto_start[1] else NA_real_
    auto_end   <- if (nrow(a)) a$auto_end[1]   else NA_real_
    eff_start  <- if (is.finite(r$fit_start)) r$fit_start else auto_start
    eff_end    <- if (is.finite(r$fit_end))   r$fit_end   else auto_end
    res <- list(ok = FALSE, eff_start = eff_start, eff_end = eff_end,
                co = NULL, rmse = NA_real_, r2 = NA_real_, n = 0L)
    if (has_minpack && is.finite(eff_start) && is.finite(eff_end) && eff_end > eff_start) {
      dw <- d[d$Time >= eff_start & d$Time <= eff_end, , drop = FALSE]
      if (nrow(dw) >= 6) {
        t0     <- dw$Time - min(dw$Time)
        slope0 <- suppressWarnings(stats::median(diff(dw$Oxygen) / diff(t0), na.rm = TRUE))
        K0     <- if (is.finite(slope0)) max(abs(slope0), 1e-6) else 1e-3
        ft <- try(minpack.lm::nlsLM(
          Oxygen ~ resp_model(r, K, Time0, O2_0),
          data  = data.frame(Time0 = t0, Oxygen = dw$Oxygen),
          start = list(r = 1e-3, K = K0, O2_0 = dw$Oxygen[1]),
          lower = c(r = 1e-6, K = 1e-10, O2_0 = min(dw$Oxygen) - 1),
          upper = c(r = 0.15, K = 1,     O2_0 = max(dw$Oxygen) + 1),
          control = minpack.lm::nls.lm.control(maxiter = 200)
        ), silent = TRUE)
        if (!inherits(ft, "try-error")) {
          co    <- coef(ft)
          preds <- resp_model(co[["r"]], co[["K"]], t0, co[["O2_0"]])
          resid <- dw$Oxygen - preds
          sstot <- sum((dw$Oxygen - mean(dw$Oxygen))^2)
          res$ok   <- TRUE
          res$co   <- co
          res$rmse <- sqrt(mean(resid^2))
          res$r2   <- if (sstot > 1e-12) 1 - sum(resid^2) / sstot else NA_real_
          res$n    <- nrow(dw)
        }
      }
    }
    res
  })

  output$plot <- renderPlot({
    k <- input$curve
    d <- long %>% dplyr::filter(key == k, is.finite(Time), is.finite(Oxygen)) %>%
      dplyr::arrange(Time)
    a <- meta[meta$key == k, , drop = FALSE]
    r <- get_row(k)
    is_excl <- k %in% excl()

    auto_start <- if (nrow(a)) a$auto_start[1] else NA_real_
    auto_end   <- if (nrow(a)) a$auto_end[1]   else NA_real_
    start_ref  <- if (is.finite(r$fit_start)) r$fit_start else auto_start
    eff_start  <- start_ref
    eff_end    <- if (is.finite(r$fit_end)) r$fit_end else auto_end

    # Suggested END — recomputed LIVE from wherever Start currently is (your
    # clicked start if set, else the peak), at the % below. Moves on every click.
    frac     <- (if (is.null(input$exp_frac) || !is.finite(input$exp_frac)) 85 else input$exp_frac) / 100
    peak_t   <- exp_phase_bounds(d$Time, d$Oxygen)[1]
    # Anchor for the suggested end: your clicked start if set, else the automatic
    # (green-line) start from 02, else the O2 peak. So it shows sensibly BEFORE
    # you click a start, then tracks your start once you do.
    anchor_s <- if (is.finite(r$fit_start)) r$fit_start else
                if (is.finite(auto_start))  auto_start  else peak_t
    sugg_end <- if (isTRUE(input$show_exp)) suggest_end_95(d$Time, d$Oxygen, anchor_s, frac) else NA_real_

    # Crop the long low-O2 plateau, but keep room past the end so the projected
    # guide curve AND the suggested-end guide stay on screen.
    o2min <- min(d$Oxygen); o2max <- max(d$Oxygen); yrng <- o2max - o2min
    thr   <- o2min + 0.03 * yrng
    below <- d$Time[d$Oxygen <= thr]
    t_bot <- if (length(below)) min(below) else max(d$Time)
    x_hi  <- min(max(d$Time),
                 max(t_bot + 60, eff_start + 60, auto_end + 60, eff_end + 90,
                     if (is.finite(sugg_end)) sugg_end + 90 else -Inf, na.rm = TRUE),
                 na.rm = TRUE)

    # Guide curve: use the shared window fit (winfit). Draw it SOLID within the
    # window and DASHED beyond it (back- and forward-projection).
    wf <- winfit()
    guide <- NULL
    if (isTRUE(input$show_guide) && isTRUE(wf$ok)) {
      co <- wf$co
      gx <- seq(min(d$Time), x_hi, length.out = 400)
      gy <- resp_model(co[["r"]], co[["K"]], gx - eff_start, co[["O2_0"]])
      guide <- data.frame(Time = gx, Oxygen = gy)
      guide <- guide[is.finite(guide$Oxygen) &
                       guide$Oxygen >= o2min - 0.12 * yrng &
                       guide$Oxygen <= o2max + 0.12 * yrng, , drop = FALSE]
    }

    p <- ggplot(d, aes(Time, Oxygen)) +
      geom_point(size = 1.1, alpha = 0.7, colour = if (is_excl) "grey65" else "black") +
      labs(title = curves$label[match(k, curves$key)],
           subtitle = if (is_excl) "EXCLUDED - not in plots" else NULL,
           x = "Time (min)", y = expression("O"[2]*" (mg/L)")) +
      coord_cartesian(xlim = c(min(d$Time), x_hi),
                      ylim = c(o2min - 0.12 * yrng, o2max + 0.05 * yrng)) +
      theme_classic(13) +
      theme(plot.subtitle = ggplot2::element_text(colour = "red", face = "bold"))
    if (!is.null(guide) && nrow(guide) > 1) {
      g_in   <- guide[guide$Time >= eff_start & guide$Time <= eff_end, , drop = FALSE]  # solid
      g_pre  <- guide[guide$Time <= eff_start, , drop = FALSE]   # back-projection (dashed)
      g_post <- guide[guide$Time >= eff_end,   , drop = FALSE]   # forward projection (dashed)
      if (nrow(g_pre) > 1)
        p <- p + geom_line(data = g_pre,  aes(Time, Oxygen), colour = "blue",
                           linewidth = 0.8, linetype = "dashed")
      if (nrow(g_post) > 1)
        p <- p + geom_line(data = g_post, aes(Time, Oxygen), colour = "blue",
                           linewidth = 0.8, linetype = "dashed")
      if (nrow(g_in) > 1)
        p <- p + geom_line(data = g_in,   aes(Time, Oxygen), colour = "blue", linewidth = 1)
    }
    # Only guide drawn: the suggested END (orange), live from the current Start,
    # with a label marking it.
    if (is.finite(sugg_end)) {
      p <- p +
        geom_vline(xintercept = sugg_end, colour = "darkorange",
                   linetype = "dashed", linewidth = 1.1) +
        annotate("label", x = sugg_end, y = o2max, hjust = 0.5, vjust = 1,
                 label = sprintf("%g%% end\n%.0f min", frac * 100, sugg_end),
                 colour = "darkorange", fill = "white", label.size = 0, size = 3.1)
    }
    if (is.finite(r$fit_start))
      p <- p + geom_vline(xintercept = r$fit_start, colour = "forestgreen", linewidth = 1.1)
    if (is.finite(r$fit_end))
      p <- p + geom_vline(xintercept = r$fit_end, colour = "red", linewidth = 1.1)
    p
  })

  output$curinfo <- renderText({
    r  <- get_row(input$curve)
    wf <- winfit()
    qual <- if (isTRUE(wf$ok)) {
      sprintf("FIT MATCH:  R2 = %.4f   |   RMSE = %.4f   (n = %d points)",
              wf$r2, wf$rmse, wf$n)
    } else {
      "FIT MATCH:  (need >= 6 points between start and end)"
    }
    # Suggested end at frac% of draw-down measured from the anchor start.
    d <- long %>% dplyr::filter(key == input$curve, is.finite(Time), is.finite(Oxygen))
    a  <- meta[meta$key == input$curve, , drop = FALSE]
    auto_start <- if (nrow(a)) a$auto_start[1] else NA_real_
    frac   <- (if (is.null(input$exp_frac) || !is.finite(input$exp_frac)) 85 else input$exp_frac) / 100
    peak_t <- exp_phase_bounds(d$Time, d$Oxygen)[1]
    anchor <- if (is.finite(r$fit_start)) r$fit_start else
              if (is.finite(auto_start))  auto_start  else peak_t
    anchor_lbl <- if (is.finite(r$fit_start)) "your start" else
                  if (is.finite(auto_start))  "auto start" else "peak"
    sugg   <- suggest_end_95(d$Time, d$Oxygen, anchor, frac)
    sugg_txt <- if (is.finite(sugg))
      sprintf("SUGGESTED END (%g%% down from %s): %.1f min",
              frac * 100, anchor_lbl, sugg)
    else "SUGGESTED END: (n/a)"
    sprintf("Start: %s    End: %s    (blank = automatic)\n%s\n%s",
            ifelse(is.finite(r$fit_start), r$fit_start, "auto"),
            ifelse(is.finite(r$fit_end),   r$fit_end,   "auto"),
            sugg_txt, qual)
  })

  win_df <- reactive({
    w <- wins()
    if (nrow(w) == 0) return(curves[0, c("Taxon", "Replicate")] %>%
                               dplyr::mutate(fit_start = numeric(0), fit_end = numeric(0)))
    curves %>% dplyr::inner_join(w, by = "key") %>%
      dplyr::arrange(Taxon, Replicate) %>%
      dplyr::select(Taxon, Replicate, fit_start, fit_end)
  })

  output$count <- renderText(sprintf("Manual windows set: %d curve(s)", nrow(win_df())))
  output$tbl <- renderTable(win_df(), digits = 1, na = "auto")

  observe({ tryCatch(readr::write_csv(win_df(), out_csv), error = function(e) NULL) })

  # Discarded samples (excluded from plots). Shared with PLOT_EXCLUDE_POINTS.
  excl_df <- reactive({
    keys <- excl()
    if (length(keys) == 0) return(curves[0, c("Taxon", "Replicate")])
    curves %>% dplyr::filter(key %in% keys) %>%
      dplyr::distinct(Taxon, Replicate) %>%
      dplyr::arrange(Taxon, Replicate)
  })
  output$excl_count <- renderText(sprintf("Excluded samples: %d", nrow(excl_df())))
  observe({ tryCatch(readr::write_csv(excl_df(), excl_csv), error = function(e) NULL) })
  output$excl_dl <- downloadHandler(
    filename = function() "plot_exclude_points.csv",
    content = function(file) readr::write_csv(excl_df(), file)
  )

  output$combo <- renderText({
    w <- win_df(); e <- excl_df()
    fmt <- function(v) ifelse(is.na(v), "NA", format(v, trim = TRUE))

    # PLOT_EXCLUDE_POINTS block - same shape as the data.frame in config.R.
    if (nrow(e) == 0) {
      excl_txt <- paste0("PLOT_EXCLUDE_POINTS <- data.frame(\n",
        "  Taxon = character(0), Replicate = character(0),\n",
        "  stringsAsFactors = FALSE\n)")
    } else {
      excl_txt <- paste0(
        "PLOT_EXCLUDE_POINTS <- data.frame(\n",
        "  Taxon     = c(", paste(sprintf('"%s"', e$Taxon), collapse = ", "), "),\n",
        "  Replicate = c(", paste(sprintf('"%s"', e$Replicate), collapse = ", "), "),\n",
        "  stringsAsFactors = FALSE\n)")
    }

    if (nrow(w) == 0) {
      win_txt <- paste0("MANUAL_FIT_WINDOWS <- data.frame(\n",
        "  Taxon = character(0), Replicate = character(0),\n",
        "  fit_start = numeric(0), fit_end = numeric(0),\n",
        "  stringsAsFactors = FALSE\n)")
    } else {
      win_txt <- paste0(
        "MANUAL_FIT_WINDOWS <- data.frame(\n",
        "  Taxon     = c(", paste(sprintf('"%s"', w$Taxon), collapse = ", "), "),\n",
        "  Replicate = c(", paste(sprintf('"%s"', w$Replicate), collapse = ", "), "),\n",
        "  fit_start = c(", paste(vapply(w$fit_start, fmt, ""), collapse = ", "), "),\n",
        "  fit_end   = c(", paste(vapply(w$fit_end,   fmt, ""), collapse = ", "), "),\n",
        "  stringsAsFactors = FALSE\n)")
    }

    paste0("# Windows + exclusions also auto-save to tables/ and are loaded by\n",
           "# config.R when USE_APP_TRIM_FILES <- TRUE. Paste below only if you\n",
           "# prefer to hard-code them into config.R:\n",
           excl_txt, "\n\n", win_txt)
  })

  output$dl <- downloadHandler(
    filename = function() "manual_fit_windows.csv",
    content = function(file) readr::write_csv(win_df(), file)
  )
}

# ---- Launch ------------------------------------------------------------------
# HEADLESS GUARD: this file must be source()-able in an unattended run
# (scripts/run_all.sh sources it to prove it still parses and that its inputs
# exist) without opening a browser or blocking. The app only starts in an
# interactive session, or when OXYMODEL_LAUNCH_APPS=1 is set explicitly.
.launch_app <- interactive() ||
  identical(Sys.getenv("OXYMODEL_LAUNCH_APPS"), "1")

if (.launch_app) {
  shinyApp(ui, server)
} else {
  message("03_trim_selector: non-interactive session - app NOT launched.")
  message("  Its committed outputs are treated as INPUTS to the pipeline:")
  message("    ", out_csv)
  message("    ", excl_csv)
  message("  To edit them, open this file in RStudio and 'Run App', or run")
  message("    OXYMODEL_LAUNCH_APPS=1 Rscript scripts/03_trim_selector.R")
  invisible(NULL)
}
