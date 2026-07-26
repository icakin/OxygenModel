# =============================================================================
# 02_trimming.R - Spline-shape trimming (hybrid start/end) + manual overrides
# =============================================================================
# For each (Taxon, Replicate) series: smooth with a spline, find the oxygen tip
# (start) and a deterministic end (whichever comes first of a time cap, the
# slope-flattening onset, a drop-fraction cap, or end of recording), and keep
# that window. Records per-curve metadata used downstream: the trim start time
# becomes delta_Ninoc_to_N0_min (the inoculation -> O2-start delay) for N0
# reconstruction in 05. Manual overrides (per curve code) take precedence.
#
#   Input : results/tables/Oxygen_All_Long.csv                (from 01)
#   Output: results/tables/Oxygen_Data_Smoothed_Trimmed.csv
#           results/tables/Oxygen_Data_Filtered.csv
#           results/tables/Skipped_Series_Log.csv
#           results/tables/Oxygen_Trimmed_Series_Metadata.csv   (-> 03 app, 05)
#           results/tables/Oxygen_Curve_Code_Key.csv
#   Figure: results/figures/oxygen_trimming_diagnostics.pdf     (multi-page)
# =============================================================================

# ---- locate scripts/ dir and source shared config ---------------------------
.this_dir <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  d <- if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE)) else
         tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (is.null(d) || is.na(d) || !nzchar(d)) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable() &&
        nzchar(rstudioapi::getActiveDocumentContext()$path))
      d <- dirname(rstudioapi::getActiveDocumentContext()$path) else d <- getwd()
  }
  d
})
source(file.path(.this_dir, "config.R"))
suppressPackageStartupMessages({
  library(zoo); library(stringr); library(scales); library(grid)
})

# ===== Trimming knobs (tune to your recording time-scale) ====================
TRIM_SPAR       <- 0.40    # smoothing-spline stiffness
MAX_TIME_MIN    <- Inf     # hard cap: drop Time > this before trimming. Inf = no cap.
MIN_ROWS_SERIES <- 8
MIN_RANGE_O2    <- 0.05    # skip a series whose smoothed O2 range is below this

# Start = argmax of the smoothed spline (the oxygen tip) within this window.
# Default [0, Inf] = global spline peak (matches "peak" onset). Narrow it if the
# rising limb or a late artefact is being picked as the tip.
START_SEARCH_MIN_TIME <- 0
START_SEARCH_MAX_TIME <- Inf

# Descending-run detection (used for run-start/end metadata + the main-run guard)
NEG_SLOPE_EPS     <- 0.00015
MIN_RUN_POINTS    <- 6
MAX_GAP_POINTS    <- 2
SEARCH_END_FRAC   <- 0.98
START_SMOOTH_K    <- 5

# End = whichever comes first of: start + END_MAX_WINDOW_MIN; the slope-recovery
# (flattening) onset; the point O2 has fallen by END_DROP_FRAC of its post-tip
# drop; end of recording.
END_MAX_WINDOW_MIN     <- Inf
END_DROP_FRAC          <- 0.85
END_SLOPE_RECOVER_FRAC <- 0.6
END_SMOOTH_K           <- 5

# Guards: reject a kept window shorter than these.
MIN_FIT_DURATION_MIN <- 6
MIN_FIT_POINTS       <- 8

# Diagnostic-plot zoom around the kept window.
PLOT_PAD_FRAC    <- 0.35
PLOT_PAD_MIN_MIN <- 5
X_MAJOR_BREAKS_N <- 12
Y_MAJOR_BREAKS_N <- 10

# ===== Manual override section ===============================================
# Run once, inspect results/figures/oxygen_trimming_diagnostics.pdf and
# results/tables/Oxygen_Curve_Code_Key.csv, then fill only the side(s) to change.
# (The 03 trim-selector app is the interactive alternative; it feeds 05 directly.)
manual_trim_overrides <- tibble::tibble(
  curve_code = character(),
  keep_start = numeric(),
  keep_end   = numeric()
)

# ===== Helpers ===============================================================
safe_spline_fit <- function(time, oxygen, spar = TRIM_SPAR) {
  fit <- tryCatch(smooth.spline(time, oxygen, spar = spar), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, length(time)))
  tryCatch(predict(fit, x = time)$y, error = function(e) rep(NA_real_, length(time)))
}
get_slopes <- function(time, y) { s <- diff(y) / diff(time); s[!is.finite(s)] <- NA_real_; s }
rle_to_runs <- function(flag_vec) {
  r <- rle(flag_vec); ends <- cumsum(r$lengths); starts <- c(1, head(ends, -1) + 1)
  tibble::tibble(value = r$values, start = starts, end = ends, len = r$lengths)
}
merge_descending_runs <- function(run_tbl, max_gap = MAX_GAP_POINTS) {
  if (nrow(run_tbl) == 0) return(run_tbl)
  desc <- run_tbl %>% dplyr::filter(value)
  if (nrow(desc) <= 1) return(desc %>% dplyr::select(start, end, len))
  out <- list(); cur_s <- desc$start[1]; cur_e <- desc$end[1]
  for (i in 2:nrow(desc)) {
    if (desc$start[i] - cur_e - 1 <= max_gap) { cur_e <- desc$end[i] } else {
      out[[length(out) + 1]] <- tibble::tibble(start = cur_s, end = cur_e)
      cur_s <- desc$start[i]; cur_e <- desc$end[i]
    }
  }
  out[[length(out) + 1]] <- tibble::tibble(start = cur_s, end = cur_e)
  dplyr::bind_rows(out) %>% dplyr::mutate(len = end - start + 1L)
}
score_descending_runs <- function(time, o2_fit, slopes, run_tbl) {
  if (nrow(run_tbl) == 0) return(run_tbl)
  run_tbl %>% dplyr::rowwise() %>%
    dplyr::mutate(
      time_start = time[start], time_end = time[min(end + 1L, length(time))],
      duration = time_end - time_start,
      y_start = o2_fit[start], y_end = o2_fit[min(end + 1L, length(o2_fit))],
      total_drop = y_start - y_end, mean_slope = mean(slopes[start:end], na.rm = TRUE),
      score = pmax(duration, 0) * pmax(total_drop, 0)) %>%
    dplyr::ungroup()
}
find_main_descending_run <- function(time, o2_fit) {
  slopes <- get_slopes(time, o2_fit)
  if (length(slopes) < 3 || all(!is.finite(slopes))) return(list(run = NULL, slopes = slopes, all_runs = NULL))
  peak_win <- which(is.finite(o2_fit) & is.finite(time) &
                      time >= START_SEARCH_MIN_TIME & time <= START_SEARCH_MAX_TIME)
  idx_peak <- if (length(peak_win) > 0) peak_win[which.max(o2_fit[peak_win])] else max(1L, floor(length(slopes) * 0.08))
  lo <- max(1L, min(idx_peak, length(slopes)))
  hi <- min(length(slopes), ceiling(length(slopes) * SEARCH_END_FRAC))
  neg_flag <- rep(FALSE, length(slopes)); idx_search <- lo:hi
  neg_flag[idx_search] <- is.finite(slopes[idx_search]) & slopes[idx_search] < -NEG_SLOPE_EPS
  run_tbl <- merge_descending_runs(rle_to_runs(neg_flag), max_gap = MAX_GAP_POINTS)
  if (nrow(run_tbl) == 0) return(list(run = NULL, slopes = slopes, all_runs = NULL))
  run_tbl <- run_tbl %>% dplyr::filter(len >= MIN_RUN_POINTS)
  if (nrow(run_tbl) == 0) return(list(run = NULL, slopes = slopes, all_runs = NULL))
  scored <- score_descending_runs(time, o2_fit, slopes, run_tbl) %>%
    dplyr::arrange(dplyr::desc(score), dplyr::desc(duration), dplyr::desc(total_drop))
  list(run = scored[1, ], slopes = slopes, all_runs = scored)
}
find_start_simple <- function(time, o2_fit, t_min = START_SEARCH_MIN_TIME, t_max = START_SEARCH_MAX_TIME) {
  n <- length(o2_fit); if (n == 0) return(NA_integer_)
  in_win <- which(is.finite(o2_fit) & is.finite(time) & time >= t_min & time <= t_max)
  if (length(in_win) == 0) {
    first_ok <- which(is.finite(o2_fit) & is.finite(time))
    if (length(first_ok) == 0) return(NA_integer_)
    return(as.integer(first_ok[1]))
  }
  as.integer(in_win[which.max(o2_fit[in_win])])
}
find_end_hybrid <- function(time, o2_fit, idx_start, max_window_min = END_MAX_WINDOW_MIN,
                            drop_frac = END_DROP_FRAC, recover_frac = END_SLOPE_RECOVER_FRAC,
                            smooth_k = END_SMOOTH_K) {
  n <- length(o2_fit)
  if (!is.finite(idx_start) || idx_start < 1L || idx_start >= n)
    return(list(idx_end = NA_integer_, end_method = "invalid_start"))
  after_idx <- (idx_start + 1L):n
  after_idx <- after_idx[is.finite(o2_fit[after_idx]) & is.finite(time[after_idx])]
  if (length(after_idx) == 0) return(list(idx_end = NA_integer_, end_method = "no_points_after_start"))
  y_start <- o2_fit[idx_start]; y_min <- min(o2_fit[after_idx], na.rm = TRUE)
  total_drop <- y_start - y_min; start_time <- time[idx_start]
  cap_time <- start_time + max_window_min
  idx_cap_c <- after_idx[time[after_idx] <= cap_time]
  idx_cap <- if (length(idx_cap_c) > 0) max(idx_cap_c) else after_idx[1]
  idx_flat <- NA_integer_
  slopes <- diff(o2_fit) / diff(time); slopes[!is.finite(slopes)] <- NA_real_
  if (length(slopes) >= smooth_k) {
    sm <- zoo::rollmedian(slopes, k = smooth_k, fill = NA, align = "center")
    sm[is.na(sm)] <- slopes[is.na(sm)]; slopes <- sm
  }
  slope_idx <- idx_start:(n - 1L); slope_idx <- slope_idx[is.finite(slopes[slope_idx])]
  if (length(slope_idx) > 0) {
    idx_steep <- slope_idx[which.min(slopes[slope_idx])]; min_slope <- slopes[idx_steep]
    if (is.finite(min_slope) && min_slope < 0) {
      thresh <- recover_frac * min_slope
      after_steep <- (idx_steep + 1L):(n - 1L); after_steep <- after_steep[is.finite(slopes[after_steep])]
      rec <- after_steep[slopes[after_steep] >= thresh]
      if (length(rec) > 0) idx_flat <- as.integer(min(rec[1] + 1L, n))
    }
  }
  idx_drop <- NA_integer_
  if (is.finite(total_drop) && total_drop > 0) {
    target_y <- y_start - drop_frac * total_drop
    hit <- after_idx[o2_fit[after_idx] <= target_y]
    if (length(hit) > 0) idx_drop <- as.integer(hit[1])
  }
  idx_eor <- as.integer(n)
  candidates <- c(idx_cap = idx_cap, idx_flat = idx_flat, idx_drop = idx_drop, idx_eor = idx_eor)
  candidates <- candidates[is.finite(candidates)]
  idx_end <- as.integer(min(candidates))
  end_method <- switch(names(which.min(candidates))[1],
                       idx_cap = "hybrid_time_cap", idx_flat = "hybrid_slope_flatten",
                       idx_drop = "hybrid_drop_fraction", idx_eor = "hybrid_end_of_recording",
                       "hybrid_unknown")
  list(idx_end = idx_end, end_method = end_method)
}

trim_one_series <- function(df) {
  df <- df %>% dplyr::arrange(Time)
  if (nrow(df) < MIN_ROWS_SERIES) return(list(ok = FALSE, reason = "Too few rows"))
  o2_fit <- safe_spline_fit(df$Time, df$Oxygen, spar = TRIM_SPAR)
  if (all(!is.finite(o2_fit))) return(list(ok = FALSE, reason = "Spline failed"))
  if ((max(o2_fit, na.rm = TRUE) - min(o2_fit, na.rm = TRUE)) < MIN_RANGE_O2)
    return(list(ok = FALSE, reason = "Too flat"))
  df <- df %>% dplyr::mutate(O2_fit = o2_fit)

  run_out <- find_main_descending_run(df$Time, df$O2_fit)
  main_run <- run_out$run
  if (is.null(main_run) || nrow(main_run) == 0) return(list(ok = FALSE, reason = "No main descending run found"))
  idx_run_start <- as.integer(main_run$start[1]); idx_run_end <- as.integer(main_run$end[1])

  idx_peak <- find_start_simple(df$Time, df$O2_fit)
  if (!is.finite(idx_peak) || is.na(idx_peak)) idx_peak <- idx_run_start
  idx_peak_raw <- idx_peak

  end_out <- find_end_hybrid(df$Time, df$O2_fit, idx_start = idx_peak)
  idx_end <- end_out$idx_end; end_reason <- end_out$end_method
  if (!is.finite(idx_end) || is.na(idx_end)) { idx_end <- min(nrow(df), idx_run_end + 1L); end_reason <- "fallback_main_run_end" }

  idx_steepest <- {
    sl <- diff(df$O2_fit) / diff(df$Time); sl[!is.finite(sl)] <- NA_real_
    after <- seq(idx_peak, length(sl)); after <- after[is.finite(sl[after])]
    if (length(after)) as.integer(after[which.min(sl[after])]) else NA_integer_
  }

  idx_peak <- max(1L, min(idx_peak, nrow(df)))
  idx_end  <- max(idx_peak + 1L, min(idx_end, nrow(df)))
  if (!is.finite(idx_steepest) || is.na(idx_steepest))
    idx_steepest <- max(idx_peak, min(nrow(df) - 1L, idx_run_start))
  if (idx_end <= idx_peak) return(list(ok = FALSE, reason = "Invalid peak/end order"))

  n_kept <- idx_end - idx_peak + 1L; dur_kept <- df$Time[idx_end] - df$Time[idx_peak]
  if (n_kept < MIN_FIT_POINTS)
    return(list(ok = FALSE, reason = sprintf("Kept window has %d points (< %d)", n_kept, MIN_FIT_POINTS)))
  if (!is.finite(dur_kept) || dur_kept < MIN_FIT_DURATION_MIN)
    return(list(ok = FALSE, reason = sprintf("Kept window %.2f min (< %.2f)", dur_kept, MIN_FIT_DURATION_MIN)))

  trimmed <- df[idx_peak:idx_end, ] %>%
    dplyr::mutate(curve_code = df$curve_code[1], series_id = df$series_id[1],
                  peak_time = df$Time[idx_peak], end_time_chosen = df$Time[idx_end],
                  end_reason = end_reason)

  meta <- tibble::tibble(
    curve_code = df$curve_code[1], series_id = df$series_id[1],
    Taxon = df$Taxon[1], Replicate = df$Replicate[1],
    n_points_total = nrow(df), n_points_trimmed = nrow(trimmed),
    peak_time = df$Time[idx_peak],
    main_run_start_time = df$Time[idx_run_start],
    main_run_end_time = df$Time[min(idx_run_end + 1L, nrow(df))],
    steepest_drop_time = df$Time[min(idx_steepest + 1L, nrow(df))],
    chosen_end_time = df$Time[idx_end],
    fit_duration_min = df$Time[idx_end] - df$Time[idx_peak],
    delta_Ninoc_to_N0_min = df$Time[idx_peak],   # inoculation -> O2-start delay
    end_reason = end_reason,
    main_run_duration = main_run$duration[1], main_run_drop = main_run$total_drop[1],
    main_run_score = main_run$score[1])

  list(ok = TRUE, reason = NA_character_, data = trimmed, meta = meta,
       idx_peak_raw = idx_peak_raw, idx_peak = idx_peak, idx_run_start = idx_run_start,
       idx_run_end = idx_run_end, idx_steepest = idx_steepest, idx_end = idx_end,
       all_runs = run_out$all_runs, full_df = df, manual_override = FALSE,
       manual_keep_start = NA_real_, manual_keep_end = NA_real_)
}

# ===== Load data =============================================================
if (!file.exists(LONG_CSV)) stop("Input file not found: ", LONG_CSV)
raw_before_cut <- readr::read_csv(LONG_CSV, show_col_types = FALSE) %>%
  dplyr::mutate(Taxon = as.character(Taxon), Replicate = as.character(Replicate),
                Time = as.numeric(Time), Oxygen = as.numeric(Oxygen),
                series_id = paste0(Taxon, " | Rep=", Replicate))
missing_cols <- setdiff(c("Time", "Taxon", "Replicate", "Oxygen"), names(raw_before_cut))
if (length(missing_cols) > 0) stop("Missing required columns in Oxygen_All_Long.csv: ",
                                   paste(missing_cols, collapse = ", "))
raw <- raw_before_cut %>% dplyr::filter(Time <= MAX_TIME_MIN)
if (nrow(raw) == 0) stop("No rows remain after MAX_TIME_MIN = ", MAX_TIME_MIN)

# ===== Curve codes ===========================================================
curve_key <- raw %>%
  dplyr::distinct(series_id, Taxon, Replicate) %>%
  dplyr::mutate(rep_num = suppressWarnings(as.integer(stringr::str_remove(Replicate, "^R")))) %>%
  dplyr::arrange(Taxon, rep_num, Replicate) %>%
  dplyr::mutate(curve_code = paste0("C", stringr::str_pad(dplyr::row_number(), width = 3, pad = "0"))) %>%
  dplyr::select(curve_code, series_id, Taxon, Replicate)
raw <- raw %>% dplyr::left_join(curve_key, by = c("series_id", "Taxon", "Replicate"))
readr::write_csv(curve_key, CURVE_KEY_CSV)

# ===== Validate manual overrides =============================================
if (nrow(manual_trim_overrides) > 0) {
  bad <- setdiff(manual_trim_overrides$curve_code, curve_key$curve_code)
  if (length(bad) > 0) stop("Unknown curve_code(s) in manual_trim_overrides: ", paste(bad, collapse = ", "))
}

# ===== Trim each series ======================================================
series_ids  <- unique(raw$series_id)
trimmed_lst <- setNames(vector("list", length(series_ids)), series_ids)
meta_lst    <- setNames(vector("list", length(series_ids)), series_ids)
diag_lst    <- setNames(vector("list", length(series_ids)), series_ids)
skipped_log <- tibble::tibble(curve_code = character(), series_id = character(), reason = character())

for (sid in series_ids) {
  df <- raw %>% dplyr::filter(series_id == sid) %>% dplyr::arrange(Time)
  this_code <- df$curve_code[1]
  out <- trim_one_series(df)
  if (!isTRUE(out$ok)) {
    skipped_log <- tibble::add_row(skipped_log, curve_code = this_code, series_id = sid, reason = out$reason)
    next
  }
  ov <- manual_trim_overrides %>% dplyr::filter(curve_code == this_code)
  if (nrow(ov) > 0) {
    auto_start <- df$Time[out$idx_peak]; auto_end <- df$Time[out$idx_end]
    final_start <- if (!is.na(ov$keep_start[1])) ov$keep_start[1] else auto_start
    final_end   <- if (!is.na(ov$keep_end[1]))   ov$keep_end[1]   else auto_end
    final_end   <- min(final_end, max(df$Time, na.rm = TRUE), MAX_TIME_MIN, na.rm = TRUE)
    if (!is.finite(final_start) || !is.finite(final_end) || final_end <= final_start) {
      skipped_log <- tibble::add_row(skipped_log, curve_code = this_code, series_id = sid,
                                     reason = "Manual override gave invalid final start/end"); next
    }
    manual_trimmed <- df %>% dplyr::filter(Time >= final_start, Time <= final_end)
    if (nrow(manual_trimmed) < 2) {
      skipped_log <- tibble::add_row(skipped_log, curve_code = this_code, series_id = sid,
                                     reason = "Manual override produced < 2 rows"); next
    }
    manual_trimmed <- manual_trimmed %>%
      dplyr::mutate(O2_fit = safe_spline_fit(Time, Oxygen, spar = TRIM_SPAR),
                    curve_code = this_code, series_id = df$series_id[1],
                    peak_time = final_start, end_time_chosen = final_end,
                    end_reason = "manual_override_partial_or_full")
    manual_meta <- tibble::tibble(
      curve_code = this_code, series_id = df$series_id[1], Taxon = df$Taxon[1], Replicate = df$Replicate[1],
      n_points_total = nrow(df), n_points_trimmed = nrow(manual_trimmed),
      peak_time = final_start, main_run_start_time = df$Time[out$idx_run_start],
      main_run_end_time = df$Time[min(out$idx_run_end + 1L, nrow(df))],
      steepest_drop_time = if (is.finite(out$idx_steepest)) df$Time[min(out$idx_steepest + 1L, nrow(df))] else NA_real_,
      chosen_end_time = final_end, fit_duration_min = final_end - final_start,
      delta_Ninoc_to_N0_min = final_start, end_reason = "manual_override_partial_or_full",
      main_run_duration = if (!is.null(out$all_runs)) out$all_runs$duration[1] else NA_real_,
      main_run_drop = if (!is.null(out$all_runs)) out$all_runs$total_drop[1] else NA_real_,
      main_run_score = if (!is.null(out$all_runs)) out$all_runs$score[1] else NA_real_)
    out$data <- manual_trimmed; out$meta <- manual_meta; out$manual_override <- TRUE
    out$manual_keep_start <- if (!is.na(ov$keep_start[1])) ov$keep_start[1] else NA_real_
    out$manual_keep_end   <- if (!is.na(ov$keep_end[1]))   ov$keep_end[1]   else NA_real_
    out$final_start <- final_start; out$final_end <- final_end
  } else {
    out$final_start <- df$Time[out$idx_peak]; out$final_end <- min(df$Time[out$idx_end], MAX_TIME_MIN)
  }
  trimmed_lst[[sid]] <- out$data; meta_lst[[sid]] <- out$meta; diag_lst[[sid]] <- out
}

# ===== Save trimmed data + metadata ==========================================
trimmed <- dplyr::bind_rows(trimmed_lst)
trim_meta <- dplyr::bind_rows(meta_lst)
if (nrow(trimmed) == 0) stop("No curves were successfully trimmed. See Skipped_Series_Log.csv.")
trimmed <- trimmed %>% dplyr::filter(Time <= MAX_TIME_MIN)
readr::write_csv(trimmed, TRIMMED_CSV)
filtered <- trimmed %>% dplyr::select(curve_code, Taxon, Replicate, Time, Oxygen, O2_fit)
readr::write_csv(filtered, FILTERED_CSV)
readr::write_csv(skipped_log, SKIPPED_CSV)
trim_meta <- trim_meta %>%
  dplyr::mutate(rep_num = suppressWarnings(as.integer(stringr::str_remove(Replicate, "^R")))) %>%
  dplyr::arrange(Taxon, rep_num, Replicate) %>% dplyr::select(-rep_num)
readr::write_csv(trim_meta, TRIM_META_CSV)

# ===== Diagnostics PDF =======================================================
pdf_file <- fig("oxygen_trimming_diagnostics.pdf")
if (file.exists(pdf_file)) file.remove(pdf_file)
grDevices::pdf(pdf_file, width = 8.2, height = 6.2)
for (sid in names(diag_lst)) {
  out <- diag_lst[[sid]]; if (is.null(out)) next
  df <- out$full_df %>% dplyr::filter(Time <= MAX_TIME_MIN)
  this_code <- df$curve_code[1]
  xmin_show <- max(min(out$final_start, MAX_TIME_MIN), min(df$Time, na.rm = TRUE))
  xmax_show <- min(out$final_end, MAX_TIME_MIN, max(df$Time, na.rm = TRUE))
  kept_span <- max(xmax_show - xmin_show, 0)
  plot_pad  <- max(PLOT_PAD_MIN_MIN, PLOT_PAD_FRAC * kept_span)
  x_plot_lo <- max(min(df$Time, na.rm = TRUE), xmin_show - plot_pad)
  x_plot_hi <- min(MAX_TIME_MIN, max(df$Time, na.rm = TRUE), xmax_show + plot_pad)
  rect_df <- tibble::tibble(xmin = xmin_show, xmax = xmax_show, ymin = -Inf, ymax = Inf)
  subtitle_text <- if (isTRUE(out$manual_override)) {
    paste0("MANUAL OVERRIDE | kept ", round(xmin_show, 1), " to ", round(xmax_show, 1), " min | spar = ", TRIM_SPAR)
  } else {
    paste0("Black = start | Magenta = chosen end | End reason = ", out$meta$end_reason[1], " | spar = ", TRIM_SPAR)
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(Time, Oxygen)) +
    ggplot2::geom_rect(data = rect_df, inherit.aes = FALSE,
                       ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = "orange", alpha = 0.12) +
    ggplot2::geom_line(colour = "grey60", linewidth = 0.8) +
    ggplot2::geom_line(ggplot2::aes(y = O2_fit), colour = "blue", linewidth = 1.1) +
    ggplot2::coord_cartesian(xlim = c(x_plot_lo, x_plot_hi)) +
    ggplot2::geom_vline(xintercept = xmin_show, colour = "black", linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = xmax_show, colour = "magenta", linewidth = 0.9) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = X_MAJOR_BREAKS_N)) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = Y_MAJOR_BREAKS_N)) +
    ggplot2::labs(title = paste0(this_code, " | ", sid), subtitle = subtitle_text,
                  x = "Time (min)", y = expression(O[2] ~ "(mg L"^{-1} * ")")) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   panel.grid.major = ggplot2::element_line(colour = "grey88", linewidth = 0.25))
  print(p)
}
grDevices::dev.off()

message("02_trimming: kept ", nrow(trim_meta), " curves; skipped ", nrow(skipped_log),
        ". Metadata -> ", TRIM_META_CSV, " ; diagnostics -> ", pdf_file)
