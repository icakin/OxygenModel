# =============================================================================
# 01_longdata.R - Load and tidy the raw long-format oxygen table
# =============================================================================
# Reads data/Oxygen_Data_Long.csv, coerces types, sorts, adds a series_id, and
# writes results/tables/Oxygen_All_Long.csv for the rest of the pipeline.
#
#   Input : data/Oxygen_Data_Long.csv  (columns: Taxon, Replicate, Time, Oxygen)
#   Output: results/tables/Oxygen_All_Long.csv
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

# ---- load raw oxygen --------------------------------------------------------
raw <- readr::read_csv(OXYGEN_LONG_CSV, show_col_types = FALSE) %>%
  dplyr::mutate(
    Taxon     = as.character(Taxon),
    Replicate = as.character(Replicate),
    Time      = as.numeric(Time),
    Oxygen    = as.numeric(Oxygen)
  ) %>%
  dplyr::arrange(Taxon, Replicate, Time)

stopifnot(all(c("Taxon", "Replicate", "Time", "Oxygen") %in% names(raw)))

raw <- raw %>%
  dplyr::mutate(series_id = paste(Taxon, Replicate, sep = " | "))

readr::write_csv(raw, LONG_CSV)

message("01_longdata: ", nrow(raw), " rows across ",
        dplyr::n_distinct(raw$series_id), " series -> ", LONG_CSV)
