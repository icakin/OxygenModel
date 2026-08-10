# =============================================================================
# D5 PART 0 - does the pipeline still reproduce results/ after 12 was added?
# =============================================================================
# Compares every regenerated results/tables/*.csv against the committed copy and
# records, per file, whether it is byte-identical and its maximum absolute and
# relative numeric difference.
#
# Run IMMEDIATELY after `bash scripts/run_all.sh`, while results/ still holds the
# regenerated tree. Pass the committed tree as argument 1.
#
#   Rscript .../d5_part0_reproduction.R <committed_tables_dir>
#
# OUTPUT  runs/D5_analysis/P0_results_reproduction.csv
# =============================================================================

source(file.path(dirname(sub("^--file=", "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "d5_common.R"))

args <- commandArgs(TRUE)
ref  <- if (length(args)) args[1] else stop("give the committed results/tables dir")
new  <- tables_dir

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
csvs <- function(d) list.files(d, pattern = "[.]csv$")
fs <- sort(intersect(csvs(ref), csvs(new)))

out <- do.call(rbind, lapply(fs, function(f) {
  pa <- file.path(ref, f); pb <- file.path(new, f)
  ident <- identical(unname(tools::md5sum(pa)), unname(tools::md5sum(pb)))
  a <- try(read.csv(pa, stringsAsFactors = FALSE, check.names = FALSE), silent = TRUE)
  b <- try(read.csv(pb, stringsAsFactors = FALSE, check.names = FALSE), silent = TRUE)
  ma <- NA_real_; mr <- NA_real_
  if (!inherits(a, "try-error") && !inherits(b, "try-error") &&
      identical(dim(a), dim(b)) && identical(names(a), names(b))) {
    ma <- 0; mr <- 0
    for (cn in names(a)) {
      an <- num(a[[cn]]); bn <- num(b[[cn]])
      ok <- is.finite(an) & is.finite(bn)
      if (any(ok)) {
        d   <- abs(an[ok] - bn[ok])
        den <- pmax(abs(an[ok]), .Machine$double.eps)
        ma  <- max(ma, max(d)); mr <- max(mr, max(d / den))
      }
    }
  }
  data.frame(file = f, byte_identical = ident,
             max_abs_diff = ma, max_rel_diff = mr, stringsAsFactors = FALSE)
}))

out <- out[order(-ifelse(is.na(out$max_rel_diff), Inf, out$max_rel_diff)), ]
out$verdict <- ifelse(out$byte_identical, "byte-identical",
                ifelse(!is.na(out$max_rel_diff) & out$max_rel_diff <= 1e-9,
                       "numerically identical (<=1e-9)", "differs"))

d5_write(tibble::as_tibble(out), "P0_results_reproduction.csv")

cat(sprintf("\n  %d files | %d byte-identical | %d numerically identical | %d differ\n",
            nrow(out), sum(out$byte_identical),
            sum(!out$byte_identical & out$verdict == "numerically identical (<=1e-9)"),
            sum(out$verdict == "differs")))
d <- out[out$verdict == "differs", c("file", "max_abs_diff", "max_rel_diff")]
if (nrow(d)) print(d, row.names = FALSE)
