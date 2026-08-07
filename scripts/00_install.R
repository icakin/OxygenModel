# =============================================================================
# 00_install.R - one command to pin and verify the whole environment
# =============================================================================
# Idempotent. Run from the project root:
#
#     Rscript scripts/00_install.R
#
# What it does, in order:
#   1. Makes sure renv is available and restores the project library from
#      renv.lock (the exact package versions the analysis was run with).
#   2. Checks that every package the pipeline actually loads can be loaded.
#   3. Checks the C++ toolchain that stage 12 (rstan) needs in order to compile
#      its Stan model. This is a WARNING, not a failure: stages 01-11 and 13 do
#      not need a compiler, and run_all.sh can skip 12 with OXYMODEL_SKIP_STAN=1.
#   4. Writes env/versions.json recording R, platform, renv and the git commit.
#
# It touches nothing in data/, results/ or scripts/original_scripts/.
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ---- locate project root ----------------------------------------------------
.this_dir <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  d <- if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE)) else
         tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (is.null(d) || is.na(d) || !nzchar(d)) d <- getwd()
  d
})
base_dir <- dirname(normalizePath(.this_dir, mustWork = FALSE))
env_dir  <- file.path(base_dir, "env")
dir.create(env_dir, showWarnings = FALSE, recursive = TRUE)

hr   <- function() message(strrep("=", 78))
step <- function(...) { hr(); message("00_install: ", ...); hr() }

die <- function(title, ...) {
  msg <- paste0(
    "\n", strrep("!", 78), "\n",
    "FATAL: ", title, "\n", strrep("!", 78), "\n",
    paste(c(...), collapse = "\n"), "\n", strrep("!", 78), "\n"
  )
  stop(msg, call. = FALSE)
}

# ---- 1. renv ----------------------------------------------------------------
step("1/4  restoring the R library from renv.lock")

if (!requireNamespace("renv", quietly = TRUE)) {
  message("  renv not found - installing it.")
  utils::install.packages("renv")
}
if (!requireNamespace("renv", quietly = TRUE)) {
  die("renv could not be installed.",
      "Install it by hand and re-run:",
      "    R -e 'install.packages(\"renv\")'")
}

lockfile <- file.path(base_dir, "renv.lock")
if (!file.exists(lockfile)) {
  die("renv.lock is missing.",
      paste0("Expected at: ", lockfile),
      "This file pins every package version. Restore it from git:",
      "    git checkout renv.lock")
}

# Activate the project library, then restore. `prompt = FALSE` keeps it
# non-interactive; restore is a no-op when the library already matches.
renv::activate(project = base_dir)
renv::restore(project = base_dir, prompt = FALSE)

# ---- 2. can every package the pipeline loads actually load? -----------------
step("2/4  checking the packages the pipeline loads")

# One entry per package that a script in scripts/01..13 or config.R attaches.
PIPELINE_PKGS <- c(
  # core / tidyverse
  "tidyverse", "readr", "tibble", "tidyr", "purrr", "dplyr", "ggplot2",
  "stringr", "scales", "glue",
  # modelling
  "zoo", "minpack.lm", "lme4", "lmerTest", "multcomp",
  # figures
  "patchwork", "gridExtra", "ggsignif",
  # apps (sourced headless, but must be installed)
  "shiny",
  # stage 12: joint r-K estimator
  "rstan",
  # environment bookkeeping / IDE detection
  "jsonlite", "rstudioapi"
)

bad <- character(0)
for (p in PIPELINE_PKGS) {
  ok <- suppressWarnings(suppressMessages(
    requireNamespace(p, quietly = TRUE)
  ))
  message(sprintf("  %-12s %s", p, if (ok) "ok" else "MISSING"))
  if (!ok) bad <- c(bad, p)
}
if (length(bad)) {
  die("These packages are in the lockfile but will not load:",
      paste0("    ", paste(bad, collapse = ", ")),
      "Try:",
      "    R -e 'renv::restore(prompt = FALSE)'",
      "If a package fails to compile, install its system dependencies first",
      "(see SETUP.md).")
}

# ---- 3. C++ toolchain for stage 12 (rstan) ----------------------------------
# NOTE ON RcppParallel: the lockfile pins it to 5.1.10 ON PURPOSE. RcppParallel
# 6.x ships a TBB release that dropped `tbb::task_scheduler_init`, which the
# StanHeaders 2.32.x code compiled into every Stan model still calls. With 6.x
# installed, stage 12 dies at stan_model() with
#     symbol not found in flat namespace '__ZN3tbb19task_scheduler_init...'
# Do not bump RcppParallel without checking that stage 12 still compiles.
step("3/4  checking the C++ toolchain that stage 12 needs")

have_tools <- isTRUE(tryCatch(
  pkgbuild::has_build_tools(debug = FALSE), error = function(e) FALSE))

if (have_tools) {
  message("  C++ toolchain  ok  (stage 12 can compile its Stan model)")
} else {
  message("  C++ toolchain  NOT FOUND")
  message("  Stage 12 (12_joint_rK_estimator.R) compiles a Stan model and will")
  message("  fail without one. Everything else runs fine. To install:")
  message("    macOS  : xcode-select --install")
  message("    Linux  : sudo apt-get install -y build-essential")
  message("    Windows: install Rtools matching your R version")
  message("  Or skip that one stage:")
  message("    OXYMODEL_SKIP_STAN=1 bash scripts/run_all.sh")
}

# ---- 4. record the environment ---------------------------------------------
step("4/4  writing env/versions.json")

git_field <- function(args) {
  out <- tryCatch(
    suppressWarnings(system2("git", c("-C", shQuote(base_dir), args),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) NA_character_)
  if (length(out)) trimws(out[1]) else NA_character_
}

lock <- jsonlite::fromJSON(lockfile)
pkg_versions <- vapply(lock$Packages, function(p) as.character(p$Version), character(1))

info <- list(
  recorded_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  git = list(
    commit       = git_field("rev-parse HEAD"),
    short_commit = git_field("rev-parse --short HEAD"),
    branch       = git_field("rev-parse --abbrev-ref HEAD")
  ),
  R = list(
    version   = R.version.string,
    major     = R.version$major,
    minor     = R.version$minor,
    platform  = R.version$platform,
    arch      = R.version$arch,
    os        = utils::sessionInfo()$running
  ),
  renv = list(
    version    = as.character(utils::packageVersion("renv")),
    lockfile   = "renv.lock",
    lockfile_R = lock$R$Version,
    n_packages = length(pkg_versions)
  ),
  cxx_toolchain = list(
    available = have_tools,
    needed_by = "scripts/12_joint_rK_estimator.R (rstan)"
  ),
  key_package_versions = as.list(pkg_versions[intersect(PIPELINE_PKGS, names(pkg_versions))]),
  notes = paste(
    "Written by scripts/00_install.R - do not hand-edit.",
    "Stage 12 needs rstan and a C++ toolchain; run_all.sh can skip it with",
    "OXYMODEL_SKIP_STAN=1."
  )
)

versions_json <- file.path(env_dir, "versions.json")
writeLines(jsonlite::toJSON(info, auto_unbox = TRUE, pretty = TRUE), versions_json)
message("  wrote ", versions_json)

hr()
message("00_install: environment is pinned and verified.")
message("Next:  bash scripts/run_all.sh")
hr()
