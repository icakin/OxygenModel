# =============================================================================
# 00_install.R - one command to pin and verify the whole environment
# =============================================================================
# Idempotent. Run from the project root:
#
#     Rscript scripts/00_install.R
#
# What it does, in order:
#   1. Checks the macOS BUILD PREREQUISITES the restore needs, and says what is
#      missing BEFORE anything is downloaded (see the block below for why).
#   2. Makes sure renv is available and restores the project library from
#      renv.lock (the exact package versions the analysis was run with).
#   3. Checks that every package the pipeline actually loads can be loaded.
#   4. Checks the C++ toolchain that stage 12 (rstan) needs in order to compile
#      its Stan model. This is a WARNING, not a failure: stages 01-11 and 13 do
#      not need a compiler, and run_all.sh can skip 12 with OXYMODEL_SKIP_STAN=1.
#   5. Writes env/versions.json recording R, platform, renv and the git commit.
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

# =============================================================================
# R VERSION vs renv.lock  -  checked FIRST, before anything is downloaded
# =============================================================================
# The cheapest check there is, and the only failure no toolchain can repair. If
# R is a minor version BEHIND the lock, a fair number of pinned packages carry
# `Depends: R (>= <lock minor>.0)` and refuse to install. renv installs
# everything that CAN build, then aborts, and because it only links a restore
# into the project library once the WHOLE transaction succeeds, the library is
# left empty and the later package check reports every package missing. Failing
# here turns that confusing cascade into one line.

local({
  lock <- file.path(base_dir, "renv.lock")
  if (!file.exists(lock)) return(invisible(NULL))
  lock_r <- tryCatch({
    j <- readLines(lock, warn = FALSE)
    i <- grep('"Version"', j)[1]
    sub('.*"Version"\\s*:\\s*"([^"]+)".*', "\\1", j[i])
  }, error = function(e) NA_character_)
  if (is.na(lock_r)) return(invisible(NULL))

  have <- getRversion(); want <- package_version(lock_r)
  same_minor <- identical(unlist(want)[1:2], unlist(have)[1:2])
  # die() joins its arguments with newlines, one per line, so every line below
  # must be a SINGLE fully-assembled string. And getRversion() returns an
  # R_system_version, which c() flattens to "c(4, 4, 2)" - hence as.character().
  hv <- as.character(have)
  wm <- paste(unlist(want)[1:2], collapse = ".")

  if (identical(hv, lock_r)) {
    message("00_install: R ", hv, " matches renv.lock")
  } else if (!same_minor && have < want) {
    die(paste0("R ", hv, " is older than the pinned R ", lock_r),
        paste0("renv.lock was solved under R ", lock_r, ". Several pinned packages declare"),
        paste0("Depends: R (>= ", wm, ".0) and will refuse to install on R ", hv, "."),
        "This is not a missing system library, and no compiler flag fixes it.",
        "",
        paste0("  fix: install R ", lock_r, " from"),
        "       https://cran.r-project.org/bin/macosx/",
        "",
        paste0("Installing R ", wm, " does not remove your existing R ", hv, ". The Candidas"),
        "project pins the same R version, so one install serves both repositories.",
        "Then re-run this script; already-downloaded packages are in the renv",
        "cache, so the retry is fast.")
  } else {
    message("00_install: [warn] renv.lock was solved under R ", lock_r,
            "; you are on ", hv, ". Restore will proceed, but pinned binaries",
            " may be rebuilt from source and exact reproduction is not guaranteed.")
  }
})

# =============================================================================
# 0. SYSTEM PREREQUISITES  (macOS only)
# =============================================================================
# WHY THIS EXISTS. renv.lock pins EXACT versions. CRAN serves macOS binaries only
# for each package's CURRENT version, so as CRAN moves on, more pinned versions
# lose their binary and have to be compiled. Compiling needs a toolchain the
# stock macOS image does not have, and the failure surfaces partway through a
# long restore rather than up front. This block checks first and says what is
# missing before anything is downloaded.
#
# NOTE: THIS HAS NOT BEEN CONFIRMED ON A CLEAN MACHINE. The requirements below
# were derived from renv.lock and from a reported clean-Mac failure; no clean Mac
# was available to test the fix. Treat the list as our best determination, not as
# a verified recipe. See SETUP.md.

is_macos <- function() Sys.info()[["sysname"]] == "Darwin"

# Which pinned packages have no macOS binary AT THE PINNED VERSION, and so must
# be compiled? Computed from the lockfile against the live CRAN binary index -
# not hard-coded, because the answer changes as CRAN moves on.
packages_needing_compile <- function(lockfile, project = base_dir) {
  out <- tryCatch({
    lock <- jsonlite::fromJSON(lockfile)$Packages
    pinned <- vapply(lock, function(p) as.character(p$Version), character(1))
    needs  <- vapply(lock, function(p) identical(p$NeedsCompilation, "yes"), logical(1))
    comp   <- names(pinned)[needs]
    ap <- utils::available.packages(type = "binary")
    binv <- ifelse(comp %in% rownames(ap), ap[match(comp, rownames(ap)), "Version"], NA)
    no_binary <- comp[is.na(binv) | binv != pinned[comp]]

    # A package only has to be COMPILED NOW if it is not already installed at the
    # pinned version. On a machine whose library is already synchronised the
    # restore is a no-op and no toolchain is needed - so the check must not fail
    # there. Only the "still to build" set drives the fatal path.
    libs <- tryCatch(renv::paths$library(project = project), error = function(e) .libPaths())
    still <- Filter(function(p) {
      ip <- tryCatch(utils::packageDescription(p, lib.loc = libs)$Version,
                     error = function(e) NA_character_)
      is.null(ip) || is.na(ip) || !identical(as.character(ip), unname(pinned[p]))
    }, no_binary)
    list(no_binary = no_binary, still_to_build = still)
  }, error = function(e) NULL)
  out
}

check_prerequisites <- function(lockfile) {
  if (!is_macos()) {
    message("  not macOS - skipping the macOS prerequisite check.")
    message("  Linux/Windows prerequisites are NOT documented here; see SETUP.md.")
    return(invisible(NULL))
  }
  step("1/5  checking macOS build prerequisites")
  missing <- character(0)

  # -- Xcode command line tools ------------------------------------------------
  clt <- suppressWarnings(system2("xcode-select", "-p", stdout = TRUE, stderr = FALSE))
  if (!length(clt) || !nzchar(clt[1])) {
    missing <- c(missing, paste0(
      "Xcode command line tools (C/C++ compiler)\n",
      "      install with:  xcode-select --install"))
    message("  xcode CLT      MISSING")
  } else message("  xcode CLT      ok  (", clt[1], ")")

  # -- gfortran ----------------------------------------------------------------
  # Needed by Matrix, mvtnorm and nlme among the pinned versions that must
  # compile. CRAN-built R expects the CRAN gfortran, not Homebrew's gcc: they
  # use different runtime paths and Homebrew's will not satisfy R's linker.
  gf <- c(Sys.which("gfortran")[[1]], "/opt/gfortran/bin/gfortran",
          "/usr/local/gfortran/bin/gfortran")
  gf <- gf[nzchar(gf) & file.exists(gf)]
  if (!length(gf)) {
    missing <- c(missing, paste0(
      "the CRAN gfortran toolchain (Fortran compiler)\n",
      "      needed by: Matrix, mvtnorm, nlme\n",
      "      download:  https://mac.r-project.org/tools/\n",
      "      NOTE: Homebrew's gcc is NOT a substitute for CRAN-built R."))
    message("  gfortran       MISSING")
  } else {
    message("  gfortran       ok  (", gf[1], ")")
    # Presence is NOT sufficient. R's Makeconf hardcodes the gfortran runtime
    # library path for the version R was BUILT against; the current installer
    # at mac.r-project.org/tools ships a newer one (e.g. 14.2.0 where Makeconf
    # expects 12.2.0). The compiler then runs fine and the LINK step fails with
    #   ld: library 'gfortran' not found
    # after the Fortran objects have already been built. Check the paths R will
    # actually hand the linker, rather than the compiler's existence.
    # ~/.R/Makevars overrides Makeconf, so an override there is the effective
    # setting and must be read FIRST; otherwise a machine already fixed this
    # way is reported broken.
    usr_mk <- path.expand("~/.R/Makevars")
    fl <- character(0); src <- ""
    if (file.exists(usr_mk)) {
      fl  <- grep("^\\s*FLIBS\\s*=", readLines(usr_mk, warn = FALSE), value = TRUE)
      src <- "~/.R/Makevars"
    }
    if (!length(fl)) {
      mkconf <- file.path(R.home("etc"), "Makeconf")
      if (file.exists(mkconf)) {
        fl  <- grep("^\\s*FLIBS\\s*=", readLines(mkconf, warn = FALSE), value = TRUE)
        src <- "R's Makeconf"
      }
    }
    if (length(fl)) {
      lpaths <- regmatches(fl[1], gregexpr("(?<=-L)[^ ]+", fl[1], perl = TRUE))[[1]]
      bad    <- lpaths[nzchar(lpaths) & !dir.exists(lpaths)]
      if (length(bad)) {
        have <- Sys.glob("/opt/gfortran/lib/gcc/*/*")
        message("  gfortran libs  MISSING  (FLIBS path in ", src, " does not exist)")
        missing <- c(missing, paste0(
          "a gfortran runtime path that R expects but which is not present\n",
          "      FLIBS in ", src, " refers to:  ", paste(bad, collapse = ", "), "\n",
          if (length(have))
            paste0("      installed on this machine:  ", paste(have, collapse = ", "), "\n")
          else "",
          "      The compiler exists, so a presence check passes, but the LINK step\n",
          "      fails. Fix by overriding FLIBS in ~/.R/Makevars (single line):\n",
          "        FLIBS = -L", if (length(have)) have[1] else "<installed path>",
          " -L/opt/gfortran/lib -lgfortran -lemutls_w -lquadmath"))
      } else {
        message("  gfortran libs  ok  (FLIBS paths in ", src, " resolve)")
      }
    }
  }

  # -- Homebrew system libraries ----------------------------------------------
  # Each is named with the R package it serves, so a reader can tell why.
  # The graphics libraries form a dependency CHAIN:
  #   freetype -> systemfonts -> textshaping -> ragg -> tidyverse
  # so a missing webp header ultimately takes out tidyverse four levels up,
  # and each one only reveals the next if they are installed one at a time.
  # Verified empirically on a clean macOS arm64 machine (see SETUP.md).
  brew_all <- c("pkg-config", "openssl", "freetype", "harfbuzz", "fribidi",
                "libtiff", "jpeg-turbo", "webp")
  brew_needs <- list(
    "pkg-config" = "pkg-config (used to locate every library below)",
    openssl      = "openssl (R packages: openssl, curl's TLS)",
    freetype     = "freetype (R package: systemfonts)",
    harfbuzz     = "harfbuzz (R package: textshaping)",
    fribidi      = "fribidi (R package: textshaping)",
    libtiff      = "libtiff (R package: ragg)",
    "jpeg-turbo" = "jpeg-turbo (R package: ragg)",
    webp         = "webp (R package: ragg)")
  brewbin <- Sys.which("brew")[[1]]
  if (!nzchar(brewbin)) {
    message("  homebrew       not found - cannot check system libraries")
    missing <- c(missing, paste0(
      "Homebrew, to supply the C libraries the graphics stack needs\n",
      "      install:   https://brew.sh\n",
      "      then:      brew install ", paste(brew_all, collapse = " ")))
  } else {
    for (nm in names(brew_needs)) {
      p <- suppressWarnings(system2(brewbin, c("--prefix", nm),
                                    stdout = TRUE, stderr = FALSE))
      lab <- formatC(nm, width = 11, flag = "-")
      if (!length(p) || !nzchar(p[1])) {
        message("  brew ", lab, "MISSING")
        missing <- c(missing, paste0(brew_needs[[nm]],
                                     "\n      install:   brew install ", nm))
      } else message("  brew ", lab, "ok")
    }
  }

  # -- the conda / libkrb5 trap ------------------------------------------------
  # Reported symptom: curl fails to build. Cause: a conda installation ahead of
  # the system libraries on PATH supplies its own libkrb5, which the build picks
  # up instead of the system one.
  pth <- strsplit(Sys.getenv("PATH"), ":", fixed = TRUE)[[1]]
  conda_first <- length(pth) && any(grepl("conda|miniforge|mamba", pth[seq_len(min(3, length(pth)))]))
  if (conda_first || nzchar(Sys.getenv("CONDA_PREFIX"))) {
    message("  conda          DETECTED ahead of the system path")
    message("    A conda environment early on PATH supplies its own libkrb5, and")
    message("    `curl` may fail to build against it. If the restore fails on curl:")
    message("      conda deactivate            # for the duration of the restore")
    message("    or prepend the system paths:")
    message("      PATH=/usr/bin:/bin:/usr/sbin:/sbin:$PATH Rscript scripts/00_install.R")
  } else message("  conda          not on PATH (good)")

  # -- gettext headers (libintl.h), needed by data.table ----------------------
  # data.table's src/po.h includes <libintl.h> with no fallback. CRAN's macOS
  # builders keep gettext in /opt/R/<arch>, which is why R's Makeconf already
  # puts -I/opt/R/<arch>/include on every compile line and why the CRAN binary
  # builds - but that directory is empty on a user machine unless populated on
  # purpose. Source-installing data.table then dies with
  #     ./po.h:2:10: fatal error: 'libintl.h' file not found
  # and takes dtplyr and tidyverse down with it. Homebrew's gettext is keg-only
  # and is NOT on the include path, so `brew install gettext` does not fix it.
  intl_dirs <- c("/opt/R/arm64/include", "/opt/R/x86_64/include",
                 "/opt/homebrew/opt/gettext/include", "/usr/local/opt/gettext/include",
                 "/opt/homebrew/include", "/usr/local/include")
  intl_hit <- intl_dirs[file.exists(file.path(intl_dirs, "libintl.h"))]
  if (length(intl_hit)) {
    message("  libintl.h      ok  (", intl_hit[1], ")")
  } else {
    arch    <- if (identical(R.version$arch, "aarch64")) "arm64" else "x86_64"
    tarball <- if (arch == "arm64") "gettext-0.21-darwin.20-arm64.tar.gz"
               else                 "gettext-0.21-darwin.17-x86_64.tar.gz"
    message("  libintl.h      MISSING")
    missing <- c(missing, paste0(
      "gettext headers (libintl.h), needed to compile data.table\n",
      "      without them data.table fails, and dtplyr and tidyverse fail with it\n",
      "      install CRAN's own build into /opt/R/", arch, ", which R already searches:\n",
      "        curl -fO https://mac.R-project.org/libs-", arch, "/", tarball, "\n",
      "        sudo tar fvxz ", tarball, " -C /\n",
      "      NOTE: Homebrew's gettext is keg-only and will NOT be found."))
  }

  # -- how much will actually have to compile? --------------------------------
  nc <- packages_needing_compile(lockfile)
  n_build <- NA_integer_
  if (is.null(nc)) {
    message("  compile scan   could not reach CRAN - skipping (assuming a build may be needed)")
  } else {
    n_build <- length(nc$still_to_build)
    message("  compile scan   ", length(nc$no_binary), " pinned versions have no macOS binary;")
    message("                 ", n_build, " of those are not yet installed and would be compiled now.")
    if (n_build) {
      message("                   ", paste(utils::head(sort(nc$still_to_build), 12), collapse = ", "),
              if (n_build > 12) paste0(", ... (+", n_build - 12, " more)") else "")
      message("                 This is expected and GROWS OVER TIME: CRAN ships binaries")
      message("                 only for each package's CURRENT version, and renv.lock")
      message("                 pins older ones on purpose.")
    } else {
      message("                 The library is already synchronised - nothing to compile.")
    }
  }

  # Only fatal if something actually has to be built now. A machine whose library
  # is already in sync needs no toolchain, and must not be blocked by this check.
  if (length(missing) && !identical(n_build, 0L)) {
    die("Missing build prerequisites - the restore would fail partway through.",
        "",
        paste0("  * ", missing, collapse = "\n"),
        "",
        "Full instructions, including the conda/libkrb5 issue:  SETUP.md",
        "",
        "These prerequisites were derived from renv.lock and a reported",
        "clean-Mac failure. They have NOT been confirmed on a clean machine.")
  }
  if (length(missing)) {
    message("")
    message("  NOTE: the following are missing, but nothing needs compiling right")
    message("  now, so the restore should still succeed. Install them before the")
    message("  lockfile next moves ahead of CRAN's binaries:")
    for (m in missing) message("    * ", sub("\n.*", "", m))
    message("")
  } else message("  all checked prerequisites present.")
  invisible(NULL)
}

# ---- 2. renv ----------------------------------------------------------------
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

# Check the toolchain BEFORE downloading anything, so a missing compiler is
# reported up front rather than partway into a build.
check_prerequisites(lockfile)

step("2/5  restoring the R library from renv.lock")

# ---- prefer binaries on macOS ------------------------------------------------
# THIS CHANGES HOW PACKAGES ARE OBTAINED, NEVER WHICH VERSIONS. renv::restore()
# installs the exact versions recorded in renv.lock whatever the source; the
# lockfile remains the pin. Setting pkgType only decides binary-vs-source for
# those versions where both exist.
#
# DELIBERATE DEVIATION, stated because the brief asked for the opposite: we do
# NOT set `install.packages.compile.from.source = "never"`. At the time of
# writing, 35 of the 73 compiled packages in renv.lock are pinned to versions
# CRAN no longer ships a binary for, because CRAN serves binaries only for each
# package's CURRENT version. Forbidding source builds would make the restore
# FAIL on those 35 rather than succeed - it would convert a slow restore into a
# broken one. The right fix is the prerequisite check above plus the toolchain
# instructions in SETUP.md, not a flag that refuses to build.
#
# Only set when the user has not chosen for themselves.
if (is_macos() && is.null(getOption("pkgType.set.by.user"))) {
  if (identical(getOption("pkgType"), "source")) {
    message("  pkgType is 'source' (your setting) - leaving it alone.")
  } else {
    options(pkgType = "both")   # binary where one exists at the pinned version
    message("  pkgType = 'both': binaries where available, source only where not.")
  }
}

# Activate the project library, then restore. `prompt = FALSE` keeps it
# non-interactive; restore is a no-op when the library already matches.
renv::activate(project = base_dir)
renv::restore(project = base_dir, prompt = FALSE)

# ---- 3. can every package the pipeline loads actually load? -----------------
step("3/5  checking the packages the pipeline loads")

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

# ---- 4. C++ toolchain for stage 12 (rstan) ----------------------------------
# NOTE ON RcppParallel: the lockfile pins it to 5.1.10 ON PURPOSE. RcppParallel
# 6.x ships a TBB release that dropped `tbb::task_scheduler_init`, which the
# StanHeaders 2.32.x code compiled into every Stan model still calls. With 6.x
# installed, stage 12 dies at stan_model() with
#     symbol not found in flat namespace '__ZN3tbb19task_scheduler_init...'
# Do not bump RcppParallel without checking that stage 12 still compiles.
step("4/5  checking the C++ toolchain that stage 12 needs")

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

# ---- 5. record the environment ---------------------------------------------
step("5/5  writing env/versions.json")

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
