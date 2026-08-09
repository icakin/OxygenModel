#!/usr/bin/env bash
# =============================================================================
# run_all.sh - run the whole OxygenModel pipeline unattended
# =============================================================================
# Usage, from the project root:
#
#     bash scripts/run_all.sh
#
# Zero user interaction. Every stage is timed, every stage's console output is
# tee'd into logs/<run-id>/ and into one combined log. A non-zero exit from any
# stage aborts the run (set -e) - there are no silent skips.
#
# This is the shell companion to scripts/run_all.R. It differs in three ways:
#   * it runs each stage in its OWN R process, so one stage cannot leak objects
#     or options into the next (run_all.R source()s them into one session);
#   * it covers stages 11, 12 and 13, which run_all.R does not list;
#   * it records per-stage wall time and exit status to a CSV.
#
# The two Shiny apps (03 trim-selector, 04 experiment-inputs) are SOURCED but
# never launched: their committed outputs
#     results/tables/manual_fit_windows.csv
#     results/tables/plot_exclude_points.csv
#     data/taxon_cell_sizes.csv
# are manual curation and are treated as INPUTS to this run. Sourcing them
# proves they still parse and that their inputs are present.
#
# Stage 12 needs rstan and compiles a Stan model on first use; it is the
# slowest stage by a wide margin. Set OXYMODEL_SKIP_STAN=1 to skip it.
# =============================================================================

set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$BASE_DIR"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_DIR="$BASE_DIR/logs/$RUN_ID"
mkdir -p "$LOG_DIR"
MAIN_LOG="$LOG_DIR/run_all.log"
TIMING_CSV="$LOG_DIR/stage_timings.csv"

# Millisecond clock; falls back to whole seconds if perl is unavailable.
now_ms() {
  if command -v perl >/dev/null 2>&1; then
    perl -MTime::HiRes=time -e 'printf "%.0f\n", time()*1000'
  else
    echo $(( $(date +%s) * 1000 ))
  fi
}

# Everything from here on is tee'd into the combined log.
exec > >(tee -a "$MAIN_LOG") 2>&1

echo "============================================================================="
echo "OxygenModel pipeline - unattended run"
echo "  run id       : $RUN_ID"
echo "  project      : $BASE_DIR"
echo "  git commit   : $(git -C "$BASE_DIR" rev-parse HEAD 2>/dev/null || echo 'not a git repo')"
echo "  git status   : $(git -C "$BASE_DIR" status --porcelain 2>/dev/null | wc -l | tr -d ' ') modified/untracked paths"
echo "  R            : $(Rscript -e 'cat(R.version.string)' 2>/dev/null)"
echo "  logs         : $LOG_DIR"
echo "  started      : $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
echo "============================================================================="

echo "stage,script,mode,seconds,status" > "$TIMING_CSV"

RUN_START="$(now_ms)"
FAILED_STAGE=""

on_exit() {
  local rc=$?
  local total=$(( ($(now_ms) - RUN_START) / 1000 ))
  echo
  echo "============================================================================="
  if [ $rc -eq 0 ]; then
    echo "PIPELINE COMPLETE - total wall time ${total}s"
  else
    echo "PIPELINE FAILED at stage '${FAILED_STAGE:-unknown}' (exit $rc) after ${total}s"
    echo "See $LOG_DIR for the per-stage logs."
  fi
  echo "  finished : $(date -u '+%Y-%m-%d %H:%M:%S UTC')"
  echo "  timings  : $TIMING_CSV"
  echo "  log      : $MAIN_LOG"
  echo "============================================================================="
  # Give the tee subshell a moment to flush before the shell exits.
  sleep 0.2
}
trap on_exit EXIT

# run_stage <label> <script> <mode>
#   mode = "run"    : execute the script; a failure aborts the pipeline
#   mode = "source" : parse/source-only check (the Shiny apps); must not launch
#   mode = "skip"   : recorded as intentionally skipped; not executed
run_stage() {
  local label="$1" script="$2" mode="$3"
  local stage_log="$LOG_DIR/${label}.log"
  FAILED_STAGE="$label"

  if [ "$mode" = "skip" ]; then
    echo
    echo ">>> [$label] $script  (SKIPPED on request)"
    echo "$label,$script,skip,0.00,skipped" >> "$TIMING_CSV"
    FAILED_STAGE=""
    return 0
  fi

  echo
  echo "-----------------------------------------------------------------------------"
  echo ">>> [$label] $script  (mode=$mode)"
  echo "-----------------------------------------------------------------------------"

  local t0 t1 secs status
  t0="$(now_ms)"

  if [ "$mode" = "source" ]; then
    # Explicitly non-interactive; the apps' headless guard keys off this.
    OXYMODEL_LAUNCH_APPS=0 Rscript -e \
      "cat('sourcing (headless check):', '$script', '\n'); source('$script', echo = FALSE)" \
      2>&1 | tee "$stage_log"
  else
    Rscript "$script" 2>&1 | tee "$stage_log"
  fi
  status="${PIPESTATUS[0]}"

  t1="$(now_ms)"
  secs="$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", (b-a)/1000}')"

  echo "$label,$script,$mode,$secs,$([ "$status" = "0" ] && echo ok || echo FAILED)" >> "$TIMING_CSV"
  echo "<<< [$label] finished in ${secs}s (exit $status)"

  if [ "$status" != "0" ]; then
    echo "ERROR: stage $label failed."
    exit "$status"
  fi
}

STAN_MODE="run"
if [ "${OXYMODEL_SKIP_STAN:-0}" = "1" ]; then STAN_MODE="skip"; fi

run_stage 01_longdata                 scripts/01_longdata.R                 run
run_stage 02_trimming                 scripts/02_trimming.R                 run
run_stage 03_trim_selector            scripts/03_trim_selector.R            source
run_stage 04_experiment_inputs        scripts/04_experiment_inputs.R        source
run_stage 05_oxygen_fits              scripts/05_oxygen_fits.R              run
run_stage 06_main_figures             scripts/06_main_figures.R             run
run_stage 07_cutoff_sensitivity       scripts/07_cutoff_sensitivity.R       run
run_stage 08_window_sensitivity       scripts/08_window_sensitivity.R       run
run_stage 09_montecarlo_N0            scripts/09_montecarlo_N0.R            run
run_stage 10_simulation_recovery      scripts/10_simulation_recovery.R      run
run_stage 11_temperature_cue          scripts/11_temperature_cue.R          run
run_stage 12_joint_rK_estimator       scripts/12_joint_rK_estimator.R       "$STAN_MODE"
run_stage 13_depletion_frac_sensitivity scripts/13_depletion_frac_sensitivity.R run

FAILED_STAGE=""

echo
echo "Per-stage wall time:"
column -s, -t "$TIMING_CSV" 2>/dev/null || cat "$TIMING_CSV"

# Surface anything that looks like a problem, so nothing is buried in the log.
echo
echo "Grep of the combined log for Error / Warning / skipped:"
grep -nE "^Error|Error in|Error:|Execution halted|skipped|Skipped|Warning message" "$MAIN_LOG" \
  | grep -v "run_all.sh" || echo "  (none)"
