#!/usr/bin/env bash
# =============================================================================
# run_originals.sh - run scripts/original_scripts/ in an isolated scratch tree
# =============================================================================
# Comparison D2 needs the numbers the ORIGINAL (pre-restructure) scripts
# produce, on the same deposited inputs. Those scripts write to `data/`,
# `plots/` and `Tables/` RELATIVE TO THE WORKING DIRECTORY, so they are run
# from a scratch directory:
#
#     runs/D1_originals/
#       data/     <- copies of the deposited inputs (read-only to the originals)
#       Tables/   <- their tables land here
#       plots/    <- their figures land here
#       src/      <- COPIES of the original scripts (the originals are never touched)
#
# scripts/original_scripts/ is the provenance record of the published analysis
# and is NEVER modified or executed in place.
#
# One data-provisioning fix is applied, and only in the scratch tree:
# OxygenModel.R hard-codes `data/Ninoc_and_deltaTime_to_N0.csv`, while the
# deposited file is named `data/Ninoc.csv`. The scratch copy is made under BOTH
# names. The file contents are identical; no script is edited.
#
# Usage:  bash reports/D1_baseline/scripts/run_originals.sh
# =============================================================================

set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
SCRATCH="$BASE_DIR/runs/D1_originals"
LOG_DIR="$SCRATCH/logs"

rm -rf "$SCRATCH"
mkdir -p "$SCRATCH/data" "$SCRATCH/Tables" "$SCRATCH/plots" "$SCRATCH/src" "$LOG_DIR"

# ---- inputs -----------------------------------------------------------------
cp "$BASE_DIR"/data/*.csv "$SCRATCH/data/"
# legacy filename the original pipeline expects (identical contents)
cp "$BASE_DIR/data/Ninoc.csv" "$SCRATCH/data/Ninoc_and_deltaTime_to_N0.csv"

# ---- scripts (copies; originals untouched) ----------------------------------
cp "$BASE_DIR"/scripts/original_scripts/*.R "$SCRATCH/src/"

now_ms() {
  if command -v perl >/dev/null 2>&1; then
    perl -MTime::HiRes=time -e 'printf "%.0f\n", time()*1000'
  else
    echo $(( $(date +%s) * 1000 )); fi
}

echo "script,seconds,status" > "$LOG_DIR/timings.csv"

# OxygenModel.R must run first (MC.R consumes its oxygen_results_with_R.csv).
for s in OxygenModel.R MC.R 05CUTOFF.R Simulation.R CUE.R; do
  echo
  echo "-----------------------------------------------------------------------------"
  echo ">>> original: $s"
  echo "-----------------------------------------------------------------------------"
  t0="$(now_ms)"
  set +e
  ( cd "$SCRATCH" && Rscript "src/$s" ) > "$LOG_DIR/${s%.R}.log" 2>&1
  status=$?
  set -e
  t1="$(now_ms)"
  secs="$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", (b-a)/1000}')"
  if [ $status -eq 0 ]; then
    echo "    ok   (${secs}s)"
  else
    echo "    FAILED (exit $status, ${secs}s) - last 25 lines:"
    tail -25 "$LOG_DIR/${s%.R}.log" | sed 's/^/      /'
  fi
  echo "$s,$secs,$([ $status -eq 0 ] && echo ok || echo FAILED)" >> "$LOG_DIR/timings.csv"
done

echo
echo "Original-script run finished. Outputs:"
echo "  tables : $SCRATCH/Tables"
echo "  plots  : $SCRATCH/plots"
echo "  logs   : $LOG_DIR"
ls -1 "$SCRATCH/Tables" 2>/dev/null | sed 's/^/    /' || echo "    (no tables produced)"
