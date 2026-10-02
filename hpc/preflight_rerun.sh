#!/bin/bash
# Verify that a clean Torch project has one current R script for each 00-43 step.

set -euo pipefail

PROJECT_DIR="${SIM_PROJECT_DIR:-$PWD}"
ANALYTIC_RDS="${SIM_ANALYTIC_RDS:-$PROJECT_DIR/data/MIDUS_discrimination_analysis.rds}"

cd "$PROJECT_DIR"

PROBLEMS=0

echo "MIDUS two-phase clean-rerun preflight"
echo "Project: $PROJECT_DIR"
echo

for NUMBER in $(seq -w 0 43); do
  shopt -s nullglob
  MATCHES=("${NUMBER}"_*.R)
  shopt -u nullglob

  if [[ "${#MATCHES[@]}" -eq 1 ]]; then
    printf 'PASS  %s\n' "${MATCHES[0]}"
  elif [[ "${#MATCHES[@]}" -eq 0 ]]; then
    printf 'FAIL  no %s_*.R script found\n' "$NUMBER"
    PROBLEMS=$((PROBLEMS + 1))
  else
    printf 'FAIL  multiple %s_*.R scripts found: %s\n' "$NUMBER" "${MATCHES[*]}"
    PROBLEMS=$((PROBLEMS + 1))
  fi
done

echo

if [[ -f "$ANALYTIC_RDS" ]]; then
  echo "PASS  analysis input: $ANALYTIC_RDS"
else
  echo "FAIL  analysis input missing: $ANALYTIC_RDS"
  PROBLEMS=$((PROBLEMS + 1))
fi

REQUIRED_HPC_FILES=(
  "hpc/env.sh"
  "hpc/setup_project.sh"
  "hpc/calibrate.sbatch"
  "hpc/submit_calibrate.sh"
  "hpc/smoke.sbatch"
  "hpc/submit_smoke.sh"
  "hpc/cache.sbatch"
  "hpc/submit_cache.sh"
  "hpc/primary_array.sbatch"
  "hpc/submit_primary.sh"
  "hpc/resubmit_missing.sh"
  "hpc/content_qc.sbatch"
  "hpc/combine.sbatch"
  "hpc/summary.sbatch"
  "hpc/submit_postprocess.sh"
)

for FILE in "${REQUIRED_HPC_FILES[@]}"; do
  if [[ -f "$FILE" ]]; then
    echo "PASS  $FILE"
  else
    echo "FAIL  missing $FILE"
    PROBLEMS=$((PROBLEMS + 1))
  fi
done

echo

if [[ "$PROBLEMS" -ne 0 ]]; then
  echo "PREFLIGHT FAILED: $PROBLEMS problem(s) found."
  exit 2
fi

echo "PREFLIGHT PASSED. The project is ready for setup/calibration."
