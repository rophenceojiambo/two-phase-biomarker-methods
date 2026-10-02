#!/bin/bash
# Submit the MIDUS calibration job for a clean rerun.

set -euo pipefail

: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT before submitting.}"

PROJECT_DIR="${SIM_PROJECT_DIR:-$PWD}"
ANALYTIC_RDS="${SIM_ANALYTIC_RDS:-$PROJECT_DIR/data/MIDUS_discrimination_analysis.rds}"

if [[ ! -f "$PROJECT_DIR/01_calibrate_midus.R" ]]; then
  echo "Required script not found: $PROJECT_DIR/01_calibrate_midus.R"
  exit 2
fi

if [[ ! -f "$ANALYTIC_RDS" ]]; then
  echo "Analysis-ready input data not found: $ANALYTIC_RDS"
  exit 2
fi

mkdir -p "$PROJECT_DIR/logs"

JOB_ID=$(
  sbatch \
    --parsable \
    --account="$SLURM_ACCOUNT" \
    --chdir="$PROJECT_DIR" \
    --export="ALL,SIM_PROJECT_DIR=$PROJECT_DIR,SIM_ANALYTIC_RDS=$ANALYTIC_RDS" \
    hpc/calibrate.sbatch
)

JOB_ID="${JOB_ID%%;*}"
echo "Submitted calibration job: $JOB_ID"
echo "Monitor with: squeue -j $JOB_ID"
