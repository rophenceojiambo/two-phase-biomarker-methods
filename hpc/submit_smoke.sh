#!/bin/bash
# Submit the final-setting HPC smoke test.

set -euo pipefail

: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT before submitting.}"

PROJECT_DIR="${SIM_PROJECT_DIR:-$PWD}"
CALIBRATION_RDS="${SIM_CALIBRATION_RDS:-$PROJECT_DIR/calibration/midus_calibration.rds}"

if [[ ! -f "$CALIBRATION_RDS" ]]; then
  echo "Calibration object not found: $CALIBRATION_RDS"
  echo "Run hpc/submit_calibrate.sh first."
  exit 2
fi

mkdir -p "$PROJECT_DIR/logs"

JOB_ID=$(
  sbatch \
    --parsable \
    --account="$SLURM_ACCOUNT" \
    --chdir="$PROJECT_DIR" \
    --export="ALL,SIM_PROJECT_DIR=$PROJECT_DIR,SIM_CALIBRATION_RDS=$CALIBRATION_RDS" \
    hpc/smoke.sbatch
)

JOB_ID="${JOB_ID%%;*}"
echo "Submitted smoke-test job: $JOB_ID"
echo "Monitor with: squeue -j $JOB_ID"
