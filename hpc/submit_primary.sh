#!/bin/bash
# Submit the 324-scenario primary simulation array after preflight checks.

set -euo pipefail

: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT before submitting.}"

PROJECT_DIR="${SIM_PROJECT_DIR:-$PWD}"
MAX_CONCURRENT="${MAX_CONCURRENT:-2000}"
N_TASKS="${N_TASKS:-6480}"
CALIBRATION_RDS="${SIM_CALIBRATION_RDS:-$PROJECT_DIR/calibration/midus_calibration.rds}"
SMOKE_PASS_FILE="$PROJECT_DIR/results/HPC_SMOKE_TEST_PASS.txt"
CACHE_DIR="$PROJECT_DIR/calibration/dgm_cache"

for VALUE_NAME in MAX_CONCURRENT N_TASKS; do
  VALUE="${!VALUE_NAME}"
  if ! [[ "$VALUE" =~ ^[1-9][0-9]*$ ]]; then
    echo "$VALUE_NAME must be a positive integer."
    exit 2
  fi
done

if [[ "$N_TASKS" -ne 6480 ]]; then
  echo "The primary production design requires N_TASKS=6480."
  exit 2
fi

REQUIRED_FILES=(
  "00_config.R"
  "02_dgm.R"
  "03_methods.R"
  "05_run_primary_chunk.R"
  "hpc/primary_array.sbatch"
)

for REQUIRED_FILE in "${REQUIRED_FILES[@]}"; do
  if [[ ! -f "$PROJECT_DIR/$REQUIRED_FILE" ]]; then
    echo "Required file not found: $PROJECT_DIR/$REQUIRED_FILE"
    exit 2
  fi
done

if [[ ! -f "$CALIBRATION_RDS" ]]; then
  echo "Calibration object not found: $CALIBRATION_RDS"
  exit 2
fi

if [[ ! -f "$SMOKE_PASS_FILE" ]] || \
   ! grep -q "HPC SMOKE TEST: PASS" "$SMOKE_PASS_FILE"; then
  echo "Smoke-test PASS marker is missing or invalid:"
  echo "$SMOKE_PASS_FILE"
  echo "Run and review hpc/submit_smoke.sh first."
  exit 2
fi

if [[ ! -d "$CACHE_DIR" ]]; then
  echo "DGM cache directory not found: $CACHE_DIR"
  exit 2
fi

CACHE_COUNT=$(find "$CACHE_DIR" -maxdepth 1 -type f -name 'dgm_f*_rA*_rY*_th*.rds' | wc -l)

if [[ "$CACHE_COUNT" -ne 108 ]]; then
  echo "Expected 108 cached DGM objects, found $CACHE_COUNT."
  echo "Run and review hpc/submit_cache.sh first."
  exit 2
fi

mkdir -p "$PROJECT_DIR/logs"

echo "Primary production submission"
echo "Project: $PROJECT_DIR"
echo "Account: $SLURM_ACCOUNT"
echo "Array tasks: $N_TASKS"
echo "Maximum concurrent tasks: $MAX_CONCURRENT"
echo "Cached DGM objects: $CACHE_COUNT"

JOB_ID=$(
  sbatch \
    --parsable \
    --account="$SLURM_ACCOUNT" \
    --chdir="$PROJECT_DIR" \
    --array="1-${N_TASKS}%${MAX_CONCURRENT}" \
    --export="ALL,SIM_PROJECT_DIR=$PROJECT_DIR,SIM_CALIBRATION_RDS=$CALIBRATION_RDS" \
    hpc/primary_array.sbatch
)

JOB_ID="${JOB_ID%%;*}"
echo "Submitted primary array: $JOB_ID"
echo "Monitor with: squeue -j $JOB_ID"
