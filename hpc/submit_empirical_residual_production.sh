#!/bin/bash
# Submit empirical residual-resampling production, combine, and summary jobs.

set -euo pipefail

: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT before submitting.}"

PROJECT_DIR="${SIM_PROJECT_DIR:-$PWD}"
SLURM_QOS="${SLURM_QOS:-${SBATCH_QOS:-cpu48}}"
MAX_CONCURRENT="${MAX_CONCURRENT:-800}"
SIM_EMPIRICAL_NSIM="${SIM_EMPIRICAL_NSIM:-2000}"
SIM_EMPIRICAL_REPS_PER_CHUNK="${SIM_EMPIRICAL_REPS_PER_CHUNK:-10}"
SIM_EMPIRICAL_PRODUCTION_VALIDATION_N="${SIM_EMPIRICAL_PRODUCTION_VALIDATION_N:-100000}"
SIM_EMPIRICAL_SELECTION_CALIBRATION_N="${SIM_EMPIRICAL_SELECTION_CALIBRATION_N:-50000}"

for VALUE_NAME in \
  MAX_CONCURRENT \
  SIM_EMPIRICAL_NSIM \
  SIM_EMPIRICAL_REPS_PER_CHUNK \
  SIM_EMPIRICAL_PRODUCTION_VALIDATION_N \
  SIM_EMPIRICAL_SELECTION_CALIBRATION_N; do
  VALUE="${!VALUE_NAME}"
  if ! [[ "$VALUE" =~ ^[1-9][0-9]*$ ]]; then
    echo "$VALUE_NAME must be a positive integer."
    exit 2
  fi
done

if [[ "$SIM_EMPIRICAL_NSIM" -ne 2000 ]] || \
   [[ "$SIM_EMPIRICAL_REPS_PER_CHUNK" -ne 10 ]]; then
  echo "The prespecified production run requires SIM_EMPIRICAL_NSIM=2000 and SIM_EMPIRICAL_REPS_PER_CHUNK=10."
  exit 2
fi

if [[ "$SIM_EMPIRICAL_PRODUCTION_VALIDATION_N" -lt 100000 ]]; then
  echo "SIM_EMPIRICAL_PRODUCTION_VALIDATION_N must be at least 100000."
  exit 2
fi

CHUNKS_PER_SCENARIO=$((
  (SIM_EMPIRICAL_NSIM + SIM_EMPIRICAL_REPS_PER_CHUNK - 1) /
  SIM_EMPIRICAL_REPS_PER_CHUNK
))

N_SCENARIOS=4
N_TASKS=$((N_SCENARIOS * CHUNKS_PER_SCENARIO))

REQUIRED_FILES=(
  "00_config.R"
  "02_dgm.R"
  "03_methods.R"
  "22_empirical_residual_sensitivity_helpers.R"
  "23_run_empirical_residual_canary.R"
  "24_build_empirical_residual_production_cache.R"
  "25_run_empirical_residual_production_chunk.R"
  "26_combine_empirical_residual_production.R"
  "27_summarize_empirical_residual_production.R"
  "hpc/empirical_residual_production.sbatch"
  "hpc/empirical_residual_combine.sbatch"
  "hpc/empirical_residual_summary.sbatch"
)

for REQUIRED_FILE in "${REQUIRED_FILES[@]}"; do
  if [[ ! -f "$PROJECT_DIR/$REQUIRED_FILE" ]]; then
    echo "Required file not found: $PROJECT_DIR/$REQUIRED_FILE"
    exit 2
  fi
done

CANARY_PASS_FILE="$PROJECT_DIR/results/empirical_residual_sensitivity/canary/EMPIRICAL_RESIDUAL_CANARY_PASS.txt"
CACHE_PASS_FILE="$PROJECT_DIR/results/empirical_residual_sensitivity/production/cache/EMPIRICAL_RESIDUAL_PRODUCTION_CACHE_PASS.txt"

if [[ ! -f "$CANARY_PASS_FILE" ]] || \
   ! grep -q "CANARY: PASS" "$CANARY_PASS_FILE"; then
  echo "Empirical-residual canary PASS marker is missing or invalid:"
  echo "$CANARY_PASS_FILE"
  exit 2
fi

if [[ ! -f "$CACHE_PASS_FILE" ]] || \
   ! grep -q "PRODUCTION CACHE: PASS" "$CACHE_PASS_FILE"; then
  echo "Empirical-residual production cache PASS marker is missing or invalid:"
  echo "$CACHE_PASS_FILE"
  exit 2
fi

mkdir -p "$PROJECT_DIR/logs"

EXPORT_SETTINGS="ALL,SIM_PROJECT_DIR=$PROJECT_DIR,SIM_EMPIRICAL_NSIM=$SIM_EMPIRICAL_NSIM,SIM_EMPIRICAL_REPS_PER_CHUNK=$SIM_EMPIRICAL_REPS_PER_CHUNK,SIM_EMPIRICAL_PRODUCTION_VALIDATION_N=$SIM_EMPIRICAL_PRODUCTION_VALIDATION_N,SIM_EMPIRICAL_SELECTION_CALIBRATION_N=$SIM_EMPIRICAL_SELECTION_CALIBRATION_N"

echo "Empirical residual-resampling production submission"
echo "Project: $PROJECT_DIR"
echo "Account: $SLURM_ACCOUNT"
echo "QoS: $SLURM_QOS"
echo "Scenarios: $N_SCENARIOS"
echo "Repetitions per scenario: $SIM_EMPIRICAL_NSIM"
echo "Repetitions per array task: $SIM_EMPIRICAL_REPS_PER_CHUNK"
echo "Array tasks: $N_TASKS"
echo "Maximum concurrent tasks: $MAX_CONCURRENT"

PRODUCTION_JOB=$(
  sbatch \
    --parsable \
    --account="$SLURM_ACCOUNT" \
    --qos="$SLURM_QOS" \
    --chdir="$PROJECT_DIR" \
    --array="1-${N_TASKS}%${MAX_CONCURRENT}" \
    --export="$EXPORT_SETTINGS" \
    hpc/empirical_residual_production.sbatch
)
PRODUCTION_JOB="${PRODUCTION_JOB%%;*}"
echo "Submitted production array: $PRODUCTION_JOB"

COMBINE_JOB=$(
  sbatch \
    --parsable \
    --account="$SLURM_ACCOUNT" \
    --qos="$SLURM_QOS" \
    --chdir="$PROJECT_DIR" \
    --dependency="afterok:${PRODUCTION_JOB}" \
    --export="$EXPORT_SETTINGS" \
    hpc/empirical_residual_combine.sbatch
)
COMBINE_JOB="${COMBINE_JOB%%;*}"
echo "Submitted combine job: $COMBINE_JOB"

SUMMARY_JOB=$(
  sbatch \
    --parsable \
    --account="$SLURM_ACCOUNT" \
    --qos="$SLURM_QOS" \
    --chdir="$PROJECT_DIR" \
    --dependency="afterok:${COMBINE_JOB}" \
    --export="$EXPORT_SETTINGS" \
    hpc/empirical_residual_summary.sbatch
)
SUMMARY_JOB="${SUMMARY_JOB%%;*}"
echo "Submitted summary job: $SUMMARY_JOB"

echo
echo "Monitor with:"
echo "squeue -j ${PRODUCTION_JOB},${COMBINE_JOB},${SUMMARY_JOB}"
