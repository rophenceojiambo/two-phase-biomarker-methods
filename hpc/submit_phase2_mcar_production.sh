#!/bin/bash
# Submit Phase-2 MCAR production, combine, and summary jobs.

set -euo pipefail

: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT before submitting.}"

PROJECT_DIR="${SIM_PROJECT_DIR:-$PWD}"
SLURM_QOS="${SLURM_QOS:-${SBATCH_QOS:-cpu48}}"
MAX_CONCURRENT="${MAX_CONCURRENT:-800}"
SIM_MCAR_NSIM="${SIM_MCAR_NSIM:-2000}"
SIM_MCAR_REPS_PER_CHUNK="${SIM_MCAR_REPS_PER_CHUNK:-10}"
SIM_MCAR_PRODUCTION_VALIDATION_N="${SIM_MCAR_PRODUCTION_VALIDATION_N:-100000}"

for VALUE_NAME in \
  MAX_CONCURRENT \
  SIM_MCAR_NSIM \
  SIM_MCAR_REPS_PER_CHUNK \
  SIM_MCAR_PRODUCTION_VALIDATION_N; do
  VALUE="${!VALUE_NAME}"
  if ! [[ "$VALUE" =~ ^[1-9][0-9]*$ ]]; then
    echo "$VALUE_NAME must be a positive integer."
    exit 2
  fi
done

if [[ "$SIM_MCAR_NSIM" -ne 2000 ]] || \
   [[ "$SIM_MCAR_REPS_PER_CHUNK" -ne 10 ]]; then
  echo "The prespecified production run requires SIM_MCAR_NSIM=2000 and SIM_MCAR_REPS_PER_CHUNK=10."
  exit 2
fi

if [[ "$SIM_MCAR_PRODUCTION_VALIDATION_N" -lt 100000 ]]; then
  echo "SIM_MCAR_PRODUCTION_VALIDATION_N must be at least 100000."
  exit 2
fi

CHUNKS_PER_SCENARIO=$((
  (SIM_MCAR_NSIM + SIM_MCAR_REPS_PER_CHUNK - 1) /
  SIM_MCAR_REPS_PER_CHUNK
))

N_SCENARIOS=4
N_TASKS=$((N_SCENARIOS * CHUNKS_PER_SCENARIO))

REQUIRED_FILES=(
  "00_config.R"
  "01_calibrate_midus.R"
  "02_dgm.R"
  "03_methods.R"
  "28_phase2_mcar_sensitivity_helpers.R"
  "29_run_phase2_mcar_canary.R"
  "30_build_phase2_mcar_production_cache.R"
  "31_run_phase2_mcar_production_chunk.R"
  "32_combine_phase2_mcar_results.R"
  "33_summarize_phase2_mcar_results.R"
  "hpc/phase2_mcar_production.sbatch"
  "hpc/phase2_mcar_combine.sbatch"
  "hpc/phase2_mcar_summary.sbatch"
)

for REQUIRED_FILE in "${REQUIRED_FILES[@]}"; do
  if [[ ! -f "$PROJECT_DIR/$REQUIRED_FILE" ]]; then
    echo "Required file not found: $PROJECT_DIR/$REQUIRED_FILE"
    exit 2
  fi
done

EMPIRICAL_PASS_FILE="$PROJECT_DIR/results/empirical_residual_sensitivity/production/summary/EMPIRICAL_RESIDUAL_PRODUCTION_SUMMARY_PASS.txt"
CANARY_PASS_FILE="$PROJECT_DIR/results/phase2_mcar_sensitivity/canary/PHASE2_MCAR_CANARY_PASS.txt"
CACHE_PASS_FILE="$PROJECT_DIR/results/phase2_mcar_sensitivity/production/cache/PHASE2_MCAR_PRODUCTION_CACHE_PASS.txt"

if [[ ! -f "$EMPIRICAL_PASS_FILE" ]] || \
   ! grep -q "PRODUCTION SUMMARY: PASS" "$EMPIRICAL_PASS_FILE"; then
  echo "Empirical-residual production summary PASS marker is missing or invalid:"
  echo "$EMPIRICAL_PASS_FILE"
  exit 2
fi

if [[ ! -f "$CANARY_PASS_FILE" ]] || \
   ! grep -q "MCAR CANARY: PASS" "$CANARY_PASS_FILE"; then
  echo "Phase-2 MCAR canary PASS marker is missing or invalid:"
  echo "$CANARY_PASS_FILE"
  exit 2
fi

if [[ ! -f "$CACHE_PASS_FILE" ]] || \
   ! grep -q "MCAR PRODUCTION CACHE: PASS" "$CACHE_PASS_FILE"; then
  echo "Phase-2 MCAR production cache PASS marker is missing or invalid:"
  echo "$CACHE_PASS_FILE"
  exit 2
fi

mkdir -p "$PROJECT_DIR/logs"

EXPORT_SETTINGS="ALL,SIM_PROJECT_DIR=$PROJECT_DIR,SIM_MCAR_NSIM=$SIM_MCAR_NSIM,SIM_MCAR_REPS_PER_CHUNK=$SIM_MCAR_REPS_PER_CHUNK,SIM_MCAR_PRODUCTION_VALIDATION_N=$SIM_MCAR_PRODUCTION_VALIDATION_N"

echo "Phase-2 MCAR sensitivity production submission"
echo "Project: $PROJECT_DIR"
echo "Account: $SLURM_ACCOUNT"
echo "QoS: $SLURM_QOS"
echo "Scenarios: $N_SCENARIOS"
echo "Repetitions per scenario: $SIM_MCAR_NSIM"
echo "Repetitions per array task: $SIM_MCAR_REPS_PER_CHUNK"
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
    hpc/phase2_mcar_production.sbatch
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
    hpc/phase2_mcar_combine.sbatch
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
    hpc/phase2_mcar_summary.sbatch
)

SUMMARY_JOB="${SUMMARY_JOB%%;*}"
echo "Submitted summary job: $SUMMARY_JOB"

echo
echo "Monitor the full chain with:"
echo "squeue -j ${PRODUCTION_JOB},${COMBINE_JOB},${SUMMARY_JOB}"
echo
echo "Combine and summary run only after their upstream jobs succeed."
