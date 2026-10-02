#!/bin/bash
# Submit IPW weight-capping production, combine, and summary jobs.

set -euo pipefail

: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT before submitting.}"

PROJECT_DIR="${SIM_PROJECT_DIR:-$PWD}"
SLURM_QOS="${SLURM_QOS:-${SBATCH_QOS:-cpu48}}"
MAX_CONCURRENT="${MAX_CONCURRENT:-80}"
SIM_IPW_CAP_NSIM="${SIM_IPW_CAP_NSIM:-2000}"
SIM_IPW_CAP_REPS_PER_CHUNK="${SIM_IPW_CAP_REPS_PER_CHUNK:-100}"

for VALUE_NAME in \
  MAX_CONCURRENT \
  SIM_IPW_CAP_NSIM \
  SIM_IPW_CAP_REPS_PER_CHUNK; do
  VALUE="${!VALUE_NAME}"
  if ! [[ "$VALUE" =~ ^[1-9][0-9]*$ ]]; then
    echo "$VALUE_NAME must be a positive integer."
    exit 2
  fi
done

if [[ "$SIM_IPW_CAP_NSIM" -ne 2000 ]] || \
   [[ "$SIM_IPW_CAP_REPS_PER_CHUNK" -ne 100 ]]; then
  echo "The prespecified run requires SIM_IPW_CAP_NSIM=2000 and SIM_IPW_CAP_REPS_PER_CHUNK=100."
  exit 2
fi

CHUNKS_PER_SCENARIO=$((
  (SIM_IPW_CAP_NSIM + SIM_IPW_CAP_REPS_PER_CHUNK - 1) /
  SIM_IPW_CAP_REPS_PER_CHUNK
))

N_SCENARIOS=4
N_TASKS=$((N_SCENARIOS * CHUNKS_PER_SCENARIO))

REQUIRED_FILES=(
  "00_config.R"
  "02_dgm.R"
  "03_methods.R"
  "39_ipw_weight_capping_helpers.R"
  "40_run_ipw_weight_capping_canary.R"
  "41_run_ipw_weight_capping_production_chunk.R"
  "42_combine_ipw_weight_capping_results.R"
  "43_summarize_ipw_weight_capping_results.R"
  "hpc/ipw_weight_capping_production.sbatch"
  "hpc/ipw_weight_capping_combine.sbatch"
  "hpc/ipw_weight_capping_summary.sbatch"
)

for REQUIRED_FILE in "${REQUIRED_FILES[@]}"; do
  if [[ ! -f "$PROJECT_DIR/$REQUIRED_FILE" ]]; then
    echo "Required file not found: $PROJECT_DIR/$REQUIRED_FILE"
    exit 2
  fi
done

FCS_PMM_PASS_FILE="$PROJECT_DIR/results/fcs_pmm_sensitivity/production/summary/FCS_PMM_PRODUCTION_SUMMARY_PASS.txt"
CANARY_PASS_FILE="$PROJECT_DIR/results/ipw_weight_capping_sensitivity/canary/IPW_WEIGHT_CAPPING_CANARY_PASS.txt"
CANARY_QC_FILE="$PROJECT_DIR/results/ipw_weight_capping_sensitivity/canary/ipw_weight_capping_canary_qc.csv"

if [[ ! -f "$FCS_PMM_PASS_FILE" ]] || \
   ! grep -q "FCS-PMM PRODUCTION SUMMARY: PASS" "$FCS_PMM_PASS_FILE"; then
  echo "FCS-PMM production summary PASS marker is missing or invalid:"
  echo "$FCS_PMM_PASS_FILE"
  exit 2
fi

if [[ ! -f "$CANARY_PASS_FILE" ]] || \
   ! grep -q "IPW WEIGHT CAPPING CANARY: PASS" "$CANARY_PASS_FILE"; then
  echo "IPW weight-capping canary PASS marker is missing or invalid:"
  echo "$CANARY_PASS_FILE"
  exit 2
fi

if [[ ! -f "$CANARY_QC_FILE" ]]; then
  echo "IPW weight-capping canary QC file is missing:"
  echo "$CANARY_QC_FILE"
  exit 2
fi

mkdir -p "$PROJECT_DIR/logs"

EXPORT_SETTINGS="ALL,SIM_PROJECT_DIR=$PROJECT_DIR,SIM_IPW_CAP_NSIM=$SIM_IPW_CAP_NSIM,SIM_IPW_CAP_REPS_PER_CHUNK=$SIM_IPW_CAP_REPS_PER_CHUNK"

echo "IPW weight-capping sensitivity production submission"
echo "Project: $PROJECT_DIR"
echo "Account: $SLURM_ACCOUNT"
echo "QoS: $SLURM_QOS"
echo "Scenarios: $N_SCENARIOS"
echo "Repetitions per scenario: $SIM_IPW_CAP_NSIM"
echo "Repetitions per array task: $SIM_IPW_CAP_REPS_PER_CHUNK"
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
    hpc/ipw_weight_capping_production.sbatch
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
    hpc/ipw_weight_capping_combine.sbatch
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
    hpc/ipw_weight_capping_summary.sbatch
)

SUMMARY_JOB="${SUMMARY_JOB%%;*}"
echo "Submitted summary job: $SUMMARY_JOB"

echo
echo "Monitor the full chain with:"
echo "squeue -j ${PRODUCTION_JOB},${COMBINE_JOB},${SUMMARY_JOB}"
echo
echo "Combine and summary run only after their upstream jobs succeed."
