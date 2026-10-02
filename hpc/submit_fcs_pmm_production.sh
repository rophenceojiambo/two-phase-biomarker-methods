#!/bin/bash
# Submit FCS-PMM production, combine, and summary jobs.

set -euo pipefail

: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT before submitting.}"

PROJECT_DIR="${SIM_PROJECT_DIR:-$PWD}"
SLURM_QOS="${SLURM_QOS:-${SBATCH_QOS:-cpu48}}"
MAX_CONCURRENT="${MAX_CONCURRENT:-800}"
SIM_FCS_PMM_NSIM="${SIM_FCS_PMM_NSIM:-2000}"
SIM_FCS_PMM_REPS_PER_CHUNK="${SIM_FCS_PMM_REPS_PER_CHUNK:-10}"

for VALUE_NAME in \
  MAX_CONCURRENT \
  SIM_FCS_PMM_NSIM \
  SIM_FCS_PMM_REPS_PER_CHUNK; do
  VALUE="${!VALUE_NAME}"
  if ! [[ "$VALUE" =~ ^[1-9][0-9]*$ ]]; then
    echo "$VALUE_NAME must be a positive integer."
    exit 2
  fi
done

if [[ "$SIM_FCS_PMM_NSIM" -ne 2000 ]] || \
   [[ "$SIM_FCS_PMM_REPS_PER_CHUNK" -ne 10 ]]; then
  echo "The prespecified production run requires SIM_FCS_PMM_NSIM=2000 and SIM_FCS_PMM_REPS_PER_CHUNK=10."
  exit 2
fi

CHUNKS_PER_SCENARIO=$((
  (SIM_FCS_PMM_NSIM + SIM_FCS_PMM_REPS_PER_CHUNK - 1) /
  SIM_FCS_PMM_REPS_PER_CHUNK
))

N_SCENARIOS=4
N_TASKS=$((N_SCENARIOS * CHUNKS_PER_SCENARIO))

REQUIRED_FILES=(
  "00_config.R"
  "02_dgm.R"
  "03_methods.R"
  "34_fcs_pmm_sensitivity_helpers.R"
  "35_run_fcs_pmm_canary.R"
  "36_run_fcs_pmm_production_chunk.R"
  "37_combine_fcs_pmm_results.R"
  "38_summarize_fcs_pmm_results.R"
  "hpc/fcs_pmm_production.sbatch"
  "hpc/fcs_pmm_combine.sbatch"
  "hpc/fcs_pmm_summary.sbatch"
)

for REQUIRED_FILE in "${REQUIRED_FILES[@]}"; do
  if [[ ! -f "$PROJECT_DIR/$REQUIRED_FILE" ]]; then
    echo "Required file not found: $PROJECT_DIR/$REQUIRED_FILE"
    exit 2
  fi
done

MCAR_PASS_FILE="$PROJECT_DIR/results/phase2_mcar_sensitivity/production/summary/PHASE2_MCAR_PRODUCTION_SUMMARY_PASS.txt"
CANARY_PASS_FILE="$PROJECT_DIR/results/fcs_pmm_sensitivity/canary/FCS_PMM_CANARY_PASS.txt"
CANARY_QC_FILE="$PROJECT_DIR/results/fcs_pmm_sensitivity/canary/fcs_pmm_canary_qc.csv"

if [[ ! -f "$MCAR_PASS_FILE" ]] || \
   ! grep -q "MCAR PRODUCTION SUMMARY: PASS" "$MCAR_PASS_FILE"; then
  echo "Phase-2 MCAR production summary PASS marker is missing or invalid:"
  echo "$MCAR_PASS_FILE"
  exit 2
fi

if [[ ! -f "$CANARY_PASS_FILE" ]] || \
   ! grep -q "FCS-PMM CANARY: PASS" "$CANARY_PASS_FILE"; then
  echo "FCS-PMM canary PASS marker is missing or invalid:"
  echo "$CANARY_PASS_FILE"
  exit 2
fi

if [[ ! -f "$CANARY_QC_FILE" ]]; then
  echo "FCS-PMM canary QC file is missing:"
  echo "$CANARY_QC_FILE"
  exit 2
fi

mkdir -p "$PROJECT_DIR/logs"

EXPORT_SETTINGS="ALL,SIM_PROJECT_DIR=$PROJECT_DIR,SIM_FCS_PMM_NSIM=$SIM_FCS_PMM_NSIM,SIM_FCS_PMM_REPS_PER_CHUNK=$SIM_FCS_PMM_REPS_PER_CHUNK"

echo "FCS-PMM sensitivity production submission"
echo "Project: $PROJECT_DIR"
echo "Account: $SLURM_ACCOUNT"
echo "QoS: $SLURM_QOS"
echo "Scenarios: $N_SCENARIOS"
echo "Repetitions per scenario: $SIM_FCS_PMM_NSIM"
echo "Repetitions per array task: $SIM_FCS_PMM_REPS_PER_CHUNK"
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
    hpc/fcs_pmm_production.sbatch
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
    hpc/fcs_pmm_combine.sbatch
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
    hpc/fcs_pmm_summary.sbatch
)

SUMMARY_JOB="${SUMMARY_JOB%%;*}"
echo "Submitted summary job: $SUMMARY_JOB"

echo
echo "Monitor the full chain with:"
echo "squeue -j ${PRODUCTION_JOB},${COMBINE_JOB},${SUMMARY_JOB}"
echo
echo "Combine and summary run only after their upstream jobs succeed."
