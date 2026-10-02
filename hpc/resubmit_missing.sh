#!/bin/bash
# Check the primary array on a Torch compute node and resubmit only missing tasks.
#
# This wrapper is intended to be launched from a Torch login node. It does NOT
# run R on the login node. Instead, it submits the completion checker as a small
# Slurm job, waits for that checker to finish, then reads the missing-task list.

set -euo pipefail

: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT before submitting.}"

PROJECT_DIR="${SIM_PROJECT_DIR:-$PWD}"
SLURM_QOS="${SLURM_QOS:-cpu48}"
MAX_CONCURRENT="${MAX_CONCURRENT:-2000}"
R_MODULE="${R_MODULE:-r/4.5.1}"
R_LIBS_USER="${R_LIBS_USER:-${SCRATCH:-/scratch/$USER}/R/4.5.1/library}"
SCRATCH_DIR="${SCRATCH:-/scratch/$USER}"
CHECK_TMPDIR="$SCRATCH_DIR/tmp"
MISSING_FILE="$PROJECT_DIR/results/primary_missing_array_ids.txt"

if ! [[ "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]]; then
  echo "MAX_CONCURRENT must be a positive integer."
  exit 2
fi

if [[ ! -d "$PROJECT_DIR" ]]; then
  echo "Project directory not found: $PROJECT_DIR"
  exit 2
fi

for REQUIRED_FILE in \
  "hpc/env.sh" \
  "hpc/02_check_primary_completion.R" \
  "hpc/primary_array.sbatch"; do
  if [[ ! -f "$PROJECT_DIR/$REQUIRED_FILE" ]]; then
    echo "Required file not found: $PROJECT_DIR/$REQUIRED_FILE"
    exit 2
  fi
done

export SIM_PROJECT_DIR="$PROJECT_DIR"
export SLURM_QOS
export R_MODULE
export R_LIBS_USER

cd "$PROJECT_DIR"
mkdir -p logs results "$CHECK_TMPDIR"

CHECK_OUT="$PROJECT_DIR/logs/primary_completion_%j.out"
CHECK_ERR="$PROJECT_DIR/logs/primary_completion_%j.err"

echo "Primary completion check"
echo "Project: $PROJECT_DIR"
echo "Account: $SLURM_ACCOUNT"
echo "QOS:     $SLURM_QOS"
echo

echo "Submitting the completion checker to a compute node..."

# Use sbatch --wait so this login-node wrapper continues only after the R checker
# has finished on a compute node. Explicit scratch temp paths prevent the
# container-backed R module from inheriting an unusable login-node TMPDIR.
set +e
CHECK_JOB_RAW=$(
  sbatch \
    --parsable \
    --wait \
    --account="$SLURM_ACCOUNT" \
    --qos="$SLURM_QOS" \
    --chdir="$PROJECT_DIR" \
    --job-name=primary_completion_check \
    --time=00:10:00 \
    --mem=1G \
    --output="$CHECK_OUT" \
    --error="$CHECK_ERR" \
    --export="ALL,SIM_PROJECT_DIR=$PROJECT_DIR,R_MODULE=$R_MODULE,R_LIBS_USER=$R_LIBS_USER,TMPDIR=$CHECK_TMPDIR,APPTAINER_TMPDIR=$CHECK_TMPDIR,SINGULARITY_TMPDIR=$CHECK_TMPDIR" \
    --wrap='source hpc/env.sh; module purge; module load "$R_MODULE"; Rscript hpc/02_check_primary_completion.R'
)
CHECK_STATUS=$?
set -e

CHECK_JOB_ID="${CHECK_JOB_RAW%%;*}"

if [[ -z "$CHECK_JOB_ID" ]]; then
  echo "Could not determine the completion-check Slurm job ID."
  exit 1
fi

echo "Completion-check job: $CHECK_JOB_ID"
echo "Stdout: logs/primary_completion_${CHECK_JOB_ID}.out"
echo "Stderr: logs/primary_completion_${CHECK_JOB_ID}.err"

if [[ "$CHECK_STATUS" -ne 0 ]]; then
  echo
  echo "The completion-check job failed (sbatch --wait exit status $CHECK_STATUS)."
  echo "Inspect:"
  echo "  cat logs/primary_completion_${CHECK_JOB_ID}.out"
  echo "  cat logs/primary_completion_${CHECK_JOB_ID}.err"
  exit "$CHECK_STATUS"
fi

if [[ ! -f "$MISSING_FILE" ]]; then
  echo "Completion checker did not create: $MISSING_FILE"
  exit 1
fi

MISSING_IDS="$(tr -d '[:space:]' < "$MISSING_FILE")"

if [[ -z "$MISSING_IDS" ]]; then
  echo
  echo "No missing primary array tasks."
  echo "Primary completion check passed; nothing was resubmitted."
  exit 0
fi

# Count comma-separated task IDs for an informative message.
MISSING_COUNT=$(
  awk -F',' '{print NF}' <<< "$MISSING_IDS"
)

echo
echo "Missing primary array tasks: $MISSING_COUNT"
echo "Resubmitting missing task IDs:"
echo "$MISSING_IDS"

JOB_ID=$(
  sbatch \
    --parsable \
    --account="$SLURM_ACCOUNT" \
    --qos="$SLURM_QOS" \
    --chdir="$PROJECT_DIR" \
    --array="${MISSING_IDS}%${MAX_CONCURRENT}" \
    --export="ALL,SIM_PROJECT_DIR=$PROJECT_DIR" \
    hpc/primary_array.sbatch
)

JOB_ID="${JOB_ID%%;*}"

echo
echo "Submitted missing-task array: $JOB_ID"
echo "Maximum concurrent tasks: $MAX_CONCURRENT"
echo "Monitor with: squeue -j $JOB_ID"
echo
echo "After that array finishes, run this script again."
