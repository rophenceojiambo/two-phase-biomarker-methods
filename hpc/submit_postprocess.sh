#!/bin/bash
# Run primary content QC, combine validated chunks, summarize, and make figures.

set -euo pipefail

: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT before submitting.}"

PROJECT_DIR="${SIM_PROJECT_DIR:-$PWD}"

REQUIRED_FILES=(
  "06_combine_primary_results.R"
  "07_summarize_primary_rsimsum.R"
  "08_make_manuscript_tables.R"
  "09_make_simulation_figures.R"
  "hpc/03_validate_primary_content.R"
  "hpc/content_qc.sbatch"
  "hpc/combine.sbatch"
  "hpc/summary.sbatch"
)

for REQUIRED_FILE in "${REQUIRED_FILES[@]}"; do
  if [[ ! -f "$PROJECT_DIR/$REQUIRED_FILE" ]]; then
    echo "Required file not found: $PROJECT_DIR/$REQUIRED_FILE"
    exit 2
  fi
done

mkdir -p "$PROJECT_DIR/logs"

QC_JOB=$(
  sbatch \
    --parsable \
    --account="$SLURM_ACCOUNT" \
    --chdir="$PROJECT_DIR" \
    --export="ALL,SIM_PROJECT_DIR=$PROJECT_DIR" \
    hpc/content_qc.sbatch
)
QC_JOB="${QC_JOB%%;*}"
echo "Submitted primary content QC: $QC_JOB"

COMBINE_JOB=$(
  sbatch \
    --parsable \
    --dependency="afterok:${QC_JOB}" \
    --account="$SLURM_ACCOUNT" \
    --chdir="$PROJECT_DIR" \
    --export="ALL,SIM_PROJECT_DIR=$PROJECT_DIR" \
    hpc/combine.sbatch
)
COMBINE_JOB="${COMBINE_JOB%%;*}"
echo "Submitted combine job after QC: $COMBINE_JOB"

SUMMARY_JOB=$(
  sbatch \
    --parsable \
    --dependency="afterok:${COMBINE_JOB}" \
    --account="$SLURM_ACCOUNT" \
    --chdir="$PROJECT_DIR" \
    --export="ALL,SIM_PROJECT_DIR=$PROJECT_DIR" \
    hpc/summary.sbatch
)
SUMMARY_JOB="${SUMMARY_JOB%%;*}"
echo "Submitted summary/figure job after combine: $SUMMARY_JOB"

echo
echo "Monitor with:"
echo "squeue -j ${QC_JOB},${COMBINE_JOB},${SUMMARY_JOB}"
