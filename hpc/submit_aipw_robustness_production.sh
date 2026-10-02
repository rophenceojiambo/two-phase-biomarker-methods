#!/bin/bash
# Submit AIPW double-robustness production and postprocessing jobs.
set -euo pipefail

: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT first.}"
PROJECT_DIR="${SIM_PROJECT_DIR:-$PWD}"
QOS="${SLURM_QOS:-cpu48}"

PROD_JOB=$(
  sbatch \
    --parsable \
    --account="$SLURM_ACCOUNT" \
    --qos="$QOS" \
    --chdir="$PROJECT_DIR" \
    hpc/aipw_robustness_production.sbatch
)

echo "Submitted AIPW robustness production array: $PROD_JOB"

POST_JOB=$(
  sbatch \
    --parsable \
    --dependency="afterok:${PROD_JOB}" \
    --account="$SLURM_ACCOUNT" \
    --qos="$QOS" \
    --chdir="$PROJECT_DIR" \
    hpc/aipw_robustness_postprocess.sbatch
)

echo "Submitted postprocess after successful production: $POST_JOB"
