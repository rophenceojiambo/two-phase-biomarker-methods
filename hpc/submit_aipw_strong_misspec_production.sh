#!/bin/bash
# Submit strong AIPW misspecification production and postprocessing jobs.
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
    hpc/aipw_strong_misspec_production.sbatch
)

echo "Submitted strong AIPW production array: $PROD_JOB"

POST_JOB=$(
  sbatch \
    --parsable \
    --dependency="afterok:${PROD_JOB}" \
    --account="$SLURM_ACCOUNT" \
    --qos="$QOS" \
    --chdir="$PROJECT_DIR" \
    hpc/aipw_strong_misspec_postprocess.sbatch
)

echo "Submitted strong AIPW postprocess after production: $POST_JOB"
