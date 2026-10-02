#!/bin/bash
#SBATCH --job-name=aipw_dr_canary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/aipw_dr_canary_%j.out
#SBATCH --error=logs/aipw_dr_canary_%j.err

# Run the AIPW double-robustness canary and its QC check.

set -euo pipefail

PROJECT_DIR="${SIM_PROJECT_DIR:-$PWD}"
R_MODULE="${R_MODULE:-r/4.5.1}"
R_LIBS_USER="${R_LIBS_USER:-${SCRATCH:-/scratch/$USER}/R/4.5.1/library}"
TMPDIR="${TMPDIR:-${SCRATCH:-/scratch/$USER}/tmp}"

export SIM_PROJECT_DIR="$PROJECT_DIR"
export R_MODULE
export R_LIBS_USER
export TMPDIR
export APPTAINER_TMPDIR="$TMPDIR"
export SINGULARITY_TMPDIR="$TMPDIR"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mkdir -p "$R_LIBS_USER" "$TMPDIR" "$PROJECT_DIR/logs"

cd "$PROJECT_DIR"

module purge
module load "$R_MODULE"

Rscript 11_run_aipw_robustness_canary.R
Rscript 12_check_aipw_robustness_canary.R
