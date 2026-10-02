#!/bin/bash
#SBATCH --job-name=ipw_cap_prod
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --output=logs/ipw_cap_prod_%A_%a.out
#SBATCH --error=logs/ipw_cap_prod_%A_%a.err

# Run one IPW weight-capping sensitivity production chunk.

set -euo pipefail

PROJECT_DIR="${SIM_PROJECT_DIR:-$PWD}"
R_MODULE="${R_MODULE:-r/4.5.1}"
R_LIBS_USER="${R_LIBS_USER:-${SCRATCH:-/scratch/$USER}/R/4.5.1/library}"
TMPDIR="${SCRATCH:-/scratch/$USER}/tmp"

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

Rscript 41_run_ipw_weight_capping_production_chunk.R
