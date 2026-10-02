#!/bin/bash
#SBATCH --job-name=midus_cache
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/cache_%j.out
#SBATCH --error=logs/cache_%j.err

# Build the 108 cached DGM parameter objects used by the primary simulation.

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

Rscript 04_build_dgm_cache.R
