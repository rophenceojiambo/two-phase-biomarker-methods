#!/bin/bash
#SBATCH --job-name=emp_resid_canary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/emp_resid_canary_%j.out
#SBATCH --error=logs/emp_resid_canary_%j.err

# Run the empirical residual-resampling canary under the same Torch environment as production.

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

Rscript -e 'cat("R version:", R.version.string, "\n"); cat("R_LIBS_USER:", Sys.getenv("R_LIBS_USER"), "\n"); cat("TMPDIR:", Sys.getenv("TMPDIR"), "\n"); cat("rsimsum:", as.character(packageVersion("rsimsum")), "\n")'

Rscript 23_run_empirical_residual_canary.R
