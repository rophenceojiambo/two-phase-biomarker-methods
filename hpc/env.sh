#!/bin/bash
# Common Torch environment for the MIDUS two-phase simulation project.
# Source after setting SIM_PROJECT_DIR when running commands interactively.

export R_MODULE="${R_MODULE:-r/4.5.1}"
export R_LIBS_USER="${R_LIBS_USER:-${SCRATCH:-/scratch/$USER}/R/4.5.1/library}"
export TMPDIR="${TMPDIR:-${SCRATCH:-/scratch/$USER}/tmp}"
export APPTAINER_TMPDIR="$TMPDIR"
export SINGULARITY_TMPDIR="$TMPDIR"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mkdir -p "$R_LIBS_USER" "$TMPDIR"
