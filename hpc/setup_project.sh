#!/bin/bash
# Initialize directories and paths for a clean MIDUS two-phase simulation rerun.
#
# This script is safe to run on a Torch login node.
# It does NOT load or execute R because Torch's container-backed R module
# should be used only on an allocated compute node.

set -euo pipefail

PROJECT_DIR="${SIM_PROJECT_DIR:-$PWD}"
export SIM_PROJECT_DIR="$PROJECT_DIR"

cd "$PROJECT_DIR"

# Load shared project environment variables.
source hpc/env.sh

# Create directories used by calibration, production, and postprocessing.
mkdir -p \
  data \
  logs \
  calibration/dgm_cache \
  results/chunks \
  results/combined \
  results/summary \
  results/figures \
  results/tables

# Required analysis-ready MIDUS input.
ANALYTIC_RDS="${SIM_ANALYTIC_RDS:-$PROJECT_DIR/data/MIDUS_discrimination_analysis.rds}"

if [[ ! -f "$ANALYTIC_RDS" ]]; then
  echo "ERROR: Required analysis input not found:"
  echo "  $ANALYTIC_RDS"
  exit 1
fi

echo "MIDUS two-phase project setup complete."
echo
echo "Project directory: $PROJECT_DIR"
echo "R module:          $R_MODULE"
echo "R user library:    $R_LIBS_USER"
echo "Temporary files:   $TMPDIR"
echo
echo "Required analysis input:"
echo "  $ANALYTIC_RDS"
echo
echo "PASS: required analysis input found."
echo
echo "R is not executed by this script."
echo "Verify R and packages on a Torch compute node using:"
echo "  Rscript hpc/install_required_packages.R"
