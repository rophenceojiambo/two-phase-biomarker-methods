#!/bin/bash
# Print basic Torch account, storage, module, and scratch checks.

echo "============================================================"
echo "NYU Torch account / environment checks"
echo "============================================================"
echo

echo "1. Slurm accounts:"
my_slurm_accounts || true
echo

echo "2. Storage quota:"
myquota || true
echo

echo "3. Available R modules matching r/gcc:"
module avail r/gcc 2>&1 || true
echo

echo "4. Current directory:"
pwd
echo

echo "5. Scratch:"
echo "${SCRATCH:-/scratch/$USER}"
echo
