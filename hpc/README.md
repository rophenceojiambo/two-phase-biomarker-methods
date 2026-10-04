# HPC / Torch workflow

This folder contains the SLURM scripts used to run the simulation study on NYU Torch.

## Setup

Run commands from the repository root and set the Torch Slurm account:

```bash
export SIM_PROJECT_DIR=$PWD
export SLURM_ACCOUNT=<your-torch-slurm-account>
```

The scripts use the Torch R module `r/4.5.1` by default and use scratch space for the R library and temporary files. These settings can be changed with `R_MODULE`, `R_LIBS_USER`, and `TMPDIR`.

The MIDUS analysis file is not included. By default, the scripts look for:

```
data/MIDUS_discrimination_analysis.rds
```

Use `SIM_ANALYTIC_RDS` to set source of the file.

## Running the workflow again

Run the workflow in stages:

1. `bash hpc/preflight_rerun.sh`
2. `bash hpc/setup_project.sh`
3. Install/check the required R packages on a compute node with `hpc/install_required_packages.R`.
4. `bash hpc/submit_calibrate.sh`
5. `bash hpc/submit_smoke.sh`
6. `bash hpc/submit_cache.sh`
7. `bash hpc/submit_primary.sh`
8. Check for incomplete primary tasks and resubmit them if needed with `hpc/resubmit_missing.sh`.
9. `bash hpc/submit_postprocess.sh`
10. Run the sensitivity analyses in their documented order, including the canary/QC checks.

Do not move to the next stage if a required canary, cache, combination, summary, or PASS/QC check fails.

## Files saved to GitHub

SLURM logs, temporary files, checkpoints, chunk-level simulation results, and participant-level MIDUS data are not included.

After jobs finish, the combined and summarized results used for the paper can be added to `results/`.

## Reproducibility checks

The HPC scripts:

- use the repository root as the project directory;
- pass project paths through environment variables;
- use scratch space for R packages and temporary files;
- preserve the simulation random-number setup;
- check required files and PASS markers before later jobs run.
