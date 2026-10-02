# HPC / Torch workflow

This directory contains the SLURM submission scripts and orchestration code used to run the simulation study on NYU Torch.

## Environment

Run commands from the repository root and set the site-specific Slurm account explicitly:

```bash
export SIM_PROJECT_DIR=$PWD
export SLURM_ACCOUNT=<your-torch-slurm-account>
```

The scripts default to the Torch R module `r/4.5.1` and use the user's scratch area for the R library and temporary files. These defaults can be overridden with `R_MODULE`, `R_LIBS_USER`, and `TMPDIR`.

The analysis-ready MIDUS file is not stored in Git. By default, the workflow expects:

```
data/MIDUS_discrimination_analysis.rds
```

An authorized alternate location may be supplied with `SIM_ANALYTIC_RDS`.

## Clean-rerun sequence

A new reproducibility run should proceed in stages rather than submitting the full simulation immediately:

1. `bash hpc/preflight_rerun.sh`
2. `bash hpc/setup_project.sh`
3. Install/verify required R packages on a compute node using `hpc/install_required_packages.R`.
4. `bash hpc/submit_calibrate.sh`
5. `bash hpc/submit_smoke.sh`
6. `bash hpc/submit_cache.sh`
7. `bash hpc/submit_primary.sh`
8. Check/resubmit incomplete primary tasks if necessary using `hpc/resubmit_missing.sh`.
9. `bash hpc/submit_postprocess.sh`
10. Run the prespecified sensitivity analyses and their canary/QC steps in the documented dependency order.

Do not advance past a canary, cache, combination, or summary stage unless its corresponding PASS/QC checks succeed.

## Output policy

HPC logs, temporary files, checkpoints, chunk-level results, and participant-level data are not intended for GitHub. Production outputs should first be combined and summarized; only disclosure-safe manuscript-facing summaries, tables, figures, and QC outputs should be considered for the public repository.

## Reproducibility safeguards

The scripts are designed to:

- use the repository root as the project directory;
- avoid user-specific hard-coded absolute paths;
- pass project paths through environment variables;
- keep R package and temporary storage in scratch;
- preserve the simulation random-number workflow;
- validate expected files and PASS markers before downstream production stages.
