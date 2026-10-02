# HPC / Torch workflow

This directory contains SLURM submission scripts and other orchestration code used to run the simulation study on NYU Torch.

The HPC scripts should:

- assume the repository root as the project root;
- avoid hard-coded user-specific absolute paths where possible;
- use environment variables for site-specific settings such as project directories and SLURM accounts;
- write logs, temporary files, checkpoints, and chunk-level outputs to ignored locations;
- preserve reproducible random-number handling across simulation scenarios and parallel jobs.

Production simulation outputs should be combined and summarized before manuscript-facing results are copied to `results/`.
