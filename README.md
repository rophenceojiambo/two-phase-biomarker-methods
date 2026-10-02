# Two-Phase Biomarker Methods

Reproducible code and manuscript-facing results for a methodological study of approaches to handling partially observed leukocyte covariates in a two-phase biomarker design, with an empirical application to MIDUS.

## Repository organization

The current analysis scripts are intentionally kept at the **repository root** as numbered files (`00_*.R` through `48_*.R`). The scripts source one another by filename and are designed to be run from the project root, so retaining that layout preserves the tested execution workflow.

- `00_*.R`–`09_*.R` — shared configuration, calibration, primary simulation, summaries, tables, and figures
- `10_*.R`–`43_*.R` — prespecified simulation sensitivity and robustness analyses
- `44_*.R`–`48_*.R` — empirical MIDUS application and manuscript outputs
- `hpc/` — SLURM/Torch submission and orchestration scripts
- `results/` — manuscript-facing summarized results, quality-control summaries, tables, and figures
- `manuscript/` — final manuscript tables and figures
- `data/` — documentation and synthetic/example data only; restricted MIDUS data are excluded
- `docs/` — simulation design and reproducibility documentation
- `R/` — reserved for future modularization/refactoring; the validated numbered scripts remain at the root for the manuscript reproducibility release

## Data availability

Participant-level MIDUS data are **not distributed in this repository**. Access to MIDUS data is governed by the applicable data-use conditions. This repository contains or will contain the analysis code, variable specifications, and allowable aggregate outputs needed to reproduce the empirical application for users with authorized data access.

The analysis scripts expect the authorized analytic dataset locally at `data/MIDUS_discrimination_analysis.rds` by default, or at a path supplied through the `SIM_ANALYTIC_RDS` environment variable. The dataset itself is ignored by Git and must not be committed.

Synthetic data and simulation code used for the methodological evaluation may be distributed here.

## Reproducibility workflow

The project is organized so that the analysis can be reproduced from a clean checkout:

1. Configure the project and authorized data path.
2. Run MIDUS calibration.
3. Build and check the primary data-generating mechanisms.
4. Run the primary simulation study.
5. Run prespecified simulation sensitivity analyses.
6. Combine repetition-level outputs.
7. Summarize performance measures and Monte Carlo uncertainty.
8. Generate manuscript figures and tables.
9. Run the authorized MIDUS empirical application.
10. Generate and check empirical manuscript outputs.

The full simulation is intended for high-performance computing. Local runs should be limited to calibration, canary, quality-control, and small-scale testing.

## Software

Analyses are implemented in R. A locked package environment will be added using `renv`.

## Public-repository safeguards

Restricted participant-level data, participant-level intermediate objects, raw SLURM outputs, scratch files, and large repetition-level simulation outputs are excluded through `.gitignore`. Before each public release, repository contents should be checked for accidental data files, local absolute paths, credentials, and disclosure-sensitive outputs.

## Status

Repository setup is in progress. The cleaned numbered scripts will be added before the clean Torch reproducibility rerun. Manuscript-facing results will be added only after that rerun and quality-control review.
