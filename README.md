# Two-Phase Biomarker Methods

Reproducible code and manuscript-facing results for a methodological study of approaches to handling partially observed leukocyte covariates in a two-phase biomarker design, with an empirical application to MIDUS.

## Repository organization

- `R/simulation/` — simulation data generation, estimators, sensitivity analyses, summaries, and figures
- `R/application/` — empirical MIDUS application code
- `hpc/` — SLURM/Torch submission and orchestration scripts
- `results/` — manuscript-facing summarized results, quality-control summaries, tables, and figures
- `manuscript/` — final manuscript tables and figures
- `data/` — documentation and synthetic/example data only
- `docs/` — simulation design and reproducibility documentation

## Data availability

Participant-level MIDUS data are **not distributed in this repository**. Access to MIDUS data is governed by the applicable data-use conditions. This repository will contain the analysis code, variable specifications, and allowable aggregate outputs needed to reproduce the empirical application for users with authorized data access.

Synthetic data and simulation code used for the methodological evaluation may be distributed here.

## Reproducibility workflow

The project is organized so that the analysis can be reproduced from a clean checkout:

1. Run simulation calibration and quality-control checks.
2. Run the primary simulation study.
3. Run prespecified simulation sensitivity analyses.
4. Combine repetition-level outputs.
5. Summarize performance measures and Monte Carlo uncertainty.
6. Generate manuscript figures and tables.
7. Prepare the authorized MIDUS analytic data locally.
8. Run the empirical application.
9. Generate empirical tables and figures.

The full simulation is intended for high-performance computing. Local runs should be limited to calibration, canary, and small-scale testing.

## Software

Analyses are implemented in R. A locked package environment will be added using `renv`.

## Status

Repository setup is in progress. Scripts and manuscript-facing outputs will be added after the clean reproducibility rerun.
