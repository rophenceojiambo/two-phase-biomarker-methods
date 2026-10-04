# Reproducibility notes

## What is in this repository

This repository contains the code and summarized results used for the manuscript. The analyses have already been completed. These notes show how the pieces of the project fit together and what is needed to rerun them.

## Analysis steps

1. **Prepare the MIDUS data:** `prepare_midus_analytic_sample.R`
2. **Set up and calibrate the simulation:** `00_config.R`, `01_calibrate_midus.R`
3. **Define the data-generating mechanism and methods:** `02_dgm.R`, `03_methods.R`, `04_build_dgm_cache.R`
4. **Run and summarize the primary simulation:** `05_run_primary_chunk.R` through `09_make_simulation_figures.R`
5. **Run the sensitivity analyses:** scripts `10` through `43`
6. **Run the MIDUS application:** scripts `44` through `48`
7. **Run jobs on Torch:** `hpc/`
8. **Check the IPW/AIPW sandwich variance calculations:** `validation/`

## MIDUS data needed

The MIDUS application uses:

```
data/MIDUS_discrimination_analysis.rds
```

If you have access to the data and keep the file elsewhere, set its path with `SIM_ANALYTIC_RDS`.

The MIDUS source data and participant-level analysis file are not included in this repository. `prepare_midus_analytic_sample.R` shows how the analysis file was created.

## Simulation results

The simulation produces several types of files, including intermediate chunks, combined results, summary measures, QC checks, tables, and figures.

The files kept in GitHub are the summarized results and QC files needed to understand the paper. Large intermediate chunks, checkpoints, and simulation caches are not included.

The main primary-simulation summaries are in:

- `results/summary/`
- `results/tables/`
- `results/qc/`
- `results/figure_data/`

Sensitivity-analysis results are kept in their own folders under `results/`.

Files used directly in the manuscript are collected in `results/manuscript/` so they are easy to find.

## MIDUS results

Scripts `44`–`48` produce the method estimates, sample summaries, weight checks, tables, and figures for the MIDUS application.

The aggregate results are in `results/real_world_application/`. Figure 6 and the data used for Tables 3–4 are also placed in `results/manuscript/main/`.

Participant-level MIDUS data are not included.

## Variance checks

Before the final production simulation, the custom IPW and AIPW sandwich standard errors were checked against:

- the same estimating equations with a numerical finite-difference Jacobian; and
- a participant-level nonparametric bootstrap.

The script and compact results are in `validation/`. These checks are separate from the final production simulation. The final simulation results in `results/summary/` are the results used for the manuscript.

## Software

The Torch analysis used R 4.5.1. Package requirements are listed in `hpc/install_required_packages.R`, and `sessionInfo.txt` records the R session used for the analysis.

## Before the final release

Remaining steps:

- finish the final repository review;
- update `results/PUBLIC_RESULTS_MANIFEST.csv` after the result files are final;
- update `CITATION.cff` with the final citation information and DOI;
- create the GitHub release and DOI.
