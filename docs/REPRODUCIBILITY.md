# Reproducibility and archival guide

## Purpose

This repository is the code and results archive accompanying the manuscript. The analyses have already been completed; the purpose of the repository is to document the workflow and make the analysis code and allowable derived outputs available to readers.

## Analysis stages

1. **Configuration and calibration:** `00_config.R`, `01_calibrate_midus.R`
2. **Primary DGM and methods:** `02_dgm.R`, `03_methods.R`, `04_build_dgm_cache.R`
3. **Primary simulation:** `05_run_primary_chunk.R` through `09_make_simulation_figures.R`
4. **Sensitivity analyses:** scripts `10` through `43`
5. **MIDUS empirical application:** scripts `44` through `48`
6. **HPC execution and QC:** `hpc/`

## Required private input

The empirical application requires the analysis-ready MIDUS file:

```
data/MIDUS_discrimination_analysis.rds
```

This file is restricted and must remain outside GitHub.

The current public code does not yet include the upstream script that constructs this analysis-ready file from the authorized MIDUS source datasets. That script should be added before the final archival release if it is distributable.

## Simulation outputs

The simulation creates several classes of files:

- cached DGM objects;
- repetition/chunk results;
- combined scenario-level estimates;
- summarized performance measures;
- Monte Carlo uncertainty summaries;
- QC diagnostics;
- manuscript tables and figures.

Only the compact summaries, QC files supporting the manuscript, and final figure/table products need to be archived publicly. Raw repetition-level objects and caches are regenerable and unnecessarily large.

## Real-world outputs

Scripts `44`–`48` generate aggregate method estimates, sample summaries, weight diagnostics, manuscript tables, figures, and QC files. Only disclosure-safe aggregate outputs should be public.

## Software environment

The Torch scripts default to R 4.5.1. Package requirements are listed in `hpc/install_required_packages.R`. If the original analysis environment is still available, capture `sessionInfo()` or create an `renv.lock` file before the final release; this does not require rerunning the simulation.

## Final archival checklist

Before creating the manuscript release:

- add the shareable MIDUS data-preparation script;
- upload the final primary simulation summaries;
- upload manuscript-relevant sensitivity summaries/QC files;
- upload disclosure-safe real-world aggregate outputs;
- upload final figures and machine-readable tables;
- capture the software environment if available;
- confirm that no participant-level data or restricted metadata are tracked;
- update `CITATION.cff` with the final manuscript author list/version/DOI as appropriate;
- create a versioned GitHub release and archive it with a DOI service such as Zenodo.
