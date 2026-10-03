# Reproducibility and archival guide

## Purpose

This repository is the code and results archive accompanying the manuscript. The analyses have already been completed; the purpose of the repository is to document the workflow and make the analysis code and allowable derived outputs available to readers.

## Analysis stages

1. **MIDUS data preparation:** `prepare_midus_analytic_sample.R`
2. **Configuration and calibration:** `00_config.R`, `01_calibrate_midus.R`
3. **Primary DGM and methods:** `02_dgm.R`, `03_methods.R`, `04_build_dgm_cache.R`
4. **Primary simulation:** `05_run_primary_chunk.R` through `09_make_simulation_figures.R`
5. **Sensitivity analyses:** scripts `10` through `43`
6. **MIDUS empirical application:** scripts `44` through `48`
7. **HPC execution and QC:** `hpc/`

## Required private input

The empirical application requires the analysis-ready MIDUS file:

```
data/MIDUS_discrimination_analysis.rds
```

This file is restricted and must remain outside GitHub.

The public script `prepare_midus_analytic_sample.R` documents how the analysis-ready file is constructed from authorized MIDUS source datasets. Those source datasets remain restricted and are not redistributed in this repository.

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

The Torch analysis environment used R 4.5.1. Package requirements are listed in `hpc/install_required_packages.R`, and the captured R session and package versions are stored in `sessionInfo.txt`.

## Final archival checklist

Before creating the manuscript release:

- upload the final primary simulation summaries;
- upload manuscript-relevant sensitivity summaries/QC files;
- upload disclosure-safe real-world aggregate outputs;
- upload final figures and machine-readable tables;
- verify that `sessionInfo.txt` reflects the archived analysis environment;
- confirm that no participant-level data or restricted metadata are tracked;
- update `CITATION.cff` with the final manuscript author list/version/DOI as appropriate;
- create a versioned GitHub release and archive it with a DOI service such as Zenodo.
