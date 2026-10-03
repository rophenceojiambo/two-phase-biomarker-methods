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

The empirical application requires the analysis-ready MIDUS file. By default the workflow looks for:

```
data/MIDUS_discrimination_analysis.rds
```

Authorized users may instead keep the restricted file outside the repository and provide its path through the `SIM_ANALYTIC_RDS` environment variable.

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

The primary simulation summaries and manuscript-facing tables are archived under `results/summary/` and `results/tables/`. Sensitivity-analysis summaries and QC outputs are archived under their corresponding subdirectories in `results/`. The complete generated simulation-figure library is represented by machine-readable figure data and a figure manifest; only the exact figures selected for the manuscript should be copied into `results/figures/03_final_manuscript/` for the frozen release.

## Real-world outputs

Scripts `44`–`48` generate aggregate method estimates, sample summaries, weight diagnostics, manuscript tables, figures, and QC files. Disclosure-safe aggregate outputs are archived under `results/real_world_application/`; participant-level MIDUS data are not included.

## Software environment

The Torch analysis environment used R 4.5.1. Package requirements are listed in `hpc/install_required_packages.R`, and the captured R session and package versions are stored in `sessionInfo.txt`.

## Archival status

Completed:

- primary simulation summaries and manuscript-facing tables are archived;
- manuscript-relevant sensitivity summaries and QC files are archived;
- disclosure-safe aggregate real-world outputs, tables, figures, and QC files are archived;
- calibration summaries used by the simulation are archived;
- `sessionInfo.txt` records the analysis software environment;
- participant-level MIDUS data, raw Monte Carlo chunks, checkpoints, and caches are excluded from the public repository.

Remaining before creating the manuscript release:

- select and archive the exact primary simulation figures used in the manuscript under `results/figures/03_final_manuscript/`;
- complete the final disclosure and repository-integrity review, including confirmation that public calibration and aggregate MIDUS outputs comply with applicable MIDUS requirements;
- regenerate `results/PUBLIC_RESULTS_MANIFEST.csv` after the archive contents are final so its checksums describe the frozen release;
- update `CITATION.cff` with the final manuscript author list, version, and archival DOI as appropriate;
- create a versioned GitHub release and archive it with a DOI service such as Zenodo.
