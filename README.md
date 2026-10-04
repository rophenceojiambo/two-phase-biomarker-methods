# Incorporating Partially Observed RNA-Based Leukocyte Proxies for Cellular Heterogeneity in Epigenetic Aging Studies: A Two-Phase Simulation Study

This repository contains the R code and results for a simulation study comparing methods for handling partially observed leukocyte-marker covariates in a two-phase biomarker design, along with an application to MIDUS data.

Participant-level MIDUS data are not included.

## Methods compared

The simulation compares six approaches:

1. Naive Phase-1 analysis that omits the Phase-2 RNA markers
2. Complete-case analysis restricted to the Phase-2 sample
3. Fully conditional specification multiple imputation (FCS-MI) using `mice`
4. Joint-model multiple imputation (JM-MI) using `jomo`
5. Inverse probability weighting (IPW) with stacked sandwich variance
6. Augmented inverse probability weighting (AIPW) with stacked sandwich variance

The MIDUS application uses the same six methods.

## Primary simulation design

The primary simulation includes 324 scenarios defined by:

- Phase-1 sample size: 500, 800, or 1,500
- Phase-2 sampling fraction: 0.25, 0.50, or 0.65
- Marker contribution to the exposure: R² = 0.01, 0.05, or 0.10
- Marker contribution to the outcome: R² = 0.00, 0.02, 0.05, or 0.10
- True exposure effect: 0.00, 0.15, or 0.30

Each scenario uses 2,000 Monte Carlo repetitions. Multiple-imputation analyses use 20 imputations. The production JM-MI setup uses 1,000 burn-in iterations and 1,000 iterations between imputations.

## Repository layout

The numbered R scripts stay at the repository root because they source one another by filename and were run from the project root.

- `00_config.R` — shared settings, paths, simulation grid, random-number settings, and reporting thresholds
- `prepare_midus_analytic_sample.R` — builds the analysis-ready MIDUS dataset from authorized source files
- `01_calibrate_midus.R` — calibration using the MIDUS analytic data
- `02_dgm.R` — primary data-generating mechanism
- `03_methods.R` — the six analysis methods
- `04_build_dgm_cache.R` — primary DGM cache and checks
- `05_*.R`–`09_*.R` — primary simulation, summaries, tables, and figures
- `10_*.R`–`15_*.R` — AIPW robustness analysis
- `16_*.R`–`21_*.R` — stronger AIPW misspecification analysis
- `22_*.R`–`27_*.R` — empirical-residual sensitivity analysis
- `28_*.R`–`33_*.R` — Phase-2 MCAR sensitivity analysis
- `34_*.R`–`38_*.R` — FCS predictive mean matching sensitivity analysis
- `39_*.R`–`43_*.R` — IPW weight-capping sensitivity analysis
- `44_*.R`–`48_*.R` — MIDUS application and tables/figures
- `hpc/` — SLURM scripts used on NYU Torch
- `data/` — notes on the MIDUS files needed to run the application
- `results/` — simulation, sensitivity, and MIDUS results
- `docs/` — reproducibility notes
- `validation/` — checks of the IPW/AIPW sandwich variance calculations completed before the final production simulation

## MIDUS application

The application uses standardized everyday discrimination as the exposure and standardized GrimAge2 and DunedinPACE as the outcomes. Phase-1 covariates are standardized age, sex, and race/ethnicity. The Phase-2 covariates are eight standardized RNA leukocyte-marker transcripts. BMI and assay plate variables are not included in the final analysis.

If you have access to the MIDUS data, place the analysis-ready file at:

```
data/MIDUS_discrimination_analysis.rds
```

or set its location with the `SIM_ANALYTIC_RDS` environment variable. The data file is not included in this repository.

## Data

Participant-level MIDUS data are not included. Researchers who want to reproduce the MIDUS application need to obtain the relevant MIDUS files through the appropriate MIDUS data-access process.

The script `prepare_midus_analytic_sample.R` shows how the analysis-ready file was created from the MIDUS source files. The source data and the resulting participant-level analysis file are not included.

## Code

The repository includes the simulation code, sensitivity analyses, MIDUS analysis, tables and figures, and the HPC scripts used for the study.

The full simulation was run on NYU Torch using the SLURM scripts in `hpc/`. The summarized results used in the paper are included, so readers do not need to rerun all Monte Carlo repetitions to inspect the reported results.

## Results

For the files used directly in the manuscript, see:

```
results/manuscript/
```

This folder contains Figures 1–6, Supplementary Figures S1–S10, Tables 3–4 data, and Supplementary Tables S1–S6. The detailed simulation, sensitivity, QC, and MIDUS result files remain in their analysis-specific folders under `results/`.

Raw Monte Carlo chunk files, checkpoints, caches, and participant-level MIDUS data are not included.

## Software

The analyses were run in R. The Torch workflow used R 4.5.1. Required packages are listed in the analysis scripts and `hpc/install_required_packages.R`.

The R session used for the analysis is recorded in `sessionInfo.txt`.

## Citation

Use `CITATION.cff` or GitHub's **Cite this repository** option. Once a DOI is created for the release, it should also be added to the README and `CITATION.cff`.

## License

The code is released under the MIT License. The license does not apply to MIDUS data or other third-party materials.

## Current status

The analysis code, summarized results, manuscript figures and tables, MIDUS preparation code, software information, and variance-validation checks are included. Before creating the final release, the remaining steps are to finish the last repository review, update the result manifest, update the citation information, and create the release/DOI.
