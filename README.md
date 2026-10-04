# Two-Phase Biomarker Methods

This repository contains the analysis code supporting a methodological study of approaches for handling partially observed leukocyte-marker covariates in a two-phase biomarker design, together with a real-world application using MIDUS data.

The repository is intended to support the manuscript's **code availability** and **data availability** statements. Participant-level MIDUS data are not redistributed here.

## Methods evaluated

The simulation compares six analysis strategies:

1. Naive phase-1 analysis that omits the phase-2 RNA markers
2. Complete-case analysis restricted to the phase-2 sample
3. Fully conditional specification multiple imputation (FCS-MI) using `mice`
4. Joint-model multiple imputation (JM-MI) using `jomo`
5. Inverse probability weighting (IPW) with stacked sandwich variance
6. Augmented inverse probability weighting (AIPW) with stacked sandwich variance

The real-world MIDUS application uses the same six estimators.

## Primary simulation design

The primary factorial simulation contains 324 scenarios defined by:

- Phase-1 sample size: 500, 800, or 1,500
- Phase-2 sampling fraction: 0.25, 0.50, or 0.65
- Marker contribution to the exposure: R² = 0.01, 0.05, or 0.10
- Marker contribution to the outcome: R² = 0.00, 0.02, 0.05, or 0.10
- True exposure effect: 0.00, 0.15, or 0.30

Each primary scenario uses 2,000 Monte Carlo repetitions. Multiple-imputation analyses use 20 imputations; the production JM-MI configuration uses 1,000 burn-in iterations and 1,000 iterations between imputations.

## Repository structure

The validated analysis scripts are intentionally retained at the repository root because they source one another by filename and were executed from the project root.

- `00_config.R` — shared configuration, paths, simulation grid, random-number settings, and reporting thresholds
- `prepare_midus_analytic_sample.R` — constructs the analysis-ready MIDUS dataset from authorized source files
- `01_calibrate_midus.R` — calibration from the authorized MIDUS analytic data
- `02_dgm.R` — primary data-generating mechanism
- `03_methods.R` — implementations of the six analysis methods
- `04_build_dgm_cache.R` — primary DGM cache and checks
- `05_*.R`–`09_*.R` — primary production, combination, summaries, manuscript tables, and figures
- `10_*.R`–`15_*.R` — AIPW robustness sensitivity analysis
- `16_*.R`–`21_*.R` — stronger AIPW misspecification sensitivity analysis
- `22_*.R`–`27_*.R` — empirical-residual sensitivity analysis
- `28_*.R`–`33_*.R` — phase-2 MCAR sensitivity analysis
- `34_*.R`–`38_*.R` — FCS predictive mean matching sensitivity analysis
- `39_*.R`–`43_*.R` — IPW weight-capping sensitivity analysis
- `44_*.R`–`48_*.R` — MIDUS real-world application and manuscript-output generation
- `hpc/` — SLURM/Torch submission, checking, and orchestration scripts
- `data/` — data-access documentation only; restricted MIDUS data are excluded
- `results/` — disclosure-safe manuscript-facing results to be archived with the paper
- `docs/` — reproducibility notes and repository guidance
- `validation/` — pre-production implementation validation of the IPW/AIPW stacked-sandwich variance estimators

## MIDUS empirical application

The empirical application evaluates standardized everyday discrimination as the exposure and standardized GrimAge2 and DunedinPACE as outcomes. Phase-1 covariates are standardized age, sex, and race/ethnicity; the phase-2 covariates are eight standardized RNA leukocyte-marker transcripts. BMI and assay plate variables are not included in the final analysis.

Authorized users should place the analysis-ready file at:

```
data/MIDUS_discrimination_analysis.rds
```

or provide its location through the `SIM_ANALYTIC_RDS` environment variable. This file is intentionally excluded from Git.

## Data availability

Participant-level MIDUS data are not distributed in this repository. Researchers wishing to reproduce the empirical application must obtain the relevant MIDUS data under the applicable MIDUS access and use conditions.

The public repository may contain only disclosure-safe derived summaries, manuscript tables, figures, and code. Raw data, analysis-ready participant-level data, identifiers, participant-level intermediate files, and restricted metadata must not be committed.

The shareable preparation script `prepare_midus_analytic_sample.R` documents how the analysis-ready MIDUS file is constructed from the authorized source datasets. The source datasets themselves remain restricted and are not included.

## Code availability

All simulation, sensitivity-analysis, real-world analysis, manuscript-output, and HPC orchestration code used for the study is maintained in this repository. The numbered files preserve the execution structure used for the completed analyses.

The full simulation was run on NYU Torch using the SLURM workflow in `hpc/`. Re-running the full simulation is not required to use this repository as the manuscript code archive; the archived summaries and manuscript-facing outputs allow readers to inspect the reported findings without reproducing all Monte Carlo runs.

## Results archive

The public results archive contains summarized, disclosure-safe primary simulation results, sensitivity-analysis summaries and QC files, calibration summaries, and aggregate real-world MIDUS outputs. Raw Monte Carlo repetitions, caches, checkpoints, and participant-level data are intentionally excluded.

The final manuscript-facing outputs are organized under `results/manuscript/`. Main-text figures are numbered exactly as Figures 1–6 in the manuscript, and supplementary simulation figures are numbered S1–S10. The underlying machine-readable simulation, sensitivity, QC, and real-world outputs remain in their analysis-specific folders under `results/`. See `results/README.md` for the archive layout.

## Software

The analyses are implemented in R. The Torch workflow uses R 4.5.1 by default. Required packages are documented in the analysis scripts and `hpc/install_required_packages.R`.

The exact software environment captured from the analysis setup is recorded in `sessionInfo.txt`.

## Citation

Use the repository's `CITATION.cff` file or GitHub's **Cite this repository** function. Once a versioned archival DOI is available, the DOI should be added to both this README and `CITATION.cff`.

## License

The source code is released under the MIT License. This license applies to the repository code, not to MIDUS data or other third-party materials.

## Repository status

The cleaned analysis and HPC code, shareable MIDUS preparation script, software-environment record, disclosure-safe primary simulation summaries, sensitivity-analysis results, calibration summaries, aggregate real-world manuscript outputs, and pre-production IPW/AIPW variance-validation evidence are archived. Before the manuscript repository is frozen as a versioned release, the remaining tasks are to complete the final disclosure and integrity review, update citation metadata as appropriate, regenerate the public-results manifest for the frozen archive, and create the tagged archival release.
