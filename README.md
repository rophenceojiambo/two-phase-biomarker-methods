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

**Important reproducibility note:** scripts `44`–`48` begin from the analysis-ready file `MIDUS_discrimination_analysis.rds`. The script that constructs that file from the authorized MIDUS source datasets is not currently included in this repository. That preparation script should be added before the final archival release if it can be shared under the relevant data-use terms.

## Code availability

All simulation, sensitivity-analysis, real-world analysis, manuscript-output, and HPC orchestration code used for the study is maintained in this repository. The numbered files preserve the execution structure used for the completed analyses.

The full simulation was run on NYU Torch using the SLURM workflow in `hpc/`. Re-running the full simulation is not required to use this repository as the manuscript code archive; the final release will also contain the compact manuscript-facing results needed to inspect the reported findings.

## Results archive

The public results archive should contain summarized, disclosure-safe outputs rather than repetition-level simulation files. See `results/README.md` for the files recommended for archival.

## Software

The analyses are implemented in R. The Torch workflow uses R 4.5.1 by default. Required packages are documented in the analysis scripts and `hpc/install_required_packages.R`.

For the final archived release, an exact environment record such as `sessionInfo.txt` or an `renv.lock` file is recommended if available from the analysis environment.

## Citation

Use the repository's `CITATION.cff` file or GitHub's **Cite this repository** function. Once a versioned archival DOI is available, the DOI should be added to both this README and `CITATION.cff`.

## License

The source code is released under the MIT License. This license applies to the repository code, not to MIDUS data or other third-party materials.

## Repository status

The cleaned analysis and HPC code are archived. Before the manuscript repository is frozen as a versioned release, the remaining tasks are to add the disclosure-safe final results, add the MIDUS data-preparation script if distributable, record the final software environment if available, and create a tagged release.
