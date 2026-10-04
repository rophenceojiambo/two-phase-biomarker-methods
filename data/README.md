# Data

Participant-level MIDUS data are not included in this repository.

## File needed for the MIDUS analysis

Scripts `44`–`48` use the analysis-ready file:

```
data/MIDUS_discrimination_analysis.rds
```

The file's location can be supplied with the `SIM_ANALYTIC_RDS` environment variable.

This repository includes:

- the code used to prepare the analysis dataset;
- simulation settings;
- summarized simulation results;
- aggregate MIDUS tables, figures, and QC results used for the paper.

## Preparing the MIDUS analysis file

`prepare_midus_analytic_sample.R` shows how `MIDUS_discrimination_analysis.rds` was created from the MIDUS source files.

The source data and the resulting participant-level analysis file are not included.

## Local setup

With access to the MIDUS data, the local folder can look like this:

```
data/
├── MIDUS_discrimination_analysis.rds   # ignored by Git
└── README.md
```

The root `.gitignore` blocks the common data-file formats used in this project.
