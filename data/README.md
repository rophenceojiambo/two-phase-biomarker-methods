# Data

Participant-level MIDUS data are **not distributed in this public repository**.

## Required empirical-analysis input

Scripts `44`–`48` expect the authorized analysis-ready file:

```
data/MIDUS_discrimination_analysis.rds
```

The location may also be supplied through the `SIM_ANALYTIC_RDS` environment variable where supported.

Restricted source data, the analysis-ready dataset, participant identifiers, participant-level intermediate objects, and restricted metadata must remain outside GitHub.

The repository may contain:

- shareable data-preparation and variable-construction code;
- simulation scenario definitions;
- synthetic or simulated example data where appropriate;
- aggregate disclosure-safe tables, figures, and QC summaries used in the manuscript.

## Reproducibility gap to close before final release

The current repository contains the analysis code beginning from `MIDUS_discrimination_analysis.rds`, but it does not yet contain the upstream script that constructs that file from the authorized MIDUS source datasets. That script should be added to the public repository if allowed by the applicable data-use terms.

## Local data layout

An authorized user can keep the restricted input locally as:

```
data/
├── MIDUS_discrimination_analysis.rds   # ignored by Git
└── README.md
```

The root `.gitignore` prevents common restricted-data formats from being committed.
