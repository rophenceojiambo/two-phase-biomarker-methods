# Data

## MIDUS empirical application

Participant-level MIDUS data are **not distributed in this public repository**.

The empirical application code is intended to be run only by users with appropriate authorization to access the required MIDUS files. Restricted source data, derived participant-level analytic datasets, and participant-level intermediate objects must remain outside the repository.

The repository may contain:

- variable specifications;
- data-preparation code;
- analysis code;
- synthetic or simulated example data;
- aggregate, disclosure-safe tables and figures used in the manuscript.

The repository must not contain:

- raw MIDUS participant-level files;
- analysis-ready MIDUS participant-level datasets;
- participant identifiers;
- restricted metadata;
- participant-level intermediate exports.

The root `.gitignore` contains additional safeguards against committing common restricted-data file formats.

## Suggested local layout

Authorized users may keep restricted files in a local directory such as:

```
data/
├── raw/          # ignored by Git
├── processed/    # ignored by Git
└── README.md
```

Analysis scripts should use project-relative paths or environment variables rather than hard-coded personal paths.
