# Empirical application code

This directory contains code for the real-world MIDUS application.

Participant-level MIDUS data are not stored in this repository. Scripts added here should assume that authorized users have placed the required source data in ignored local directories.

Planned workflow:

1. Prepare the analytic dataset from authorized MIDUS source files.
2. Generate descriptive summaries.
3. Fit the primary analysis.
4. Run prespecified sensitivity analyses.
5. Generate manuscript tables.
6. Generate manuscript figures.

All scripts should avoid hard-coded user-specific paths and should write only disclosure-safe manuscript outputs to tracked directories.
