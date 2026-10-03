# Results archive

This directory contains disclosure-safe outputs from the completed analyses.

## Reader-facing manuscript outputs

Use `results/manuscript/` for the clean manuscript-facing archive.

- `results/manuscript/main/figures/` contains Figures 1–6 exactly as referenced in the manuscript.
- `results/manuscript/main/tables/` contains the data tables used for the MIDUS Results section (Tables 3–4). Tables 1–2 are simulation-design/methods tables embedded directly in the manuscript rather than generated empirical-result tables.
- `results/manuscript/supplementary/figures/` contains Supplementary Figures S1–S10 exactly as referenced in the manuscript.
- `results/manuscript/supplementary/tables/` contains Supplementary Tables S1–S6 derived from the canonical sensitivity-analysis summaries.

Power is a main-text performance result (Figure 5) and therefore belongs with the other main manuscript figures rather than in a separate supporting folder.

## Machine-readable primary simulation outputs

The primary simulation summaries remain under `results/summary/`, with manuscript-oriented summary tables under `results/tables/`, QC under `results/qc/`, and figure metadata under `results/figure_data/`.

These include the primary performance, Type I error, failure, runtime, Phase-2, and weight-diagnostic summaries needed to reproduce the manuscript displays.

## Sensitivity and robustness outputs

The full numerical outputs supporting Supplementary Tables S1–S6 remain in their analysis-specific directories:

- AIPW robustness: `results/robustness/aipw_double_robustness/`
- Strong AIPW misspecification: `results/sensitivity/aipw_strong_misspecification/`
- Empirical-residual sensitivity: `results/empirical_residual_sensitivity/`
- Phase-2 MCAR sensitivity: `results/phase2_mcar_sensitivity/`
- FCS predictive mean matching sensitivity: `results/fcs_pmm_sensitivity/`
- IPW weight-capping sensitivity: `results/ipw_weight_capping_sensitivity/`

Raw chunk files, checkpoints, scenario-level RDS files, and rsimsum object files are intentionally excluded.

## MIDUS real-world application

Aggregate MIDUS outputs remain under `results/real_world_application/`. The annotated method-comparison figure used as Figure 6 and the data tables used for Tables 3–4 are also exposed under `results/manuscript/main/` so that readers can find all manuscript-facing outputs in one place. These reader-facing paths intentionally point to the same underlying Git blobs when the files are identical; they do not represent separate analytical results.

Participant-level MIDUS data are not included.

## Do not archive

Do not commit raw or analysis-ready participant-level MIDUS data, identifiers, participant-level intermediate exports, raw Monte Carlo chunk files, checkpoints, cached simulated datasets, or SLURM scratch/log files.

The repository `.gitignore` blocks the most common restricted and regenerable file types.
