# Results

This folder contains the summarized results from the simulation study, sensitivity analyses, and MIDUS application.

## Files used in the manuscript

The easiest place to find the figures and tables used in the paper is:

```
results/manuscript/
```

It contains:

- `main/figures/` — Figures 1–6
- `main/tables/` — data used for Tables 3–4
- `supplementary/figures/` — Supplementary Figures S1–S10
- `supplementary/tables/` — Supplementary Tables S1–S6

Power is Figure 5 in the main text, so it is stored with the other main simulation figures.

## Primary simulation results

The detailed primary-simulation summaries are kept in:

- `results/summary/`
- `results/tables/`
- `results/qc/`
- `results/figure_data/`

These files contain the performance summaries, Type I error, failures, runtime, Phase-2 summaries, weight checks, and the data used to make the figures.

## Sensitivity analyses

The detailed sensitivity results are kept in their own folders:

- AIPW robustness: `results/robustness/aipw_double_robustness/`
- Stronger AIPW misspecification: `results/sensitivity/aipw_strong_misspecification/`
- Empirical residuals: `results/empirical_residual_sensitivity/`
- Phase-2 MCAR: `results/phase2_mcar_sensitivity/`
- FCS predictive mean matching: `results/fcs_pmm_sensitivity/`
- IPW weight capping: `results/ipw_weight_capping_sensitivity/`

Large chunk files, checkpoints, and scenario-level RDS files are not included.

## MIDUS application

Aggregate MIDUS results are in `results/real_world_application/`.

The annotated figure used as Figure 6 and the data used for Tables 3–4 are also placed in `results/manuscript/main/` so the manuscript files are easy to find. In cases where the same file appears in both places, it is the same result saved under the manuscript organization as well.

Participant-level MIDUS data are not included.

## Files that should stay out of GitHub

Do not add:

- raw or analysis-ready participant-level MIDUS data;
- participant identifiers;
- participant-level intermediate exports;
- raw Monte Carlo chunk files;
- checkpoints and simulation caches;
- SLURM scratch files and logs.

The root `.gitignore` blocks the common file types used for these files.
