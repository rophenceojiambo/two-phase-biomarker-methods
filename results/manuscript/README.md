# Manuscript-facing results

This directory is the clean reader-facing layer for outputs cited directly in the manuscript.

## Main text

### Figures

- `main/figures/Figure1_bias_MIDUS_like.pdf`
- `main/figures/Figure2_se_ratio_MIDUS_like.pdf`
- `main/figures/Figure3_coverage_MIDUS_like.pdf`
- `main/figures/Figure4_type1_MIDUS_like.pdf`
- `main/figures/Figure5_power_MIDUS_like.pdf`
- `main/figures/Figure6_real_world_method_comparison.pdf`

Figures 1–5 use the MIDUS-like auxiliary-information regime. Figure 6 is the annotated MIDUS Refresher method-comparison figure used in the manuscript.

### Tables

- `main/tables/Table3_characteristics_by_phase2.csv`
- `main/tables/Table4_phase2_correlations_display.csv`

Tables 1–2 in the manuscript describe the simulation design and performance-measure definitions and are embedded directly in the manuscript rather than treated as empirical-result outputs.

## Supplementary figures

The supplementary figure directory mirrors the numbering used in Additional File 1:

- S1–S2: empirical bias under weak and strong auxiliary information
- S3–S4: standard-error ratio under weak and strong auxiliary information
- S5–S6: 95% confidence-interval coverage under weak and strong auxiliary information
- S7–S8: null rejection rate / Type I error under weak and strong auxiliary information
- S9–S10: power under weak and strong auxiliary information

The weak regime uses \(R_A^2=0.01\), \(R_Y^2=0.02\); the strong regime uses \(R_A^2=0.10\), \(R_Y^2=0.10\).

## Supplementary tables

Supplementary Tables S1–S6 are archived under `supplementary/tables/`:

- `TableS1_AIPW_standard_misspecification.csv`
- `TableS2_AIPW_strong_misspecification.csv`
- `TableS3_empirical_residual_sensitivity.csv`
- `TableS4_phase2_MCAR_sensitivity.csv`
- `TableS5_FCS_PMM_sensitivity.csv`
- `TableS6_IPW_weight_capping_sensitivity.csv`

These are compact manuscript-facing tables derived from the canonical robustness and sensitivity summary files retained elsewhere under `results/`.

## Underlying numerical results

This manuscript-facing directory intentionally avoids duplicating the full computational archive. Primary simulation summaries, sensitivity-analysis results, QC files, and aggregate real-world outputs remain in their analysis-specific directories under `results/`.
