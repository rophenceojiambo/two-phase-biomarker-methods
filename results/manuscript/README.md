# Results used in the manuscript

This folder collects the figures and tables that are used directly in the paper.

## Main text

### Figures

- `main/figures/Figure1_bias_MIDUS_like.pdf`
- `main/figures/Figure2_se_ratio_MIDUS_like.pdf`
- `main/figures/Figure3_coverage_MIDUS_like.pdf`
- `main/figures/Figure4_type1_MIDUS_like.pdf`
- `main/figures/Figure5_power_MIDUS_like.pdf`
- `main/figures/Figure6_real_world_method_comparison.pdf`

Figures 1–5 use the MIDUS-like simulation setting. Figure 6 is the annotated MIDUS method-comparison figure used in the paper.

### Tables

- `main/tables/Table3_characteristics_by_phase2.csv`
- `main/tables/Table4_phase2_correlations_display.csv`

Tables 1–2 describe the simulation design and performance measures and are written directly in the manuscript.

## Supplementary figures

The supplementary figures use the same numbering as the manuscript supplement:

- S1–S2: bias under weak and strong auxiliary information
- S3–S4: SE ratio under weak and strong auxiliary information
- S5–S6: 95% CI coverage under weak and strong auxiliary information
- S7–S8: Type I error under weak and strong auxiliary information
- S9–S10: power under weak and strong auxiliary information

The weak setting uses (R_A^2=0.01) and (R_Y^2=0.02). The strong setting uses (R_A^2=0.10) and (R_Y^2=0.10).

## Supplementary tables

Supplementary Tables S1–S6 are in `supplementary/tables/`:

- `TableS1_AIPW_standard_misspecification.csv`
- `TableS2_AIPW_strong_misspecification.csv`
- `TableS3_empirical_residual_sensitivity.csv`
- `TableS4_phase2_MCAR_sensitivity.csv`
- `TableS5_FCS_PMM_sensitivity.csv`
- `TableS6_IPW_weight_capping_sensitivity.csv`

The more detailed result files used to create these tables remain in the corresponding sensitivity-analysis folders under `results/`.
