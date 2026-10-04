# Pre-production variance validation

This directory archives implementation-validation work completed before the final production simulation workflow was frozen.

It is intentionally separate from the numbered `00`–`48` production scripts and from `results/`. These files are **validation evidence**, not the source of the final primary simulation figures or headline numerical results.

## Purpose

The custom stacked-sandwich variance estimators used for IPW and AIPW were checked because model-based standard errors were smaller than empirical Monte Carlo standard errors in the sparse Phase-2 stress setting.

The archived diagnostic compared, on the same simulated datasets:

1. the hand-derived analytic stacked-sandwich standard error;
2. a numerical stacked-sandwich standard error using the same estimating equations with a central finite-difference Jacobian; and
3. a participant-level nonparametric bootstrap standard error.

The diagnostic focused on the stress setting:

- Phase-1 sample size: N = 500
- Phase-2 fraction: 0.25
- exposure-marker R-squared: 0.10
- marker-outcome R-squared: 0.10
- true exposure coefficient: 0.15
- 20 independently simulated verification datasets
- 250 bootstrap resamples per verification dataset

## Key finding

The analytic and numerical stacked-sandwich standard errors agreed essentially exactly for both IPW and AIPW. The maximum absolute difference between analytic and numerical point estimates was zero to numerical precision, and estimating-equation root checks were negligible.

In the stress benchmark, the earlier 2,000-repetition validation showed:

- IPW empirical SE = 0.09329, mean model SE = 0.08164, SE ratio = 0.875, coverage = 0.911;
- AIPW empirical SE = 0.06121, mean model SE = 0.05802, SE ratio = 0.948, coverage = 0.9335.

Across the 20 bootstrap-verification datasets:

- IPW mean analytic/numerical SE = 0.08049 and mean bootstrap SE = 0.09190;
- AIPW mean analytic/numerical SE = 0.05700 and mean bootstrap SE = 0.07034.

These diagnostics support the interpretation that the hand-derived Jacobian was implemented consistently with the stacked estimating equations. The remaining variance underestimation, particularly for IPW under sparse Phase-2 sampling, was treated as a finite-sample issue rather than a discrepancy between the analytic and numerical Jacobians.

## Archived files

### Script

`scripts/10_verify_sandwich_bootstrap_preproduction.R`

Historical diagnostic script as executed before the production simulation was frozen. The original internal script number is retained in the filename because it predates the final `10`–`48` numbering used by the production workflow.

### Compact results

- `results/sandwich_bootstrap_verification_results.csv` — 40 rows: IPW and AIPW results for 20 verification datasets.
- `results/sandwich_bootstrap_verification_summary.csv` — compact analytic/numerical/bootstrap comparison.
- `results/variance_validation_performance_summary.csv` — 2,000-repetition pre-production benchmark for CCA, IPW, and AIPW in the MIDUS-like and stress settings.
- `results/variance_validation_phase2_summary.csv` — realized Phase-2 sampling diagnostics for the benchmark.
- `results/variance_validation_weight_summary.csv` — IPW/AIPW weight diagnostics for the benchmark.
- `results/variance_validation_failure_details.csv` — failure-detail file; header-only because no failures were recorded.

## Files intentionally not archived

The following uploaded working files are intentionally excluded from GitHub:

- checkpoint RDS files;
- RDS copies of the validation results;
- the 12,000-row `variance_validation_all_results.csv`;
- the 6,000-row scenario-specific MIDUS-like and stress result files;
- scenario-specific RDS files.

Those files are pre-production repetition-level/checkpoint objects and are unnecessary for evaluating the implementation-validation conclusion. The compact summaries above retain the information needed for the manuscript supplement.

## Relation to final production results

The final primary simulation in `results/summary/` supersedes the earlier 2,000-repetition variance-validation benchmark for manuscript performance estimates. Therefore, values in this directory should not be substituted for the final production estimates, figures, or tables. They are retained solely to document implementation validation of the IPW/AIPW variance estimators.
