# Variance validation

This folder contains checks I ran before the final production simulation to make sure the custom IPW and AIPW sandwich variance calculations were working as intended.

These checks are separate from the final simulation results.

## What I checked

For the same simulated datasets, I compared:

1. the analytic stacked-sandwich standard error used in the main methods;
2. the same stacked estimating equations with a numerical finite-difference Jacobian; and
3. a nonparametric bootstrap that resampled Phase-1 participants and refit the models.

The check focused on the stress setting:

- N = 500
- Phase-2 fraction = 0.25
- exposure-marker R² = 0.10
- marker-outcome R² = 0.10
- true exposure coefficient = 0.15
- 20 simulated datasets
- 250 bootstrap resamples per dataset

## What the check showed

The analytic and numerical sandwich standard errors were essentially identical for both IPW and AIPW. The point estimates also matched to numerical precision, and the estimating-equation root checks were very small.

In the earlier 2,000-repetition stress check:

- IPW: empirical SE = 0.09329, mean model SE = 0.08164, SE ratio = 0.875, coverage = 0.911
- AIPW: empirical SE = 0.06121, mean model SE = 0.05802, SE ratio = 0.948, coverage = 0.9335

Across the 20 bootstrap datasets:

- IPW: mean analytic/numerical SE = 0.08049; mean bootstrap SE = 0.09190
- AIPW: mean analytic/numerical SE = 0.05700; mean bootstrap SE = 0.07034

The main conclusion from this check was that the analytic Jacobian agreed with the numerical Jacobian. The remaining SE underestimation, especially for IPW when the Phase-2 sample was small, was therefore treated as a finite-sample issue rather than a coding mismatch between the two Jacobians.

## Files

### Script

`scripts/verify_sandwich_bootstrap.R`

This is the diagnostic script that was used before the final production simulation scripts were set.

### Results

- `results/sandwich_bootstrap_verification_results.csv` — IPW and AIPW results for the 20 verification datasets
- `results/sandwich_bootstrap_verification_summary.csv` — summary of the analytic, numerical, and bootstrap SE comparison
- `results/variance_validation_performance_summary.csv` — earlier 2,000-repetition check for CCA, IPW, and AIPW in the MIDUS-like and stress settings
- `results/variance_validation_phase2_summary.csv` — Phase-2 sample-size checks
- `results/variance_validation_weight_summary.csv` — IPW/AIPW weight summaries
- `results/variance_validation_failure_details.csv` — failure file; it contains only the header because there were no failures

## Files not included

I did not include the checkpoint files, RDS copies, or the large repetition-level validation files. The summary files above contain what is needed to understand the validation check.

## Relation to the final simulation

The final simulation results in `results/summary/` are the results used for the paper. The earlier 2,000-repetition files in this folder are kept only to document the variance check and should not be used in place of the final simulation results.
