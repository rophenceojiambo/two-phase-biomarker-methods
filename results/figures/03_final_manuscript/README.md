# Final simulation figure selection

This directory contains the parsimonious figure set selected from the complete primary-simulation figure library for manuscript archiving.

The figures were generated from the validated primary simulation summaries by `09_make_simulation_figures.R`. The original PDF manuscript-candidate figures supplied from the completed simulation run are archived here directly.

## Core MIDUS-like figures

The MIDUS-like auxiliary-information regime uses \(R_A^2 = 0.05\) and \(R_Y^2 = 0.05\).

- `core/bias_MIDUS_like.pdf` — empirical bias across Phase-1 sample sizes, Phase-2 fractions, and non-null effect sizes.
- `core/empirical_se_MIDUS_like.pdf` — sampling variability of the estimators.
- `core/se_ratio_MIDUS_like.pdf` — model-based SE divided by empirical SE, directly assessing variance-estimator calibration.
- `core/coverage_MIDUS_like.pdf` — 95% confidence-interval coverage.
- `core/type1_MIDUS_like.pdf` — null rejection rates.

These five figures form the primary simulation-performance set because together they address bias, precision, standard-error calibration, interval calibration, and Type I error.

## Strong auxiliary-information stress figures

The strong auxiliary-information regime uses \(R_A^2 = 0.10\) and \(R_Y^2 = 0.10\).

- `stress/bias_strong_auxiliary_information.pdf`
- `stress/coverage_strong_auxiliary_information.pdf`
- `stress/type1_strong_auxiliary_information.pdf`

These were retained to show how performance changes when the omitted Phase-2 marker information is most consequential.

## Supporting figures

- `supporting/power_MIDUS_like.pdf` — power under the MIDUS-like regime; interpret together with the corresponding Type I error figure.
- `supporting/runtime_by_method_and_N.pdf` — computational cost by method and Phase-1 sample size.

## Candidate figures not duplicated here

The complete candidate and full-factorial figure libraries can be regenerated from `09_make_simulation_figures.R` and the archived machine-readable summaries.

Mean-estimate plots are redundant with bias; model-SE and relative-error plots are summarized more directly by the SE-ratio display; MSE and relative-precision plots are secondary summaries; bias-eliminated coverage is a diagnostic rather than a primary manuscript display; and the weak auxiliary-information regime adds little beyond the full factorial summaries for the final manuscript-facing archive.

No participant-level MIDUS data are contained in these figures.
