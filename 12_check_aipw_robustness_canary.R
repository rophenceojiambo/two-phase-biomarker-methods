# 12_check_aipw_robustness_canary.R
#
# Validates the structure and implementation of the AIPW double-robustness
# canary run produced by 11_run_aipw_robustness_canary.R.
#
# The canary is an implementation check, not a final Monte Carlo experiment.
# With only 20 repetitions per generating setting, its performance summaries
# should not be interpreted as substantive evidence about double robustness.
#
# This script verifies that:
# - every generated dataset is analyzed under the four prespecified models;
# - the four specification labels occur exactly once per dataset;
# - the both-correct AIPW analysis reproduces the validated primary result;
# - omitting Y meaningfully changes the estimated selection probabilities;
# - omitting Y meaningfully changes the marker-model augmentation quantities;
# - no unexpected AIPW failures occur; and
# - diagnostic calculations are available for every generated dataset.

# Load robustness settings and helper functions
source("10_aipw_robustness_helpers.R")

library(dplyr)
library(readr)

# 1. Input and output paths

robustness_root <- file.path(
  results_dir,
  "robustness",
  "aipw_double_robustness"
)

canary_dir <- file.path(
  robustness_root,
  "canary"
)

results_file <- file.path(
  canary_dir,
  "aipw_robustness_canary_results.csv"
)

diagnostics_file <- file.path(
  canary_dir,
  "aipw_robustness_canary_misspecification_diagnostics.csv"
)

if (!file.exists(results_file)) {
  stop(
    "Canary results not found:\n",
    results_file,
    "\nRun 11_run_aipw_robustness_canary.R first."
  )
}

if (!file.exists(diagnostics_file)) {
  stop(
    "Canary diagnostics not found:\n",
    diagnostics_file
  )
}

canary_results <- read_csv(
  results_file,
  show_col_types = FALSE
)

misspecification_diagnostics <- read_csv(
  diagnostics_file,
  show_col_types = FALSE
)

# Confirm that the required columns are available
check_columns <- function(data, required, data_name) {
  missing <- setdiff(required, names(data))
  
  if (length(missing) > 0L) {
    stop(
      data_name,
      " is missing: ",
      paste(missing, collapse = ", ")
    )
  }
}

check_columns(
  canary_results,
  c(
    "setting",
    "repetition",
    "specification",
    "selection_correct",
    "augmentation_correct",
    "status",
    "estimate",
    "se",
    "conf_low",
    "conf_high",
    "theta_true",
    "elapsed_seconds",
    "estimate_difference_from_primary",
    "se_difference_from_primary"
  ),
  "Canary results"
)

check_columns(
  misspecification_diagnostics,
  c(
    "setting",
    "repetition",
    "mean_abs_pi_difference",
    "max_abs_pi_difference",
    "cor_pi_correct_misspecified",
    "mean_abs_marker_mean_difference",
    "max_abs_marker_mean_difference",
    "relative_frobenius_sigma_difference",
    "diagnostic_message"
  ),
  "Canary misspecification diagnostics"
)

# 2. Structural QC

expected_specs <- aipw_robustness_specs$specification
expected_settings <- aipw_robustness_settings$setting
expected_n_specs <- length(expected_specs)

# Every setting-repetition should have four result rows
key_counts <- canary_results %>%
  count(
    setting,
    repetition,
    name = "n_rows"
  )

bad_key_counts <- key_counts %>%
  filter(n_rows != expected_n_specs)

# Every setting-repetition should contain the same four specification labels
specification_check <- canary_results %>%
  group_by(
    setting,
    repetition
  ) %>%
  summarise(
    specifications = paste(
      sort(unique(specification)),
      collapse = " | "
    ),
    n_specs = n_distinct(specification),
    .groups = "drop"
  )

expected_spec_string <- paste(
  sort(expected_specs),
  collapse = " | "
)

bad_spec_sets <- specification_check %>%
  filter(
    n_specs != expected_n_specs |
      specifications != expected_spec_string
  )

# Each specification should occur once within a setting-repetition
duplicate_rows <- canary_results %>%
  count(
    setting,
    repetition,
    specification,
    name = "n"
  ) %>%
  filter(n != 1L)

# Confirm that the expected generating settings are present
settings_match <- setequal(
  unique(canary_results$setting),
  expected_settings
)

# There should be one diagnostic row per generated dataset
diagnostic_key_counts <- misspecification_diagnostics %>%
  count(
    setting,
    repetition,
    name = "n_rows"
  )

bad_diagnostic_key_counts <- diagnostic_key_counts %>%
  filter(n_rows != 1L)

result_keys <- canary_results %>%
  distinct(
    setting,
    repetition
  )

diagnostic_keys <- misspecification_diagnostics %>%
  distinct(
    setting,
    repetition
  )

diagnostic_keys_match <- setequal(
  interaction(
    result_keys$setting,
    result_keys$repetition,
    drop = TRUE
  ),
  interaction(
    diagnostic_keys$setting,
    diagnostic_keys$repetition,
    drop = TRUE
  )
)

# 3. Exact comparison with the primary AIPW results

# Restrict the comparison to rows with finite primary differences
primary_comparison_rows <- canary_results %>%
  filter(
    specification == "Both correct",
    is.finite(estimate_difference_from_primary),
    is.finite(se_difference_from_primary)
  )

primary_match <- primary_comparison_rows %>%
  summarise(
    n_compared = n(),
    
    max_abs_estimate_difference = if (n() > 0L) {
      max(abs(estimate_difference_from_primary))
    } else {
      NA_real_
    },
    
    max_abs_se_difference = if (n() > 0L) {
      max(abs(se_difference_from_primary))
    } else {
      NA_real_
    }
  )

# Exact reproduction should differ only by machine-level numerical precision
primary_exact_pass <- with(
  primary_match,
  isTRUE(
    n_compared > 0L &&
      is.finite(max_abs_estimate_difference) &&
      is.finite(max_abs_se_difference) &&
      max_abs_estimate_difference < 1e-10 &&
      max_abs_se_difference < 1e-10
  )
)

# 4. Method failure QC

failure_summary <- canary_results %>%
  group_by(
    setting,
    specification
  ) %>%
  summarise(
    n = n(),
    
    n_failed = sum(
      status != "ok" |
        !is.finite(estimate) |
        !is.finite(se)
    ),
    
    failure_rate = n_failed / n,
    .groups = "drop"
  )

no_aipw_failures <- all(
  failure_summary$n_failed == 0L
)

# 5. Misspecification diagnostics

misspecification_summary <- misspecification_diagnostics %>%
  group_by(setting) %>%
  summarise(
    repetitions = n(),
    
    n_diagnostic_failures = sum(
      !is.na(diagnostic_message) &
        nzchar(diagnostic_message)
    ),
    
    avg_mean_abs_pi_difference = mean(
      mean_abs_pi_difference,
      na.rm = TRUE
    ),
    
    median_mean_abs_pi_difference = median(
      mean_abs_pi_difference,
      na.rm = TRUE
    ),
    
    avg_max_abs_pi_difference = mean(
      max_abs_pi_difference,
      na.rm = TRUE
    ),
    
    avg_cor_pi_correct_misspecified = mean(
      cor_pi_correct_misspecified,
      na.rm = TRUE
    ),
    
    avg_mean_abs_marker_mean_difference = mean(
      mean_abs_marker_mean_difference,
      na.rm = TRUE
    ),
    
    median_mean_abs_marker_mean_difference = median(
      mean_abs_marker_mean_difference,
      na.rm = TRUE
    ),
    
    avg_max_abs_marker_mean_difference = mean(
      max_abs_marker_mean_difference,
      na.rm = TRUE
    ),
    
    avg_relative_frobenius_sigma_difference = mean(
      relative_frobenius_sigma_difference,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

# Omission of Y should change both nuisance-model components by more than
# numerical rounding alone
nondegenerate_selection <- all(
  is.finite(
    misspecification_summary$avg_mean_abs_pi_difference
  ) &
    misspecification_summary$avg_mean_abs_pi_difference > 1e-6
)

nondegenerate_augmentation <- all(
  is.finite(
    misspecification_summary$
      avg_mean_abs_marker_mean_difference
  ) &
    misspecification_summary$
    avg_mean_abs_marker_mean_difference > 1e-6
)

no_diagnostic_failures <- all(
  misspecification_summary$n_diagnostic_failures == 0L
)

# 6. Descriptive performance preview

# Use a consistent indicator for successful finite results
performance_data <- canary_results %>%
  mutate(
    valid_result =
      status == "ok" &
      is.finite(estimate) &
      is.finite(se) &
      is.finite(conf_low) &
      is.finite(conf_high)
  )

performance_preview <- performance_data %>%
  group_by(
    setting,
    specification,
    selection_correct,
    augmentation_correct
  ) %>%
  summarise(
    nsim = sum(valid_result),
    
    mean_estimate = mean(
      estimate[valid_result],
      na.rm = TRUE
    ),
    
    bias =
      mean_estimate -
      first(theta_true),
    
    empirical_se = sd(
      estimate[valid_result],
      na.rm = TRUE
    ),
    
    average_model_se = mean(
      se[valid_result],
      na.rm = TRUE
    ),
    
    se_ratio =
      average_model_se /
      empirical_se,
    
    mse = mean(
      (
        estimate[valid_result] -
          first(theta_true)
      )^2,
      na.rm = TRUE
    ),
    
    coverage = mean(
      conf_low[valid_result] <= first(theta_true) &
        conf_high[valid_result] >= first(theta_true),
      na.rm = TRUE
    ),
    
    mean_elapsed_seconds = mean(
      elapsed_seconds,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

# These summaries are descriptive only because the canary sample is small
performance_preview <- performance_preview %>%
  mutate(
    interpretation = "Diagnostic only; not a final Monte Carlo result"
  )

# 7. Assemble the QC results

structural_qc <- tibble::tibble(
  check = c(
    "Expected generating settings are present",
    "Exactly four rows per setting-repetition",
    "Exactly the four prespecified specification labels",
    "No duplicate setting-repetition-specification keys",
    "One diagnostic row per setting-repetition",
    "Result and diagnostic setting-repetition keys agree",
    "Both-correct exactly reproduces primary AIPW",
    "Selection misspecification is non-degenerate",
    "Augmentation misspecification is non-degenerate",
    "No misspecification-diagnostic failures"
  ),
  
  pass = c(
    settings_match,
    nrow(bad_key_counts) == 0L,
    nrow(bad_spec_sets) == 0L,
    nrow(duplicate_rows) == 0L,
    nrow(bad_diagnostic_key_counts) == 0L,
    diagnostic_keys_match,
    primary_exact_pass,
    nondegenerate_selection,
    nondegenerate_augmentation,
    no_diagnostic_failures
  )
)

# 8. Write QC outputs

write_csv(
  structural_qc,
  file.path(
    canary_dir,
    "aipw_robustness_canary_qc.csv"
  )
)

write_csv(
  failure_summary,
  file.path(
    canary_dir,
    "aipw_robustness_canary_failures.csv"
  )
)

write_csv(
  misspecification_summary,
  file.path(
    canary_dir,
    "aipw_robustness_canary_misspecification_summary.csv"
  )
)

write_csv(
  performance_preview,
  file.path(
    canary_dir,
    "aipw_robustness_canary_performance_preview.csv"
  )
)

# Save any problematic keys for investigation
if (nrow(bad_key_counts) > 0L) {
  write_csv(
    bad_key_counts,
    file.path(
      canary_dir,
      "aipw_robustness_canary_bad_row_counts.csv"
    )
  )
}

if (nrow(bad_spec_sets) > 0L) {
  write_csv(
    bad_spec_sets,
    file.path(
      canary_dir,
      "aipw_robustness_canary_bad_specification_sets.csv"
    )
  )
}

if (nrow(duplicate_rows) > 0L) {
  write_csv(
    duplicate_rows,
    file.path(
      canary_dir,
      "aipw_robustness_canary_duplicate_keys.csv"
    )
  )
}

# 9. Console report

cat(
  "\n============================================================\n",
  "AIPW ROBUSTNESS CANARY QC\n",
  "============================================================\n\n",
  sep = ""
)

print(structural_qc)

cat("\nExact comparison with primary AIPW:\n")
print(primary_match)

cat("\nAIPW failures:\n")
print(failure_summary, n = Inf)

cat("\nMisspecification diagnostics:\n")
print(misspecification_summary, n = Inf)

cat(
  "\nPerformance preview ",
  "(DIAGNOSTIC ONLY; canary sample is small):\n",
  sep = ""
)

print(performance_preview, n = Inf)

# Stop production when structural validation fails
if (!all(structural_qc$pass)) {
  cat(
    "\nOne or more structural QC checks failed. ",
    "Do NOT submit production yet.\n"
  )
  
  quit(
    save = "no",
    status = 2L
  )
}

# Stop production when any AIPW specification fails
if (!no_aipw_failures) {
  cat(
    "\nOne or more AIPW fits failed. ",
    "Review the failure output before production.\n"
  )
  
  quit(
    save = "no",
    status = 3L
  )
}

cat(
  "\nStructural canary QC passed.\n",
  "No AIPW failures were detected.\n",
  "NEXT: inspect the misspecification magnitudes before production.\n",
  sep = ""
)
