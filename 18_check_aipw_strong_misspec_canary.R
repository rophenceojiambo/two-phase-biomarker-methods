################################################################################
# 18_check_aipw_strong_misspec_canary.R
#
# Validates the strong-misspecification AIPW canary.
#
# The checks confirm:
#   1. four specifications were run on every generated dataset
#   2. "Both correct" exactly reproduces the primary AIPW result
#   3. the X-only nuisance models differ from the correct models
#   4. diagnostic calculations and AIPW fits completed successfully
#
# The performance summaries are descriptive only because the canary is small.
################################################################################

source("16_aipw_strong_misspec_helpers.R")

library(dplyr)
library(readr)

# Return NA when a group contains no finite values.
safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else median(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else sd(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

# Set the canary input and output paths.
sensitivity_root <- file.path(
  results_dir,
  "sensitivity",
  "aipw_strong_misspecification"
)

canary_dir <- file.path(
  sensitivity_root,
  "canary"
)

results_file <- file.path(
  canary_dir,
  "aipw_strong_canary_results.csv"
)

diagnostics_file <- file.path(
  canary_dir,
  "aipw_strong_canary_misspecification_diagnostics.csv"
)

if (!file.exists(results_file)) {
  stop(
    "Strong canary results not found:\n",
    results_file,
    "\nRun 17_run_aipw_strong_misspec_canary.R first."
  )
}

if (!file.exists(diagnostics_file)) {
  stop(
    "Strong canary diagnostics not found:\n",
    diagnostics_file,
    "\nRun 17_run_aipw_strong_misspec_canary.R first."
  )
}

canary_results <- read_csv(
  results_file,
  show_col_types = FALSE
)

canary_diagnostics <- read_csv(
  diagnostics_file,
  show_col_types = FALSE
)

required_result_columns <- c(
  "setting_id",
  "setting",
  "repetition",
  "n_canary",
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
)

required_diagnostic_columns <- c(
  "setting_id",
  "setting",
  "repetition",
  "n_canary",
  "mean_abs_pi_difference",
  "max_abs_pi_difference",
  "cor_pi_correct_misspecified",
  "mean_abs_marker_mean_difference",
  "max_abs_marker_mean_difference",
  "relative_frobenius_sigma_difference",
  "diagnostic_message"
)

missing_result_columns <- setdiff(
  required_result_columns,
  names(canary_results)
)

missing_diagnostic_columns <- setdiff(
  required_diagnostic_columns,
  names(canary_diagnostics)
)

if (length(missing_result_columns) > 0L) {
  stop(
    "Canary results are missing: ",
    paste(missing_result_columns, collapse = ", "),
    "."
  )
}

if (length(missing_diagnostic_columns) > 0L) {
  stop(
    "Canary diagnostics are missing: ",
    paste(missing_diagnostic_columns, collapse = ", "),
    "."
  )
}

n_canary_values <- unique(canary_results$n_canary)

if (
  length(n_canary_values) != 1L ||
  is.na(n_canary_values)
) {
  stop("Could not determine a unique canary repetition count.")
}

n_canary <- as.integer(n_canary_values)
expected_repetitions <- seq_len(n_canary)
expected_setting_ids <- aipw_strong_settings$setting_id
expected_specifications <- sort(
  aipw_strong_specs$specification
)

expected_specification_string <- paste(
  expected_specifications,
  collapse = " | "
)

# Construct the expected setting-repetition keys.
expected_keys <- expand.grid(
  setting_id = expected_setting_ids,
  repetition = expected_repetitions,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

expected_keys$expected <- TRUE

# Check for exactly four result rows per generated dataset.
result_key_counts <- canary_results %>%
  count(setting_id, repetition, name = "n_rows")

result_key_check <- full_join(
  expected_keys,
  result_key_counts,
  by = c("setting_id", "repetition")
) %>%
  filter(
    is.na(expected) |
      is.na(n_rows) |
      n_rows != length(expected_specifications)
  )

# Check the specification labels within each dataset.
specification_check <- canary_results %>%
  group_by(setting_id, repetition) %>%
  summarise(
    n_specifications = n_distinct(specification),
    specifications = paste(
      sort(unique(specification)),
      collapse = " | "
    ),
    .groups = "drop"
  )

incorrect_specification_sets <- specification_check %>%
  filter(
    n_specifications !=
      length(expected_specifications) |
      specifications !=
      expected_specification_string
  )

duplicate_result_rows <- canary_results %>%
  count(
    setting_id,
    repetition,
    specification,
    name = "n"
  ) %>%
  filter(n != 1L)

# Check for exactly one diagnostic row per dataset.
diagnostic_key_counts <- canary_diagnostics %>%
  count(setting_id, repetition, name = "n_rows")

diagnostic_key_check <- full_join(
  expected_keys,
  diagnostic_key_counts,
  by = c("setting_id", "repetition")
) %>%
  filter(
    is.na(expected) |
      is.na(n_rows) |
      n_rows != 1L
  )

diagnostic_failures <- canary_diagnostics %>%
  filter(
    !is.na(diagnostic_message) &
      nzchar(diagnostic_message)
  )

# Check exact reproduction separately within each setting.
primary_match <- canary_results %>%
  filter(specification == "Both correct") %>%
  group_by(setting_id, setting) %>%
  summarise(
    n_compared = sum(
      is.finite(estimate_difference_from_primary) &
        is.finite(se_difference_from_primary)
    ),
    max_abs_estimate_difference = safe_max(
      abs(estimate_difference_from_primary)
    ),
    max_abs_se_difference = safe_max(
      abs(se_difference_from_primary)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    exact_reproduction = (
      n_compared == n_canary &
        is.finite(max_abs_estimate_difference) &
        is.finite(max_abs_se_difference) &
        max_abs_estimate_difference < 1e-10 &
        max_abs_se_difference < 1e-10
    )
  )

primary_exact_pass <- (
  nrow(primary_match) ==
    nrow(aipw_strong_settings) &&
    all(primary_match$exact_reproduction)
)

# Count failed AIPW fits.
canary_results <- canary_results %>%
  mutate(
    valid_result = (
      !is.na(status) &
        status == "ok" &
        is.finite(estimate) &
        is.finite(se)
    )
  )

failure_summary <- canary_results %>%
  group_by(setting_id, setting, specification) %>%
  summarise(
    n = n(),
    n_failed = sum(!valid_result),
    failure_rate = n_failed / n,
    .groups = "drop"
  )

# Summarize the size of the imposed misspecification.
misspecification_summary <- canary_diagnostics %>%
  group_by(setting_id, setting) %>%
  summarise(
    repetitions = n(),
    avg_mean_abs_pi_difference = safe_mean(
      mean_abs_pi_difference
    ),
    median_mean_abs_pi_difference = safe_median(
      mean_abs_pi_difference
    ),
    avg_max_abs_pi_difference = safe_mean(
      max_abs_pi_difference
    ),
    avg_cor_pi_correct_misspecified = safe_mean(
      cor_pi_correct_misspecified
    ),
    avg_mean_abs_marker_mean_difference = safe_mean(
      mean_abs_marker_mean_difference
    ),
    median_mean_abs_marker_mean_difference = safe_median(
      mean_abs_marker_mean_difference
    ),
    avg_max_abs_marker_mean_difference = safe_mean(
      max_abs_marker_mean_difference
    ),
    avg_relative_frobenius_sigma_difference = safe_mean(
      relative_frobenius_sigma_difference
    ),
    .groups = "drop"
  )

nondegenerate_selection <- (
  nrow(misspecification_summary) ==
    nrow(aipw_strong_settings) &&
    all(
      is.finite(
        misspecification_summary$
          avg_mean_abs_pi_difference
      )
    ) &&
    all(
      misspecification_summary$
        avg_mean_abs_pi_difference > 1e-6
    )
)

nondegenerate_augmentation <- (
  nrow(misspecification_summary) ==
    nrow(aipw_strong_settings) &&
    all(
      is.finite(
        misspecification_summary$
          avg_mean_abs_marker_mean_difference
      )
    ) &&
    all(
      misspecification_summary$
        avg_mean_abs_marker_mean_difference > 1e-6
    )
)

# Provide a descriptive performance preview only.
performance_preview <- canary_results %>%
  group_by(
    setting_id,
    setting,
    specification,
    selection_correct,
    augmentation_correct
  ) %>%
  summarise(
    nsim = sum(valid_result),
    mean_estimate = safe_mean(
      estimate[valid_result]
    ),
    bias = (
      mean_estimate -
        first(theta_true)
    ),
    empirical_se = safe_sd(
      estimate[valid_result]
    ),
    average_model_se = safe_mean(
      se[valid_result]
    ),
    se_ratio = (
      average_model_se /
        empirical_se
    ),
    mse = safe_mean(
      (
        estimate[valid_result] -
          first(theta_true)
      )^2
    ),
    coverage = safe_mean(
      as.numeric(
        conf_low[valid_result] <=
          first(theta_true) &
          conf_high[valid_result] >=
          first(theta_true)
      )
    ),
    mean_elapsed_seconds = safe_mean(
      elapsed_seconds
    ),
    .groups = "drop"
  )

# Collect all structural checks in one table.
structural_qc <- tibble::tibble(
  check = c(
    "Exactly four rows per setting-repetition",
    "Exactly the four prespecified specification labels",
    "No duplicate result keys",
    "Exactly one diagnostic row per setting-repetition",
    "No diagnostic calculation failures",
    "Both-correct exactly reproduces primary AIPW",
    "Strong selection misspecification is non-degenerate",
    "Strong augmentation misspecification is non-degenerate"
  ),
  pass = c(
    nrow(result_key_check) == 0L,
    nrow(incorrect_specification_sets) == 0L,
    nrow(duplicate_result_rows) == 0L,
    nrow(diagnostic_key_check) == 0L,
    nrow(diagnostic_failures) == 0L,
    primary_exact_pass,
    nondegenerate_selection,
    nondegenerate_augmentation
  )
)

# Save the canary QC outputs.
write_csv(
  structural_qc,
  file.path(
    canary_dir,
    "aipw_strong_canary_qc.csv"
  )
)

write_csv(
  failure_summary,
  file.path(
    canary_dir,
    "aipw_strong_canary_failures.csv"
  )
)

write_csv(
  misspecification_summary,
  file.path(
    canary_dir,
    "aipw_strong_canary_misspecification_summary.csv"
  )
)

write_csv(
  performance_preview,
  file.path(
    canary_dir,
    "aipw_strong_canary_performance_preview.csv"
  )
)

cat(
  "\nStrong AIPW misspecification canary QC\n\n"
)

print(structural_qc)

cat("\nExact comparison with primary AIPW:\n")
print(primary_match, n = Inf)

cat("\nFailures:\n")
print(failure_summary, n = Inf)

cat("\nStrong misspecification diagnostics:\n")
print(misspecification_summary, n = Inf)

cat(
  "\nPerformance preview ",
  "(diagnostic only; the canary is small):\n"
)

print(performance_preview, n = Inf)

if (!all(structural_qc$pass)) {
  cat(
    "\nOne or more structural QC checks failed. ",
    "Do not run the production sensitivity.\n"
  )
  
  quit(save = "no", status = 2L)
}

if (any(failure_summary$n_failed > 0L)) {
  cat(
    "\nOne or more AIPW fits failed. ",
    "Review the failures before production.\n"
  )
  
  quit(save = "no", status = 3L)
}

cat(
  "\nStrong-misspecification canary QC passed.\n",
  "Inspect the misspecification magnitudes before production.\n",
  sep = ""
)
