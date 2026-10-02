################################################################################
# 15_summarize_aipw_robustness.R
#
# Summarizes the AIPW double-robustness experiment using the ADEMP framework.
#
# Main performance measures:
#   - bias
#   - empirical standard error
#   - average model-based standard error
#   - model-SE to empirical-SE ratio
#   - mean squared error
#   - 95% confidence-interval coverage
#   - failure rate
#
# The script also summarizes paired differences from "Both correct",
# misspecification strength, and exact reproduction of primary AIPW results.
################################################################################

source("10_aipw_robustness_helpers.R")

library(dplyr)
library(tidyr)
library(readr)
library(rsimsum)

# Return NA instead of NaN or infinity when no finite values are available.
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

# Set the combined-input and summary-output directories.
robustness_root <- file.path(
  results_dir,
  "robustness",
  "aipw_double_robustness"
)

combined_robustness_dir <- file.path(
  robustness_root,
  "combined"
)

robustness_summary_dir <- file.path(
  robustness_root,
  "summary"
)

dir.create(
  robustness_summary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

estimate_file <- file.path(
  combined_robustness_dir,
  "aipw_robustness_estimates.rds"
)

diagnostic_file <- file.path(
  combined_robustness_dir,
  "aipw_robustness_misspecification_diagnostics.rds"
)

if (!file.exists(estimate_file)) {
  stop(
    "Combined robustness estimates not found:\n",
    estimate_file,
    "\nRun 14_combine_aipw_robustness.R first."
  )
}

if (!file.exists(diagnostic_file)) {
  stop(
    "Combined misspecification diagnostics not found:\n",
    diagnostic_file,
    "\nRun 14_combine_aipw_robustness.R first."
  )
}

robustness_results <- readRDS(estimate_file)
misspecification_diagnostics <- readRDS(diagnostic_file)

# Check the columns needed for performance summaries.
required_result_columns <- c(
  "setting_id",
  "setting",
  "repetition",
  "specification",
  "status",
  "estimate",
  "se",
  "conf_low",
  "conf_high",
  "estimate_difference_from_primary",
  "se_difference_from_primary"
)

required_diagnostic_columns <- c(
  "setting_id",
  "setting",
  "repetition",
  "mean_abs_pi_difference",
  "max_abs_pi_difference",
  "rmse_pi_difference",
  "cor_pi_correct_misspecified",
  "mean_abs_marker_mean_difference",
  "max_abs_marker_mean_difference",
  "rmse_marker_mean_difference",
  "relative_frobenius_sigma_difference"
)

missing_result_columns <- setdiff(
  required_result_columns,
  names(robustness_results)
)

missing_diagnostic_columns <- setdiff(
  required_diagnostic_columns,
  names(misspecification_diagnostics)
)

if (length(missing_result_columns) > 0L) {
  stop(
    "Combined robustness results are missing: ",
    paste(missing_result_columns, collapse = ", "),
    "."
  )
}

if (length(missing_diagnostic_columns) > 0L) {
  stop(
    "Combined diagnostics are missing: ",
    paste(missing_diagnostic_columns, collapse = ", "),
    "."
  )
}

specification_order <- c(
  "Both correct",
  "Selection correct only",
  "Augmentation correct only",
  "Both misspecified"
)

unknown_specifications <- setdiff(
  unique(as.character(robustness_results$specification)),
  specification_order
)

if (length(unknown_specifications) > 0L) {
  stop(
    "Unexpected AIPW specifications found: ",
    paste(unknown_specifications, collapse = ", "),
    "."
  )
}

robustness_results <- robustness_results %>%
  mutate(
    specification = factor(
      specification,
      levels = specification_order
    ),
    valid_result = (
      !is.na(status) &
        status == "ok" &
        is.finite(estimate) &
        is.finite(se)
    )
  )

# Exclude failed fits from performance and paired summaries.
analysis_results <- robustness_results %>%
  mutate(
    estimate = if_else(
      valid_result,
      estimate,
      NA_real_
    ),
    se = if_else(
      valid_result,
      se,
      NA_real_
    ),
    conf_low = if_else(
      valid_result & is.finite(conf_low),
      conf_low,
      NA_real_
    ),
    conf_high = if_else(
      valid_result & is.finite(conf_high),
      conf_high,
      NA_real_
    )
  )

# Summarize Monte Carlo performance separately for each setting.
statistics_to_keep <- c(
  "nsim",
  "thetamean",
  "bias",
  "empse",
  "mse",
  "modelse",
  "relerror",
  "cover",
  "becover"
)

rsimsum_objects <- vector(
  "list",
  nrow(aipw_robustness_settings)
)

rsimsum_tables <- vector(
  "list",
  nrow(aipw_robustness_settings)
)

for (i in seq_len(nrow(aipw_robustness_settings))) {
  setting <- aipw_robustness_settings[i, , drop = FALSE]
  
  setting_data <- analysis_results %>%
    filter(.data$setting_id == .env$setting$setting_id)
  
  if (nrow(setting_data) == 0L) {
    stop(
      "No robustness results were found for setting ",
      setting$setting_id,
      "."
    )
  }
  
  rsimsum_object <- rsimsum::simsum(
    data = setting_data,
    estvarname = "estimate",
    true = setting$theta,
    se = "se",
    methodvar = "specification",
    ref = "Both correct",
    ci.limits = c("conf_low", "conf_high"),
    dropbig = FALSE,
    x = FALSE,
    control = list(
      mcse = TRUE,
      level = nominal_level,
      na.rm = TRUE
    )
  )
  
  rsimsum_summary <- summary(
    rsimsum_object,
    stats = statistics_to_keep
  )
  
  rsimsum_tables[[i]] <- generics::tidy(
    rsimsum_summary
  ) %>%
    mutate(
      setting_id = .env$setting$setting_id,
      setting = .env$setting$setting,
      N = .env$setting$N,
      phase2_fraction = .env$setting$phase2_fraction,
      r2_a_marker = .env$setting$r2_a_marker,
      r2_y_marker = .env$setting$r2_y_marker,
      theta_true = .env$setting$theta
    ) %>%
    relocate(
      setting_id,
      setting,
      N,
      phase2_fraction,
      r2_a_marker,
      r2_y_marker,
      theta_true
    )
  
  rsimsum_objects[[i]] <- rsimsum_object
  
  setting_slug <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    setting$setting
  )
  
  saveRDS(
    rsimsum_object,
    file.path(
      robustness_summary_dir,
      sprintf(
        "aipw_robustness_%s_rsimsum.rds",
        setting_slug
      )
    )
  )
}

rsimsum_tidy <- bind_rows(rsimsum_tables)

# Standardize the method-column name across rsimsum versions.
if ("specification" %in% names(rsimsum_tidy)) {
  rsimsum_tidy <- rsimsum_tidy %>%
    rename(method = specification)
} else if (!("method" %in% names(rsimsum_tidy))) {
  stop(
    "Could not identify the rsimsum method column. Columns found: ",
    paste(names(rsimsum_tidy), collapse = ", "),
    "."
  )
}

write_csv(
  rsimsum_tidy,
  file.path(
    robustness_summary_dir,
    "aipw_robustness_rsimsum_tidy.csv"
  )
)

# Reshape the rsimsum results to one row per setting and specification.
performance_wide <- rsimsum_tidy %>%
  select(
    setting_id,
    setting,
    N,
    phase2_fraction,
    r2_a_marker,
    r2_y_marker,
    theta_true,
    method,
    stat,
    est,
    mcse
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = c(est, mcse),
    names_sep = "__"
  )

required_performance_columns <- c(
  "est__bias",
  "est__empse",
  "est__modelse",
  "est__relerror",
  "est__mse",
  "est__cover"
)

missing_performance_columns <- setdiff(
  required_performance_columns,
  names(performance_wide)
)

if (length(missing_performance_columns) > 0L) {
  stop(
    "The rsimsum output is missing: ",
    paste(missing_performance_columns, collapse = ", "),
    "."
  )
}

performance_wide <- performance_wide %>%
  mutate(
    se_ratio = 1 + est__relerror / 100,
    mcse_se_ratio = mcse__relerror / 100,
    method = factor(
      method,
      levels = specification_order
    )
  ) %>%
  arrange(setting_id, method)

write_csv(
  performance_wide,
  file.path(
    robustness_summary_dir,
    "aipw_robustness_performance_wide.csv"
  )
)

# Calculate failures using all attempted fits as the denominator.
failure_summary <- robustness_results %>%
  group_by(
    setting_id,
    setting,
    specification
  ) %>%
  summarise(
    nsim_total = n(),
    n_failed = sum(!valid_result),
    failure_rate = n_failed / nsim_total,
    mcse_failure = sqrt(
      failure_rate * (1 - failure_rate) / nsim_total
    ),
    .groups = "drop"
  )

write_csv(
  failure_summary,
  file.path(
    robustness_summary_dir,
    "aipw_robustness_failure_rates.csv"
  )
)

# Calculate paired differences from "Both correct" on the same datasets.
paired_results <- analysis_results %>%
  select(
    setting_id,
    setting,
    repetition,
    specification,
    estimate,
    se
  ) %>%
  pivot_wider(
    names_from = specification,
    values_from = c(estimate, se)
  )

paired_summary_rows <- vector(
  "list",
  length(specification_order) - 1L
)

for (i in seq_along(specification_order[-1])) {
  specification_name <- specification_order[-1][i]
  
  estimate_column <- paste0(
    "estimate_",
    specification_name
  )
  
  se_column <- paste0(
    "se_",
    specification_name
  )
  
  paired_summary_rows[[i]] <- paired_results %>%
    transmute(
      setting_id,
      setting,
      specification = specification_name,
      estimate_difference = (
        .data[[estimate_column]] -
          .data[["estimate_Both correct"]]
      ),
      se_difference = (
        .data[[se_column]] -
          .data[["se_Both correct"]]
      )
    ) %>%
    group_by(
      setting_id,
      setting,
      specification
    ) %>%
    summarise(
      n_estimate_pairs = sum(
        is.finite(estimate_difference)
      ),
      n_se_pairs = sum(
        is.finite(se_difference)
      ),
      mean_estimate_difference = safe_mean(
        estimate_difference
      ),
      sd_estimate_difference = safe_sd(
        estimate_difference
      ),
      median_estimate_difference = safe_median(
        estimate_difference
      ),
      mean_abs_estimate_difference = safe_mean(
        abs(estimate_difference)
      ),
      max_abs_estimate_difference = safe_max(
        abs(estimate_difference)
      ),
      mean_se_difference = safe_mean(
        se_difference
      ),
      mean_abs_se_difference = safe_mean(
        abs(se_difference)
      ),
      .groups = "drop"
    )
}

paired_summary <- bind_rows(paired_summary_rows)

write_csv(
  paired_summary,
  file.path(
    robustness_summary_dir,
    "aipw_robustness_paired_differences.csv"
  )
)

# Summarize how strongly the nuisance models were misspecified.
misspecification_summary <- misspecification_diagnostics %>%
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
    avg_rmse_pi_difference = safe_mean(
      rmse_pi_difference
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
    avg_rmse_marker_mean_difference = safe_mean(
      rmse_marker_mean_difference
    ),
    avg_relative_frobenius_sigma_difference = safe_mean(
      relative_frobenius_sigma_difference
    ),
    .groups = "drop"
  )

write_csv(
  misspecification_summary,
  file.path(
    robustness_summary_dir,
    "aipw_robustness_misspecification_summary.csv"
  )
)

# Confirm exact primary reproduction across both settings.
primary_reproduction_qc <- robustness_results %>%
  filter(.data$specification == "Both correct") %>%
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
      n_compared == nsim_primary &
        is.finite(max_abs_estimate_difference) &
        is.finite(max_abs_se_difference) &
        max_abs_estimate_difference < 1e-10 &
        max_abs_se_difference < 1e-10
    )
  )

write_csv(
  primary_reproduction_qc,
  file.path(
    robustness_summary_dir,
    "aipw_robustness_primary_reproduction_qc.csv"
  )
)

cat(
  "\nAIPW double-robustness summary complete.\n",
  "Summary directory: ", robustness_summary_dir, "\n\n",
  sep = ""
)

print(
  performance_wide %>%
    select(
      setting,
      method,
      est__bias,
      est__empse,
      est__modelse,
      se_ratio,
      est__mse,
      est__cover
    ),
  n = Inf
)

cat("\nPrimary reproduction QC:\n")
print(primary_reproduction_qc, n = Inf)

cat("\nFailure rates:\n")
print(failure_summary, n = Inf)

cat("\nMisspecification summary:\n")
print(misspecification_summary, n = Inf)
