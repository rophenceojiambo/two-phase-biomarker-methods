################################################################################
# 27_summarize_empirical_residual_production.R
#
# Summarizes the four empirical-residual production scenarios using rsimsum.
#
# The script produces ADEMP performance summaries, Type I error, failures,
# runtime and design diagnostics, and comparisons with the matched primary
# multivariate-normal DGMs.
################################################################################

source("00_config.R")
source("22_empirical_residual_sensitivity_helpers.R")

library(dplyr)
library(tidyr)
library(readr)
library(rsimsum)

all_true <- function(x) {
  isTRUE(all(x))
}

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

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    NA_real_
  } else {
    unname(
      quantile(
        x,
        probability,
        names = FALSE
      )
    )
  }
}

expected_methods <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

# Set combined-input and summary-output paths.
empirical_results_dir <- file.path(
  results_dir,
  "empirical_residual_sensitivity"
)

empirical_production_dir <- file.path(
  empirical_results_dir,
  "production"
)

empirical_combined_dir <- file.path(
  empirical_production_dir,
  "combined"
)

empirical_scenario_dir <- file.path(
  empirical_combined_dir,
  "scenarios"
)

empirical_summary_dir <- file.path(
  empirical_production_dir,
  "summary"
)

empirical_rsimsum_dir <- file.path(
  empirical_summary_dir,
  "rsimsum_objects"
)

invisible(
  lapply(
    c(
      empirical_summary_dir,
      empirical_rsimsum_dir
    ),
    dir.create,
    showWarnings = FALSE,
    recursive = TRUE
  )
)

summary_pass_file <- file.path(
  empirical_summary_dir,
  "EMPIRICAL_RESIDUAL_PRODUCTION_SUMMARY_PASS.txt"
)

# Prevent a previous pass file from surviving a failed summary rerun.
if (file.exists(summary_pass_file)) {
  unlink(summary_pass_file)
}

combine_pass_file <- file.path(
  empirical_combined_dir,
  "EMPIRICAL_RESIDUAL_PRODUCTION_COMBINE_PASS.txt"
)

production_qc_file <- file.path(
  empirical_combined_dir,
  "empirical_residual_production_qc.csv"
)

completion_manifest_file <- file.path(
  empirical_combined_dir,
  "empirical_residual_production_completion_manifest.csv"
)

canary_reproduction_file <- file.path(
  empirical_combined_dir,
  "empirical_residual_canary_reproduction_qc.csv"
)

required_combined_files <- c(
  combine_pass_file,
  production_qc_file,
  completion_manifest_file,
  canary_reproduction_file
)

missing_combined_files <- required_combined_files[
  !file.exists(required_combined_files)
]

if (length(missing_combined_files) > 0L) {
  stop(
    "Required combined-production files are missing:\n- ",
    paste(missing_combined_files, collapse = "\n- "),
    "\nRun 26_combine_empirical_residual_production.R first."
  )
}

production_qc <- read_csv(
  production_qc_file,
  show_col_types = FALSE
)

completion_manifest <- read_csv(
  completion_manifest_file,
  show_col_types = FALSE
)

canary_reproduction_qc <- read_csv(
  canary_reproduction_file,
  show_col_types = FALSE
)

if (
  !all(c(
    "passed",
    "overall_production_content_pass"
  ) %in% names(production_qc)) ||
  !all_true(production_qc$passed) ||
  !all_true(
    production_qc$
    overall_production_content_pass
  )
) {
  stop(
    "Combined empirical-production QC does not record a full pass."
  )
}

if (
  !all(c("complete", "problem") %in%
       names(completion_manifest)) ||
  nrow(completion_manifest) !=
  nrow(empirical_sensitivity_grid) ||
  !all_true(completion_manifest$complete) ||
  !all_true(
    is.na(completion_manifest$problem) |
      completion_manifest$problem == ""
  )
) {
  stop(
    "The empirical-production completion manifest is incomplete."
  )
}

if (
  !("reproduced_within_tolerance" %in%
    names(canary_reproduction_qc)) ||
  !all_true(
    canary_reproduction_qc$
    reproduced_within_tolerance
  )
) {
  stop(
    "Canary reproduction is not fully verified."
  )
}

statistics_to_keep <- c(
  "nsim",
  "thetamean",
  "bias",
  "empse",
  "mse",
  "relprec",
  "modelse",
  "relerror",
  "cover",
  "becover",
  "power"
)

n_scenarios <- nrow(empirical_sensitivity_grid)

rsimsum_tidy_list <- vector(
  "list",
  n_scenarios
)

type1_list <- list()
type1_index <- 0L

failure_list <- vector("list", n_scenarios)
runtime_list <- vector("list", n_scenarios)
phase2_list <- vector("list", n_scenarios)
weight_list <- vector("list", n_scenarios)
scenario_results_list <- vector(
  "list",
  n_scenarios
)

required_result_columns <- c(
  "repetition",
  "method",
  "status",
  "estimate",
  "se",
  "p_value",
  "conf_low",
  "conf_high",
  "elapsed_seconds",
  "realized_phase2_fraction",
  "n_phase2",
  "min_pi_true",
  "p01_pi_true",
  "median_pi_true",
  "p99_pi_true",
  "max_pi_true",
  "weight_min",
  "weight_p99",
  "weight_max",
  "weight_cv",
  "weight_ess"
)

# Summarize each empirical-residual scenario.
for (scenario_index in seq_len(n_scenarios)) {
  scenario <- empirical_sensitivity_grid[
    scenario_index,
    ,
    drop = FALSE
  ]
  
  scenario_id <- scenario$sensitivity_scenario_id
  theta_true <- scenario$theta[[1]]
  
  cat(
    "Summarizing empirical scenario ",
    scenario_id,
    " of ",
    n_scenarios,
    " with rsimsum.\n",
    sep = ""
  )
  
  scenario_file <- file.path(
    empirical_scenario_dir,
    sprintf(
      "empirical_scenario_%02d_estimates.rds",
      scenario_id
    )
  )
  
  if (!file.exists(scenario_file)) {
    stop(
      "Empirical scenario estimates not found:\n",
      scenario_file
    )
  }
  
  scenario_results <- readRDS(scenario_file)
  
  missing_result_columns <- setdiff(
    required_result_columns,
    names(scenario_results)
  )
  
  if (length(missing_result_columns) > 0L) {
    stop(
      "Scenario ",
      scenario_id,
      " is missing: ",
      paste(missing_result_columns, collapse = ", "),
      "."
    )
  }
  
  if (!setequal(
    scenario_results$method,
    expected_methods
  )) {
    stop(
      "Scenario ",
      scenario_id,
      " has unexpected or missing method labels."
    )
  }
  
  scenario_results <- scenario_results %>%
    mutate(
      valid_result = (
        !is.na(status) &
          status == "ok" &
          is.finite(estimate) &
          is.finite(se) &
          se > 0
      )
    )
  
  scenario_results_list[[scenario_index]] <-
    scenario_results
  
  # Failed fits are set to missing for performance calculations.
  analysis_results <- scenario_results %>%
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
  
  rsimsum_object <- rsimsum::simsum(
    data = analysis_results,
    estvarname = "estimate",
    true = theta_true,
    se = "se",
    methodvar = "method",
    ref = "CCA",
    ci.limits = c("conf_low", "conf_high"),
    dropbig = FALSE,
    x = FALSE,
    control = list(
      mcse = TRUE,
      level = nominal_level,
      na.rm = TRUE
    )
  )
  
  saveRDS(
    rsimsum_object,
    file.path(
      empirical_rsimsum_dir,
      sprintf(
        "empirical_scenario_%02d_rsimsum.rds",
        scenario_id
      )
    )
  )
  
  rsimsum_tidy_list[[scenario_index]] <- generics::tidy(
    summary(
      rsimsum_object,
      stats = statistics_to_keep
    )
  ) %>%
    mutate(
      sensitivity_scenario_id = scenario_id,
      primary_scenario_id =
        scenario$primary_scenario_id[[1]],
      scenario_label =
        scenario$scenario_label[[1]],
      scenario_key = scenario$scenario_key[[1]],
      N = scenario$N[[1]],
      phase2_fraction =
        scenario$phase2_fraction[[1]],
      r2_a_marker =
        scenario$r2_a_marker[[1]],
      r2_y_marker =
        scenario$r2_y_marker[[1]],
      theta_true = theta_true,
      marker_error = "empirical"
    ) %>%
    relocate(
      sensitivity_scenario_id,
      primary_scenario_id,
      scenario_label,
      scenario_key,
      N,
      phase2_fraction,
      r2_a_marker,
      r2_y_marker,
      theta_true,
      marker_error
    )
  
  # Calculate Type I error only for null scenarios.
  if (isTRUE(all.equal(theta_true, 0))) {
    type1_index <- type1_index + 1L
    
    type1_list[[type1_index]] <- scenario_results %>%
      group_by(method) %>%
      summarise(
        nsim_total = n(),
        nsim_success = sum(
          valid_result &
            is.finite(p_value)
        ),
        null_rejection_rate = if (
          nsim_success > 0L
        ) {
          mean(
            p_value[
              valid_result &
                is.finite(p_value)
            ] < alpha_level
          )
        } else {
          NA_real_
        },
        mcse_null_rejection = if (
          nsim_success > 0L &&
          is.finite(null_rejection_rate)
        ) {
          sqrt(
            null_rejection_rate *
              (1 - null_rejection_rate) /
              nsim_success
          )
        } else {
          NA_real_
        },
        .groups = "drop"
      )
  }
  
  failure_list[[scenario_index]] <- scenario_results %>%
    group_by(method) %>%
    summarise(
      nsim_total = n(),
      n_failed = sum(!valid_result),
      failure_rate = n_failed / nsim_total,
      mcse_failure = sqrt(
        failure_rate *
          (1 - failure_rate) /
          nsim_total
      ),
      .groups = "drop"
    )
  
  runtime_list[[scenario_index]] <- scenario_results %>%
    group_by(method) %>%
    summarise(
      mean_runtime_seconds = safe_mean(
        elapsed_seconds
      ),
      median_runtime_seconds = safe_median(
        elapsed_seconds
      ),
      p95_runtime_seconds = safe_quantile(
        elapsed_seconds,
        0.95
      ),
      total_cpu_hours = sum(
        elapsed_seconds[is.finite(elapsed_seconds)],
        na.rm = TRUE
      ) / 3600,
      .groups = "drop"
    )
  
  phase2_list[[scenario_index]] <- scenario_results %>%
    distinct(
      repetition,
      realized_phase2_fraction,
      n_phase2,
      min_pi_true,
      p01_pi_true,
      median_pi_true,
      p99_pi_true,
      max_pi_true
    ) %>%
    summarise(
      mean_n_phase2 = safe_mean(n_phase2),
      sd_n_phase2 = safe_sd(n_phase2),
      min_n_phase2 = safe_min(n_phase2),
      max_n_phase2 = {
        values <- n_phase2[is.finite(n_phase2)]
        if (length(values) == 0L) {
          NA_real_
        } else {
          max(values)
        }
      },
      mean_realized_phase2_fraction =
        safe_mean(realized_phase2_fraction),
      mean_min_pi_true = safe_mean(min_pi_true),
      mean_p01_pi_true = safe_mean(p01_pi_true),
      mean_median_pi_true =
        safe_mean(median_pi_true),
      mean_p99_pi_true = safe_mean(p99_pi_true),
      mean_max_pi_true = safe_mean(max_pi_true)
    )
  
  weight_list[[scenario_index]] <- scenario_results %>%
    filter(
      method %in% c("IPW", "AIPW"),
      valid_result
    ) %>%
    group_by(method) %>%
    summarise(
      mean_weight_min = safe_mean(weight_min),
      mean_weight_p99 = safe_mean(weight_p99),
      mean_weight_max = safe_mean(weight_max),
      mean_weight_cv = safe_mean(weight_cv),
      p95_weight_cv = safe_quantile(
        weight_cv,
        0.95
      ),
      mean_weight_ess = safe_mean(weight_ess),
      min_weight_ess = safe_min(weight_ess),
      .groups = "drop"
    )
  
  # Attach the scenario metadata to all descriptive summaries.
  add_scenario_metadata <- function(table) {
    table %>%
      mutate(
        sensitivity_scenario_id = scenario_id,
        primary_scenario_id =
          scenario$primary_scenario_id[[1]],
        scenario_label =
          scenario$scenario_label[[1]],
        scenario_key =
          scenario$scenario_key[[1]],
        N = scenario$N[[1]],
        phase2_fraction =
          scenario$phase2_fraction[[1]],
        r2_a_marker =
          scenario$r2_a_marker[[1]],
        r2_y_marker =
          scenario$r2_y_marker[[1]],
        theta_true = theta_true,
        marker_error = "empirical"
      ) %>%
      relocate(
        sensitivity_scenario_id,
        primary_scenario_id,
        scenario_label,
        scenario_key,
        N,
        phase2_fraction,
        r2_a_marker,
        r2_y_marker,
        theta_true,
        marker_error
      )
  }
  
  failure_list[[scenario_index]] <-
    add_scenario_metadata(
      failure_list[[scenario_index]]
    )
  
  runtime_list[[scenario_index]] <-
    add_scenario_metadata(
      runtime_list[[scenario_index]]
    )
  
  phase2_list[[scenario_index]] <-
    add_scenario_metadata(
      phase2_list[[scenario_index]]
    )
  
  weight_list[[scenario_index]] <-
    add_scenario_metadata(
      weight_list[[scenario_index]]
    )
  
  if (isTRUE(all.equal(theta_true, 0))) {
    type1_list[[type1_index]] <-
      add_scenario_metadata(
        type1_list[[type1_index]]
      )
  }
}

# Save empirical performance and diagnostic summaries.
rsimsum_tidy <- bind_rows(
  rsimsum_tidy_list
)

write_csv(
  rsimsum_tidy,
  file.path(
    empirical_summary_dir,
    "empirical_residual_rsimsum_tidy.csv"
  )
)

saveRDS(
  rsimsum_tidy,
  file.path(
    empirical_summary_dir,
    "empirical_residual_rsimsum_tidy.rds"
  ),
  compress = "xz"
)

performance_wide <- rsimsum_tidy %>%
  mutate(stat = as.character(stat)) %>%
  select(
    sensitivity_scenario_id,
    primary_scenario_id,
    scenario_label,
    scenario_key,
    N,
    phase2_fraction,
    r2_a_marker,
    r2_y_marker,
    theta_true,
    marker_error,
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
  "est__thetamean",
  "est__bias",
  "est__empse",
  "est__modelse",
  "est__relerror",
  "mcse__relerror",
  "est__cover",
  "est__becover",
  "est__mse"
)

missing_performance_columns <- setdiff(
  required_performance_columns,
  names(performance_wide)
)

if (length(missing_performance_columns) > 0L) {
  stop(
    "The empirical rsimsum output is missing: ",
    paste(
      missing_performance_columns,
      collapse = ", "
    ),
    "."
  )
}

performance_wide <- performance_wide %>%
  mutate(
    se_ratio = 1 + est__relerror / 100,
    mcse_se_ratio = mcse__relerror / 100
  )

type1_summary <- bind_rows(type1_list)
failure_summary <- bind_rows(failure_list)
runtime_summary <- bind_rows(runtime_list)
phase2_summary <- bind_rows(phase2_list)
weight_summary <- bind_rows(weight_list)

write_csv(
  performance_wide,
  file.path(
    empirical_summary_dir,
    "empirical_residual_performance_wide.csv"
  )
)

write_csv(
  type1_summary,
  file.path(
    empirical_summary_dir,
    "empirical_residual_type1_error.csv"
  )
)

write_csv(
  failure_summary,
  file.path(
    empirical_summary_dir,
    "empirical_residual_failure_rates.csv"
  )
)

write_csv(
  runtime_summary,
  file.path(
    empirical_summary_dir,
    "empirical_residual_runtime.csv"
  )
)

write_csv(
  phase2_summary,
  file.path(
    empirical_summary_dir,
    "empirical_residual_phase2_diagnostics.csv"
  )
)

write_csv(
  weight_summary,
  file.path(
    empirical_summary_dir,
    "empirical_residual_weight_diagnostics.csv"
  )
)

# Read the matched primary multivariate-normal summaries.
primary_performance_file <- file.path(
  summary_dir,
  "primary_performance_wide.csv"
)

primary_type1_file <- file.path(
  summary_dir,
  "primary_type1_error.csv"
)

primary_failure_file <- file.path(
  summary_dir,
  "primary_failure_rates.csv"
)

required_primary_files <- c(
  primary_performance_file,
  primary_type1_file,
  primary_failure_file
)

missing_primary_files <- required_primary_files[
  !file.exists(required_primary_files)
]

if (length(missing_primary_files) > 0L) {
  stop(
    "Matched primary summary files are missing:\n- ",
    paste(missing_primary_files, collapse = "\n- "),
    "\nRun 07_summarize_primary_rsimsum.R first."
  )
}

primary_performance <- read_csv(
  primary_performance_file,
  show_col_types = FALSE
) %>%
  filter(
    scenario_id %in%
      empirical_sensitivity_grid$
      primary_scenario_id
  )

primary_metric_columns <- c(
  "est__thetamean",
  "est__bias",
  "est__empse",
  "est__modelse",
  "se_ratio",
  "est__cover",
  "est__becover",
  "est__mse"
)

missing_metrics <- setdiff(
  primary_metric_columns,
  intersect(
    names(primary_performance),
    names(performance_wide)
  )
)

if (length(missing_metrics) > 0L) {
  stop(
    "Required performance metrics are missing: ",
    paste(missing_metrics, collapse = ", "),
    "."
  )
}

primary_performance_keys <- primary_performance %>%
  count(scenario_id, method, name = "n")

if (
  nrow(primary_performance_keys) !=
  n_scenarios * length(expected_methods) ||
  !all_true(primary_performance_keys$n == 1L)
) {
  stop(
    "Matched primary performance rows are incomplete or duplicated."
  )
}

empirical_comparison_input <- performance_wide %>%
  select(
    sensitivity_scenario_id,
    primary_scenario_id,
    scenario_label,
    scenario_key,
    N,
    phase2_fraction,
    r2_a_marker,
    r2_y_marker,
    theta_true,
    method,
    all_of(primary_metric_columns)
  ) %>%
  rename_with(
    ~ paste0("empirical_", .x),
    all_of(primary_metric_columns)
  )

primary_comparison_input <- primary_performance %>%
  select(
    scenario_id,
    method,
    all_of(primary_metric_columns)
  ) %>%
  rename(
    primary_scenario_id = scenario_id
  ) %>%
  rename_with(
    ~ paste0("primary_", .x),
    all_of(primary_metric_columns)
  )

performance_comparison <- empirical_comparison_input %>%
  inner_join(
    primary_comparison_input,
    by = c(
      "primary_scenario_id",
      "method"
    )
  ) %>%
  mutate(
    delta_mean_estimate = (
      empirical_est__thetamean -
        primary_est__thetamean
    ),
    delta_bias = (
      empirical_est__bias -
        primary_est__bias
    ),
    delta_absolute_bias = (
      abs(empirical_est__bias) -
        abs(primary_est__bias)
    ),
    delta_empirical_se = (
      empirical_est__empse -
        primary_est__empse
    ),
    delta_model_se = (
      empirical_est__modelse -
        primary_est__modelse
    ),
    delta_se_ratio = (
      empirical_se_ratio -
        primary_se_ratio
    ),
    delta_coverage = (
      empirical_est__cover -
        primary_est__cover
    ),
    delta_bias_eliminated_coverage = (
      empirical_est__becover -
        primary_est__becover
    ),
    delta_mse = (
      empirical_est__mse -
        primary_est__mse
    )
  )

write_csv(
  performance_comparison,
  file.path(
    empirical_summary_dir,
    "empirical_vs_primary_performance_comparison.csv"
  )
)

# Compare Type I error in the two null scenarios.
primary_type1 <- read_csv(
  primary_type1_file,
  show_col_types = FALSE
) %>%
  filter(
    scenario_id %in%
      empirical_sensitivity_grid$
      primary_scenario_id
  ) %>%
  select(
    scenario_id,
    method,
    primary_null_rejection_rate =
      null_rejection_rate,
    primary_mcse_null_rejection =
      mcse_null_rejection
  ) %>%
  rename(
    primary_scenario_id = scenario_id
  )

type1_comparison <- type1_summary %>%
  select(
    sensitivity_scenario_id,
    primary_scenario_id,
    scenario_label,
    scenario_key,
    N,
    phase2_fraction,
    r2_a_marker,
    r2_y_marker,
    theta_true,
    method,
    empirical_null_rejection_rate =
      null_rejection_rate,
    empirical_mcse_null_rejection =
      mcse_null_rejection
  ) %>%
  inner_join(
    primary_type1,
    by = c(
      "primary_scenario_id",
      "method"
    )
  ) %>%
  mutate(
    delta_null_rejection_rate = (
      empirical_null_rejection_rate -
        primary_null_rejection_rate
    )
  )

write_csv(
  type1_comparison,
  file.path(
    empirical_summary_dir,
    "empirical_vs_primary_type1_comparison.csv"
  )
)

# Compare method failure rates in all four scenarios.
primary_failures <- read_csv(
  primary_failure_file,
  show_col_types = FALSE
) %>%
  filter(
    scenario_id %in%
      empirical_sensitivity_grid$
      primary_scenario_id
  ) %>%
  select(
    scenario_id,
    method,
    primary_n_failed = n_failed,
    primary_failure_rate = failure_rate,
    primary_mcse_failure = mcse_failure
  ) %>%
  rename(
    primary_scenario_id = scenario_id
  )

failure_comparison <- failure_summary %>%
  select(
    sensitivity_scenario_id,
    primary_scenario_id,
    scenario_label,
    scenario_key,
    N,
    phase2_fraction,
    r2_a_marker,
    r2_y_marker,
    theta_true,
    method,
    empirical_n_failed = n_failed,
    empirical_failure_rate = failure_rate,
    empirical_mcse_failure = mcse_failure
  ) %>%
  inner_join(
    primary_failures,
    by = c(
      "primary_scenario_id",
      "method"
    )
  ) %>%
  mutate(
    delta_failure_rate = (
      empirical_failure_rate -
        primary_failure_rate
    )
  )

write_csv(
  failure_comparison,
  file.path(
    empirical_summary_dir,
    "empirical_vs_primary_failure_comparison.csv"
  )
)

# Apply final summary checks.
expected_method_scenario_rows <- (
  n_scenarios * length(expected_methods)
)

n_null_scenarios <- sum(
  empirical_sensitivity_grid$theta == 0
)

expected_null_method_rows <- (
  n_null_scenarios * length(expected_methods)
)

key_performance_columns <- c(
  "est__bias",
  "est__empse",
  "est__modelse",
  "se_ratio",
  "est__cover",
  "est__becover",
  "est__mse"
)

finite_key_performance <- all_true(
  vapply(
    performance_wide[key_performance_columns],
    function(x) {
      all(is.finite(x))
    },
    logical(1)
  )
)

summary_qc <- data.frame(
  check = c(
    "All empirical scenarios were summarized",
    "All scenario-method performance rows are present",
    "All null-scenario method rows are present",
    "All key performance summaries are finite",
    "Matched primary performance comparison is complete",
    "Matched Type I error comparison is complete",
    "Matched failure comparison is complete",
    "Canary reproduction remains verified"
  ),
  hard_check = rep(TRUE, 8L),
  passed = c(
    all_true(
      vapply(
        scenario_results_list,
        function(x) {
          is.data.frame(x) && nrow(x) > 0L
        },
        logical(1)
      )
    ),
    nrow(performance_wide) ==
      expected_method_scenario_rows,
    nrow(type1_summary) ==
      expected_null_method_rows,
    finite_key_performance,
    nrow(performance_comparison) ==
      expected_method_scenario_rows,
    nrow(type1_comparison) ==
      expected_null_method_rows,
    nrow(failure_comparison) ==
      expected_method_scenario_rows,
    all_true(
      canary_reproduction_qc$
        reproduced_within_tolerance
    )
  ),
  stringsAsFactors = FALSE
)

summary_qc$overall_summary_pass <- all_true(
  summary_qc$passed[summary_qc$hard_check]
)

write_csv(
  summary_qc,
  file.path(
    empirical_summary_dir,
    "empirical_residual_summary_qc.csv"
  )
)

# Flag failure rates above the presentation threshold.
review_flags <- failure_summary %>%
  mutate(
    exceeds_prespecified_failure_threshold = (
      failure_rate >
        deterioration_thresholds$failure_rate
    )
  ) %>%
  filter(
    exceeds_prespecified_failure_threshold
  )

write_csv(
  review_flags,
  file.path(
    empirical_summary_dir,
    "empirical_residual_failure_review_flags.csv"
  )
)

# Save a compact archive of all final summary objects.
saveRDS(
  list(
    scenarios = empirical_sensitivity_grid,
    rsimsum_tidy = rsimsum_tidy,
    performance_wide = performance_wide,
    type1_error = type1_summary,
    failure_rates = failure_summary,
    runtime = runtime_summary,
    phase2_diagnostics = phase2_summary,
    weight_diagnostics = weight_summary,
    performance_comparison =
      performance_comparison,
    type1_comparison = type1_comparison,
    failure_comparison = failure_comparison,
    summary_qc = summary_qc,
    failure_review_flags = review_flags
  ),
  file.path(
    empirical_summary_dir,
    "empirical_residual_production_summary_complete.rds"
  ),
  compress = "xz"
)

if (!all_true(
  summary_qc$passed[summary_qc$hard_check]
)) {
  failed_checks <- summary_qc$check[
    summary_qc$hard_check &
      !summary_qc$passed
  ]
  
  stop(
    "Empirical residual production summary failed. Failed checks:\n- ",
    paste(failed_checks, collapse = "\n- "),
    "\nSee outputs in: ",
    empirical_summary_dir
  )
}

writeLines(
  c(
    "EMPIRICAL RESIDUAL PRODUCTION SUMMARY: PASS",
    paste("Completed:", Sys.time()),
    paste("Production scenarios:", n_scenarios),
    paste(
      "Scenario-method performance rows:",
      nrow(performance_wide)
    ),
    paste(
      "Null-scenario method rows:",
      nrow(type1_summary)
    ),
    paste(
      "Failure-rate review flags:",
      nrow(review_flags)
    )
  ),
  summary_pass_file
)

cat(
  "\nEmpirical residual production summary passed.\n",
  "Main comparison: ",
  file.path(
    empirical_summary_dir,
    "empirical_vs_primary_performance_comparison.csv"
  ),
  "\nSummary directory: ",
  empirical_summary_dir,
  "\n",
  sep = ""
)

