################################################################################
# 33_summarize_phase2_mcar_results.R
#
# ADEMP performance summary for the four Phase-2 MCAR DGMs.
# Uses rsimsum as the authoritative summary engine and compares each MCAR
# scenario with its matched completed primary multivariate-normal DGM.
################################################################################

source("00_config.R")
source("28_phase2_mcar_sensitivity_helpers.R")

library(dplyr)
library(tidyr)
library(readr)
library(rsimsum)


# ------------------------------------------------------------------------------
# 1. Paths and required completed inputs
# ------------------------------------------------------------------------------

mcar_results_dir <- file.path(
  results_dir,
  "phase2_mcar_sensitivity"
)

mcar_production_dir <- file.path(
  mcar_results_dir,
  "production"
)

mcar_combined_dir <- file.path(
  mcar_production_dir,
  "combined"
)

mcar_scenario_dir <- file.path(
  mcar_combined_dir,
  "scenarios"
)

mcar_summary_dir <- file.path(
  mcar_production_dir,
  "summary"
)

mcar_rsimsum_dir <- file.path(
  mcar_summary_dir,
  "rsimsum_objects"
)

summary_pass_file <- file.path(
  mcar_summary_dir,
  "PHASE2_MCAR_PRODUCTION_SUMMARY_PASS.txt"
)

dir.create(
  mcar_summary_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

if (file.exists(summary_pass_file)) {
  unlink(summary_pass_file)
}

dir.create(
  mcar_rsimsum_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

combine_pass_file <- file.path(
  mcar_combined_dir,
  "PHASE2_MCAR_PRODUCTION_COMBINE_PASS.txt"
)

production_qc_file <- file.path(
  mcar_combined_dir,
  "phase2_mcar_production_qc.csv"
)

canary_reproduction_file <- file.path(
  mcar_combined_dir,
  "phase2_mcar_canary_reproduction_qc.csv"
)

primary_naive_reproduction_file <- file.path(
  mcar_combined_dir,
  "phase2_mcar_primary_naive_reproduction_qc.csv"
)

required_combined_files <- c(
  combine_pass_file,
  production_qc_file,
  canary_reproduction_file,
  primary_naive_reproduction_file
)

missing_combined_files <- required_combined_files[
  !file.exists(required_combined_files)
]

if (length(missing_combined_files) > 0L) {
  stop(
    "Required combined production files are missing:\n- ",
    paste(missing_combined_files, collapse = "\n- "),
    "\nRun 32_combine_phase2_mcar_results.R first."
  )
}

if (!any(grepl(
  "MCAR PRODUCTION COMBINATION: PASS",
  readLines(combine_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The Phase-2 MCAR production combination PASS marker is invalid.")
}

production_qc <- read_csv(
  production_qc_file,
  show_col_types = FALSE
)

canary_reproduction_qc <- read_csv(
  canary_reproduction_file,
  show_col_types = FALSE
)

primary_naive_reproduction_qc <- read_csv(
  primary_naive_reproduction_file,
  show_col_types = FALSE
)

if (
    !mcar_all_true(production_qc$passed) ||
    !mcar_all_true(production_qc$overall_production_content_pass) ||
    !mcar_all_true(canary_reproduction_qc$reproduced_within_tolerance) ||
    !mcar_all_true(primary_naive_reproduction_qc$reproduced_within_tolerance)
) {
  stop("Combined MCAR production outputs do not record a full pass.")
}


# ------------------------------------------------------------------------------
# 2. Scenario-level rsimsum and diagnostic summaries
# ------------------------------------------------------------------------------

stats_to_keep <- c(
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

expected_methods <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

rsimsum_tidy_list <- vector(
  "list",
  nrow(mcar_sensitivity_grid)
)

type1_list <- list()
type1_index <- 1L

failure_list <- vector(
  "list",
  nrow(mcar_sensitivity_grid)
)

runtime_list <- vector(
  "list",
  nrow(mcar_sensitivity_grid)
)

phase2_list <- vector(
  "list",
  nrow(mcar_sensitivity_grid)
)

weight_list <- vector(
  "list",
  nrow(mcar_sensitivity_grid)
)

scenario_results_list <- vector(
  "list",
  nrow(mcar_sensitivity_grid)
)

for (s in seq_len(nrow(mcar_sensitivity_grid))) {

  scenario_meta <- mcar_sensitivity_grid[s, , drop = FALSE]
  sensitivity_scenario_id <- scenario_meta$sensitivity_scenario_id
  theta_true <- scenario_meta$theta[[1]]

  cat(
    "Summarising MCAR scenario ",
    sensitivity_scenario_id,
    " / ",
    nrow(mcar_sensitivity_grid),
    " with rsimsum\n",
    sep = ""
  )

  scenario_file <- file.path(
    mcar_scenario_dir,
    sprintf(
      "mcar_scenario_%02d_estimates.rds",
      sensitivity_scenario_id
    )
  )

  if (!file.exists(scenario_file)) {
    stop(
      "MCAR scenario estimate file missing:\n",
      scenario_file
    )
  }

  dat <- readRDS(scenario_file)
  scenario_results_list[[s]] <- dat

  required_result_columns <- c(
    "repetition",
    "method",
    "status",
    "estimate",
    "se",
    "conf_low",
    "conf_high",
    "p_value",
    "elapsed_seconds"
  )

  missing_result_columns <- setdiff(
    required_result_columns,
    names(dat)
  )

  if (length(missing_result_columns) > 0L) {
    stop(
      "Scenario result file is missing columns: ",
      paste(missing_result_columns, collapse = ", ")
    )
  }

  result_keys <- dat %>%
    count(repetition, method, name = "n_rows")

  if (
    nrow(dat) != nsim_primary * 6L ||
      !setequal(unique(dat$repetition), seq_len(nsim_primary)) ||
      nrow(result_keys) != nsim_primary * 6L ||
      !mcar_all_true(result_keys$n_rows == 1L)
  ) {
    stop(
      "Scenario ",
      sensitivity_scenario_id,
      " does not contain the expected repetition-method keys."
    )
  }

  analysis_dat <- dat %>%
    filter(
      !is.na(status),
      status == "ok",
      is.finite(estimate),
      is.finite(se),
      se > 0,
      is.finite(conf_low),
      is.finite(conf_high)
    )

  if (!setequal(unique(analysis_dat$method), expected_methods)) {
    stop(
      "At least one method has no valid estimates in scenario ",
      sensitivity_scenario_id,
      "."
    )
  }

  rs_obj <- rsimsum::simsum(
    data = analysis_dat,
    estvarname = "estimate",
    true = theta_true,
    se = "se",
    methodvar = "method",
    ref = "CCA",
    ci.limits = c(
      "conf_low",
      "conf_high"
    ),
    dropbig = FALSE,
    x = FALSE,
    control = list(
      mcse = TRUE,
      level = nominal_level,
      na.rm = TRUE
    )
  )

  saveRDS(
    rs_obj,
    file.path(
      mcar_rsimsum_dir,
      sprintf(
        "mcar_scenario_%02d_rsimsum.rds",
        sensitivity_scenario_id
      )
    )
  )

  rs_summary <- summary(
    rs_obj,
    stats = stats_to_keep
  )

  rs_tidy <- generics::tidy(rs_summary) %>%
    mutate(
      sensitivity_scenario_id = sensitivity_scenario_id,
      primary_scenario_id = scenario_meta$primary_scenario_id[[1]],
      scenario_label = scenario_meta$scenario_label[[1]],
      scenario_key = scenario_meta$scenario_key[[1]],
      N = scenario_meta$N[[1]],
      phase2_fraction = scenario_meta$phase2_fraction[[1]],
      r2_a_marker = scenario_meta$r2_a_marker[[1]],
      r2_y_marker = scenario_meta$r2_y_marker[[1]],
      theta_true = theta_true,
      selection = "mcar",
      marker_error = "mvn"
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
      selection,
      marker_error
    )

  rsimsum_tidy_list[[s]] <- rs_tidy

  if (isTRUE(all.equal(theta_true, 0))) {

    type1_list[[type1_index]] <- dat %>%
      group_by(method) %>%
      summarise(
        nsim_total = n(),
        nsim_success = sum(
          status == "ok" & is.finite(p_value)
        ),
        null_rejection_rate = if (
          nsim_success > 0L
        ) {
          mean(
            p_value[
              status == "ok" & is.finite(p_value)
            ] < alpha_level
          )
        } else {
          NA_real_
        },
        mcse_null_rejection = if (
          nsim_success > 0L && is.finite(null_rejection_rate)
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
      ) %>%
      mutate(
        sensitivity_scenario_id = sensitivity_scenario_id,
        primary_scenario_id = scenario_meta$primary_scenario_id[[1]],
        scenario_label = scenario_meta$scenario_label[[1]],
        scenario_key = scenario_meta$scenario_key[[1]],
        N = scenario_meta$N[[1]],
        phase2_fraction = scenario_meta$phase2_fraction[[1]],
        r2_a_marker = scenario_meta$r2_a_marker[[1]],
        r2_y_marker = scenario_meta$r2_y_marker[[1]],
        theta_true = theta_true,
        selection = "mcar",
        marker_error = "mvn"
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
        selection,
        marker_error
      )

    type1_index <- type1_index + 1L
  }

  failure_list[[s]] <- dat %>%
    group_by(method) %>%
    summarise(
      nsim_total = n(),
      n_failed = sum(
        is.na(status) |
          status != "ok" |
          !is.finite(estimate) |
          !is.finite(se) |
          se <= 0
      ),
      failure_rate = n_failed / nsim_total,
      mcse_failure = sqrt(
        failure_rate * (1 - failure_rate) / nsim_total
      ),
      .groups = "drop"
    ) %>%
    mutate(
      sensitivity_scenario_id = sensitivity_scenario_id,
      primary_scenario_id = scenario_meta$primary_scenario_id[[1]],
      scenario_label = scenario_meta$scenario_label[[1]],
      scenario_key = scenario_meta$scenario_key[[1]],
      N = scenario_meta$N[[1]],
      phase2_fraction = scenario_meta$phase2_fraction[[1]],
      r2_a_marker = scenario_meta$r2_a_marker[[1]],
      r2_y_marker = scenario_meta$r2_y_marker[[1]],
      theta_true = theta_true,
      selection = "mcar",
      marker_error = "mvn"
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
      selection,
      marker_error
    )

  runtime_list[[s]] <- dat %>%
    group_by(method) %>%
    summarise(
      mean_runtime_seconds = mean(elapsed_seconds, na.rm = TRUE),
      median_runtime_seconds = median(elapsed_seconds, na.rm = TRUE),
      p95_runtime_seconds = unname(
        quantile(elapsed_seconds, 0.95, na.rm = TRUE)
      ),
      total_cpu_hours = sum(elapsed_seconds, na.rm = TRUE) / 3600,
      .groups = "drop"
    ) %>%
    mutate(
      sensitivity_scenario_id = sensitivity_scenario_id,
      primary_scenario_id = scenario_meta$primary_scenario_id[[1]],
      scenario_label = scenario_meta$scenario_label[[1]],
      scenario_key = scenario_meta$scenario_key[[1]],
      N = scenario_meta$N[[1]],
      phase2_fraction = scenario_meta$phase2_fraction[[1]],
      r2_a_marker = scenario_meta$r2_a_marker[[1]],
      r2_y_marker = scenario_meta$r2_y_marker[[1]],
      theta_true = theta_true,
      selection = "mcar",
      marker_error = "mvn"
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
      selection,
      marker_error
    )

  phase2_list[[s]] <- dat %>%
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
      mean_n_phase2 = mean(n_phase2),
      sd_n_phase2 = sd(n_phase2),
      min_n_phase2 = min(n_phase2),
      max_n_phase2 = max(n_phase2),
      mean_realized_phase2_fraction =
        mean(realized_phase2_fraction),
      mean_min_pi_true = mean(min_pi_true),
      mean_p01_pi_true = mean(p01_pi_true),
      mean_median_pi_true = mean(median_pi_true),
      mean_p99_pi_true = mean(p99_pi_true),
      mean_max_pi_true = mean(max_pi_true)
    ) %>%
    mutate(
      sensitivity_scenario_id = sensitivity_scenario_id,
      primary_scenario_id = scenario_meta$primary_scenario_id[[1]],
      scenario_label = scenario_meta$scenario_label[[1]],
      scenario_key = scenario_meta$scenario_key[[1]],
      N = scenario_meta$N[[1]],
      phase2_fraction = scenario_meta$phase2_fraction[[1]],
      r2_a_marker = scenario_meta$r2_a_marker[[1]],
      r2_y_marker = scenario_meta$r2_y_marker[[1]],
      theta_true = theta_true,
      selection = "mcar",
      marker_error = "mvn"
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
      selection,
      marker_error
    )

  weight_list[[s]] <- dat %>%
    filter(
      method %in% c("IPW", "AIPW"),
      status == "ok"
    ) %>%
    group_by(method) %>%
    summarise(
      mean_weight_min = mean(weight_min, na.rm = TRUE),
      mean_weight_p99 = mean(weight_p99, na.rm = TRUE),
      mean_weight_max = mean(weight_max, na.rm = TRUE),
      mean_weight_cv = mean(weight_cv, na.rm = TRUE),
      p95_weight_cv = unname(
        quantile(weight_cv, 0.95, na.rm = TRUE)
      ),
      mean_weight_ess = mean(weight_ess, na.rm = TRUE),
      min_weight_ess = min(weight_ess, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      sensitivity_scenario_id = sensitivity_scenario_id,
      primary_scenario_id = scenario_meta$primary_scenario_id[[1]],
      scenario_label = scenario_meta$scenario_label[[1]],
      scenario_key = scenario_meta$scenario_key[[1]],
      N = scenario_meta$N[[1]],
      phase2_fraction = scenario_meta$phase2_fraction[[1]],
      r2_a_marker = scenario_meta$r2_a_marker[[1]],
      r2_y_marker = scenario_meta$r2_y_marker[[1]],
      theta_true = theta_true,
      selection = "mcar",
      marker_error = "mvn"
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
      selection,
      marker_error
    )
}


# ------------------------------------------------------------------------------
# 3. Save MCAR performance and diagnostic summaries
# ------------------------------------------------------------------------------

rsimsum_tidy <- bind_rows(rsimsum_tidy_list)

write_csv(
  rsimsum_tidy,
  file.path(
    mcar_summary_dir,
    "phase2_mcar_rsimsum_tidy.csv"
  )
)

saveRDS(
  rsimsum_tidy,
  file.path(
    mcar_summary_dir,
    "phase2_mcar_rsimsum_tidy.rds"
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
    selection,
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

if ("est__relerror" %in% names(performance_wide)) {
  performance_wide <- performance_wide %>%
    mutate(
      se_ratio = 1 + est__relerror / 100,
      mcse_se_ratio = mcse__relerror / 100
    )
}

type1_summary <- bind_rows(type1_list)
failure_summary <- bind_rows(failure_list)
runtime_summary <- bind_rows(runtime_list)
phase2_summary <- bind_rows(phase2_list)
weight_summary <- bind_rows(weight_list)

write_csv(
  performance_wide,
  file.path(
    mcar_summary_dir,
    "phase2_mcar_performance_wide.csv"
  )
)

write_csv(
  type1_summary,
  file.path(
    mcar_summary_dir,
    "phase2_mcar_type1_error.csv"
  )
)

write_csv(
  failure_summary,
  file.path(
    mcar_summary_dir,
    "phase2_mcar_failure_rates.csv"
  )
)

write_csv(
  runtime_summary,
  file.path(
    mcar_summary_dir,
    "phase2_mcar_runtime.csv"
  )
)

write_csv(
  phase2_summary,
  file.path(
    mcar_summary_dir,
    "phase2_mcar_phase2_diagnostics.csv"
  )
)

write_csv(
  weight_summary,
  file.path(
    mcar_summary_dir,
    "phase2_mcar_weight_diagnostics.csv"
  )
)


# ------------------------------------------------------------------------------
# 4. Compare MCAR performance with matched primary MAR DGMs
# ------------------------------------------------------------------------------

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
    paste(missing_primary_files, collapse = "\n- ")
  )
}

primary_performance <- read_csv(
  primary_performance_file,
  show_col_types = FALSE
) %>%
  filter(
    scenario_id %in%
      mcar_sensitivity_grid$primary_scenario_id
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
    "Required performance metrics are missing from primary or MCAR summaries: ",
    paste(missing_metrics, collapse = ", ")
  )
}

mcar_comparison_input <- performance_wide %>%
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
    ~ paste0("mcar_", .x),
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

performance_comparison <- mcar_comparison_input %>%
  inner_join(
    primary_comparison_input,
    by = c(
      "primary_scenario_id",
      "method"
    )
  ) %>%
  mutate(
    delta_mean_estimate =
      mcar_est__thetamean - primary_est__thetamean,
    delta_bias = mcar_est__bias - primary_est__bias,
    delta_absolute_bias =
      abs(mcar_est__bias) - abs(primary_est__bias),
    delta_empirical_se =
      mcar_est__empse - primary_est__empse,
    delta_model_se =
      mcar_est__modelse - primary_est__modelse,
    delta_se_ratio =
      mcar_se_ratio - primary_se_ratio,
    delta_coverage =
      mcar_est__cover - primary_est__cover,
    delta_bias_eliminated_coverage =
      mcar_est__becover - primary_est__becover,
    delta_mse = mcar_est__mse - primary_est__mse
  )

write_csv(
  performance_comparison,
  file.path(
    mcar_summary_dir,
    "mcar_vs_primary_performance_comparison.csv"
  )
)

primary_type1 <- read_csv(
  primary_type1_file,
  show_col_types = FALSE
) %>%
  filter(
    scenario_id %in%
      mcar_sensitivity_grid$primary_scenario_id
  ) %>%
  select(
    scenario_id,
    method,
    primary_null_rejection_rate = null_rejection_rate,
    primary_mcse_null_rejection = mcse_null_rejection
  ) %>%
  rename(primary_scenario_id = scenario_id)

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
    mcar_null_rejection_rate = null_rejection_rate,
    mcar_mcse_null_rejection = mcse_null_rejection
  ) %>%
  inner_join(
    primary_type1,
    by = c(
      "primary_scenario_id",
      "method"
    )
  ) %>%
  mutate(
    delta_null_rejection_rate =
      mcar_null_rejection_rate -
      primary_null_rejection_rate
  )

write_csv(
  type1_comparison,
  file.path(
    mcar_summary_dir,
    "mcar_vs_primary_type1_comparison.csv"
  )
)

primary_failures <- read_csv(
  primary_failure_file,
  show_col_types = FALSE
) %>%
  filter(
    scenario_id %in%
      mcar_sensitivity_grid$primary_scenario_id
  ) %>%
  select(
    scenario_id,
    method,
    primary_n_failed = n_failed,
    primary_failure_rate = failure_rate,
    primary_mcse_failure = mcse_failure
  ) %>%
  rename(primary_scenario_id = scenario_id)

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
    mcar_n_failed = n_failed,
    mcar_failure_rate = failure_rate,
    mcar_mcse_failure = mcse_failure
  ) %>%
  inner_join(
    primary_failures,
    by = c(
      "primary_scenario_id",
      "method"
    )
  ) %>%
  mutate(
    delta_failure_rate =
      mcar_failure_rate - primary_failure_rate
  )

write_csv(
  failure_comparison,
  file.path(
    mcar_summary_dir,
    "mcar_vs_primary_failure_comparison.csv"
  )
)


# ------------------------------------------------------------------------------
# 5. Paired repetition-level comparison with the primary MAR simulation
# ------------------------------------------------------------------------------

paired_summary_list <- vector(
  "list",
  nrow(mcar_sensitivity_grid)
)

paired_repetition_list <- vector(
  "list",
  nrow(mcar_sensitivity_grid)
)

for (s in seq_len(nrow(mcar_sensitivity_grid))) {

  scenario_meta <- mcar_sensitivity_grid[s, , drop = FALSE]
  # Use a distinct scalar name so dplyr's theta_true column cannot shadow the
  # scenario-level truth in later mutate() and summarise() expressions.
  scenario_theta <- scenario_meta$theta[[1]]

  primary_scenario_file <- file.path(
    combined_dir,
    "scenarios",
    sprintf(
      "scenario_%03d_estimates.rds",
      scenario_meta$primary_scenario_id
    )
  )

  if (!file.exists(primary_scenario_file)) {
    stop("Matched primary scenario file missing: ", primary_scenario_file)
  }

  primary_repetitions <- readRDS(primary_scenario_file) %>%
    select(
      repetition,
      method,
      primary_status = status,
      primary_estimate = estimate,
      primary_se = se,
      primary_conf_low = conf_low,
      primary_conf_high = conf_high,
      primary_p_value = p_value
    )

  mcar_repetitions <- scenario_results_list[[s]] %>%
    select(
      repetition,
      method,
      mcar_status = status,
      mcar_estimate = estimate,
      mcar_se = se,
      mcar_conf_low = conf_low,
      mcar_conf_high = conf_high,
      mcar_p_value = p_value
    )

  paired_repetitions <- mcar_repetitions %>%
    inner_join(
      primary_repetitions,
      by = c("repetition", "method")
    ) %>%
    mutate(
      sensitivity_scenario_id =
        scenario_meta$sensitivity_scenario_id[[1]],
      primary_scenario_id = scenario_meta$primary_scenario_id[[1]],
      scenario_label = scenario_meta$scenario_label[[1]],
      theta_true = scenario_theta,
      pair_success =
        mcar_status == "ok" &
        primary_status == "ok" &
        is.finite(mcar_estimate) &
        is.finite(primary_estimate) &
        is.finite(mcar_se) &
        is.finite(primary_se),
      delta_estimate = mcar_estimate - primary_estimate,
      delta_model_se = mcar_se - primary_se,
      delta_absolute_error =
        abs(mcar_estimate - scenario_theta) -
        abs(primary_estimate - scenario_theta),
      delta_squared_error =
        (mcar_estimate - scenario_theta)^2 -
        (primary_estimate - scenario_theta)^2,
      mcar_covered =
        mcar_conf_low <= scenario_theta &
        mcar_conf_high >= scenario_theta,
      primary_covered =
        primary_conf_low <= scenario_theta &
        primary_conf_high >= scenario_theta,
      delta_coverage_indicator =
        as.numeric(mcar_covered) - as.numeric(primary_covered),
      delta_rejection_indicator = if (scenario_theta == 0) {
        as.numeric(mcar_p_value < alpha_level) -
          as.numeric(primary_p_value < alpha_level)
      } else {
        NA_real_
      }
    ) %>%
    relocate(
      sensitivity_scenario_id,
      primary_scenario_id,
      scenario_label,
      theta_true,
      repetition,
      method
    )

  paired_repetition_list[[s]] <- paired_repetitions

  paired_summary_list[[s]] <- paired_repetitions %>%
    filter(pair_success) %>%
    group_by(method) %>%
    summarise(
      n_pairs = n(),
      mean_delta_estimate = mean(delta_estimate),
      mcse_mean_delta_estimate =
        sd(delta_estimate) / sqrt(dplyr::n()),
      mean_delta_model_se = mean(delta_model_se),
      mean_delta_absolute_error = mean(delta_absolute_error),
      mean_delta_squared_error = mean(delta_squared_error),
      mean_delta_coverage = mean(delta_coverage_indicator),
      mean_delta_null_rejection = if (scenario_theta == 0) {
        mean(delta_rejection_indicator)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(
      sensitivity_scenario_id =
        scenario_meta$sensitivity_scenario_id[[1]],
      primary_scenario_id = scenario_meta$primary_scenario_id[[1]],
      scenario_label = scenario_meta$scenario_label[[1]],
      scenario_key = scenario_meta$scenario_key[[1]],
      N = scenario_meta$N[[1]],
      phase2_fraction = scenario_meta$phase2_fraction[[1]],
      r2_a_marker = scenario_meta$r2_a_marker[[1]],
      r2_y_marker = scenario_meta$r2_y_marker[[1]],
      theta_true = scenario_theta
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
      method
    )
}

paired_repetition_differences <- bind_rows(paired_repetition_list)
paired_difference_summary <- bind_rows(paired_summary_list)

saveRDS(
  paired_repetition_differences,
  file.path(
    mcar_summary_dir,
    "mcar_vs_primary_paired_repetition_differences.rds"
  ),
  compress = "xz"
)

write_csv(
  paired_difference_summary,
  file.path(
    mcar_summary_dir,
    "mcar_vs_primary_paired_difference_summary.csv"
  )
)


# ------------------------------------------------------------------------------
# 6. Final summary QC and compact archive
# ------------------------------------------------------------------------------

expected_method_scenario_rows <-
  nrow(mcar_sensitivity_grid) * 6L

expected_null_method_rows <- 2L * 6L

key_performance_columns <- c(
  "est__bias",
  "est__empse",
  "est__modelse",
  "se_ratio",
  "est__cover",
  "est__becover",
  "est__mse"
)

missing_key_performance <- setdiff(
  key_performance_columns,
  names(performance_wide)
)

if (length(missing_key_performance) > 0L) {
  stop(
    "Required performance summaries are missing: ",
    paste(missing_key_performance, collapse = ", ")
  )
}

finite_key_performance <- mcar_all_true(
  vapply(
    performance_wide[key_performance_columns],
    function(x) mcar_all_true(is.finite(x)),
    logical(1)
  )
)

summary_qc <- data.frame(
  check = c(
    "Four MCAR scenarios were summarized",
    "All 24 scenario-method performance rows are present",
    "All 12 null-scenario method rows are present",
    "All key performance summaries are finite",
    "Matched primary comparison contains 24 rows",
    "Matched Type I error comparison contains 12 rows",
    "Matched failure comparison contains 24 rows",
    "Paired primary comparison contains 24 rows",
    "Canary and paired Naive reproduction remain verified"
  ),
  hard_check = rep(TRUE, 9L),
  passed = c(
    length(scenario_results_list) == 4L,
    nrow(performance_wide) == expected_method_scenario_rows,
    nrow(type1_summary) == expected_null_method_rows,
    finite_key_performance,
    nrow(performance_comparison) == expected_method_scenario_rows,
    nrow(type1_comparison) == expected_null_method_rows,
    nrow(failure_comparison) == expected_method_scenario_rows,
    nrow(paired_difference_summary) == expected_method_scenario_rows,
    mcar_all_true(canary_reproduction_qc$reproduced_within_tolerance) &
      mcar_all_true(primary_naive_reproduction_qc$reproduced_within_tolerance)
  ),
  stringsAsFactors = FALSE
)

summary_qc$overall_summary_pass <- mcar_all_true(
  summary_qc$passed[summary_qc$hard_check]
)

write_csv(
  summary_qc,
  file.path(
    mcar_summary_dir,
    "phase2_mcar_summary_qc.csv"
  )
)

review_flags <- failure_summary %>%
  mutate(
    exceeds_prespecified_failure_threshold =
      failure_rate > deterioration_thresholds$failure_rate
  ) %>%
  filter(exceeds_prespecified_failure_threshold)

write_csv(
  review_flags,
  file.path(
    mcar_summary_dir,
    "phase2_mcar_failure_review_flags.csv"
  )
)

saveRDS(
  list(
    scenarios = mcar_sensitivity_grid,
    rsimsum_tidy = rsimsum_tidy,
    performance_wide = performance_wide,
    type1_error = type1_summary,
    failure_rates = failure_summary,
    runtime = runtime_summary,
    phase2_diagnostics = phase2_summary,
    weight_diagnostics = weight_summary,
    performance_comparison = performance_comparison,
    type1_comparison = type1_comparison,
    failure_comparison = failure_comparison,
    paired_difference_summary = paired_difference_summary,
    summary_qc = summary_qc,
    failure_review_flags = review_flags
  ),
  file.path(
    mcar_summary_dir,
    "phase2_mcar_production_summary_complete.rds"
  ),
  compress = "xz"
)

if (!mcar_all_true(summary_qc$passed[summary_qc$hard_check])) {

  failed_checks <- summary_qc$check[
    summary_qc$hard_check & !summary_qc$passed
  ]

  stop(
    "Phase-2 MCAR production summary FAILED. Failed checks:\n- ",
    paste(failed_checks, collapse = "\n- "),
    "\nSee outputs in: ",
    mcar_summary_dir
  )
}

mcar_atomic_write_lines(
  c(
    "PHASE-2 MCAR PRODUCTION SUMMARY: PASS",
    paste("Completed:", Sys.time()),
    paste("Production scenarios:", nrow(mcar_sensitivity_grid)),
    paste("Scenario-method performance rows:", nrow(performance_wide)),
    paste("Null-scenario method rows:", nrow(type1_summary)),
    paste("Failure-rate review flags:", nrow(review_flags))
  ),
  summary_pass_file
)

cat(
  "\nPhase-2 MCAR production summary PASSED.\n",
  "Main comparison file:\n",
  file.path(
    mcar_summary_dir,
    "mcar_vs_primary_performance_comparison.csv"
  ),
  "\nSummary directory: ",
  mcar_summary_dir,
  "\n",
  sep = ""
)
