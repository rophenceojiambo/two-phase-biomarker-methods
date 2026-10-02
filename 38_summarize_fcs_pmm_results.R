################################################################################
# 38_summarize_fcs_pmm_results.R
#
# ADEMP summary of FCS-PMM in the four targeted primary MAR/MVN scenarios.
# rsimsum is the authoritative summary engine. Aggregate and paired comparisons
# are made against the completed primary FCS-MI analysis using method = "norm".
################################################################################

source("00_config.R")
source("34_fcs_pmm_sensitivity_helpers.R")

library(dplyr)
library(tidyr)
library(readr)
library(rsimsum)


# ------------------------------------------------------------------------------
# 1. Paths and completed-input gate
# ------------------------------------------------------------------------------

fcs_pmm_results_dir <- file.path(results_dir, "fcs_pmm_sensitivity")
fcs_pmm_production_dir <- file.path(fcs_pmm_results_dir, "production")
fcs_pmm_combined_dir <- file.path(fcs_pmm_production_dir, "combined")
fcs_pmm_scenario_dir <- file.path(fcs_pmm_combined_dir, "scenarios")
fcs_pmm_summary_dir <- file.path(fcs_pmm_production_dir, "summary")
fcs_pmm_rsimsum_dir <- file.path(fcs_pmm_summary_dir, "rsimsum_objects")

dir.create(fcs_pmm_summary_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fcs_pmm_rsimsum_dir, showWarnings = FALSE, recursive = TRUE)

summary_pass_file <- file.path(
  fcs_pmm_summary_dir,
  "FCS_PMM_PRODUCTION_SUMMARY_PASS.txt"
)

if (file.exists(summary_pass_file)) {
  unlink(summary_pass_file)
}

combine_pass_file <- file.path(
  fcs_pmm_combined_dir,
  "FCS_PMM_PRODUCTION_COMBINE_PASS.txt"
)

production_qc_file <- file.path(
  fcs_pmm_combined_dir,
  "fcs_pmm_production_qc.csv"
)

canary_reproduction_file <- file.path(
  fcs_pmm_combined_dir,
  "fcs_pmm_canary_reproduction_qc.csv"
)

norm_reproduction_file <- file.path(
  fcs_pmm_results_dir,
  "canary",
  "fcs_norm_primary_reproduction_qc.csv"
)

required_files <- c(
  combine_pass_file,
  production_qc_file,
  canary_reproduction_file,
  norm_reproduction_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required FCS-PMM combined files are missing:\n- ",
    paste(missing_files, collapse = "\n- "),
    "\nRun 37_combine_fcs_pmm_results.R first."
  )
}

if (!any(grepl(
  "FCS-PMM PRODUCTION COMBINATION: PASS",
  readLines(combine_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The FCS-PMM production combination PASS marker is invalid.")
}

production_qc <- read_csv(production_qc_file, show_col_types = FALSE)
canary_reproduction_qc <- read_csv(
  canary_reproduction_file,
  show_col_types = FALSE
)
norm_reproduction_qc <- read_csv(
  norm_reproduction_file,
  show_col_types = FALSE
)

if (
  !fcs_pmm_all_true(production_qc$passed) ||
    !fcs_pmm_all_true(production_qc$overall_production_content_pass) ||
    !fcs_pmm_all_true(canary_reproduction_qc$reproduced_within_tolerance) ||
    !fcs_pmm_all_true(norm_reproduction_qc$reproduced_within_tolerance)
) {
  stop("Combined FCS-PMM outputs do not record a full pass.")
}


# ------------------------------------------------------------------------------
# 2. Scenario-level rsimsum and diagnostics
# ------------------------------------------------------------------------------

stats_to_keep <- c(
  "nsim",
  "thetamean",
  "bias",
  "empse",
  "mse",
  "modelse",
  "relerror",
  "cover",
  "becover",
  "power"
)

rsimsum_tidy_list <- vector("list", nrow(fcs_pmm_sensitivity_grid))
type1_list <- list()
type1_index <- 1L
failure_list <- vector("list", nrow(fcs_pmm_sensitivity_grid))
runtime_list <- vector("list", nrow(fcs_pmm_sensitivity_grid))
phase2_list <- vector("list", nrow(fcs_pmm_sensitivity_grid))
scenario_results_list <- vector("list", nrow(fcs_pmm_sensitivity_grid))

for (s in seq_len(nrow(fcs_pmm_sensitivity_grid))) {

  scenario_meta <- fcs_pmm_sensitivity_grid[s, , drop = FALSE]
  sensitivity_scenario_id <- scenario_meta$sensitivity_scenario_id
  theta_true <- scenario_meta$theta[[1]]

  cat(
    "Summarising FCS-PMM scenario ",
    sensitivity_scenario_id,
    " / ",
    nrow(fcs_pmm_sensitivity_grid),
    " with rsimsum\n",
    sep = ""
  )

  scenario_file <- file.path(
    fcs_pmm_scenario_dir,
    sprintf(
      "fcs_pmm_scenario_%02d_estimates.rds",
      sensitivity_scenario_id
    )
  )

  if (!file.exists(scenario_file)) {
    stop("FCS-PMM scenario estimate file missing:\n", scenario_file)
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
    nrow(dat) != nsim_primary ||
      !setequal(unique(dat$repetition), seq_len(nsim_primary)) ||
      nrow(result_keys) != nsim_primary ||
      !fcs_pmm_all_true(result_keys$n_rows == 1L) ||
      !identical(unique(dat$method), fcs_pmm_label)
  ) {
    stop(
      "Scenario ",
      sensitivity_scenario_id,
      " does not contain the expected FCS-PMM repetition keys."
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

  if (nrow(analysis_dat) == 0L) {
    stop(
      "Scenario ",
      sensitivity_scenario_id,
      " contains no valid FCS-PMM estimates."
    )
  }

  rs_obj <- rsimsum::simsum(
    data = analysis_dat,
    estvarname = "estimate",
    true = theta_true,
    se = "se",
    methodvar = "method",
    ci.limits = c("conf_low", "conf_high"),
    dropbig = FALSE,
    x = FALSE,
    control = list(
      mcse = TRUE,
      level = nominal_level,
      na.rm = TRUE
    )
  )

  fcs_pmm_atomic_save_rds(
    rs_obj,
    file.path(
      fcs_pmm_rsimsum_dir,
      sprintf("fcs_pmm_scenario_%02d_rsimsum.rds", sensitivity_scenario_id)
    ),
    compress = TRUE
  )

  rs_tidy <- generics::tidy(
    summary(rs_obj, stats = stats_to_keep)
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
      selection = "mar",
      marker_error = "mvn",
      mice_method = "pmm"
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
      marker_error,
      mice_method
    )

  rsimsum_tidy_list[[s]] <- rs_tidy

  if (isTRUE(all.equal(theta_true, 0))) {

    type1_list[[type1_index]] <- dat %>%
      summarise(
        method = first(method),
        nsim_total = n(),
        nsim_success = sum(
          !is.na(status) & status == "ok" & is.finite(p_value)
        ),
        null_rejection_rate = if (nsim_success > 0L) {
          mean(
            p_value[
              !is.na(status) & status == "ok" & is.finite(p_value)
            ] < alpha_level
          )
        } else {
          NA_real_
        },
        mcse_null_rejection = if (
          nsim_success > 0L && is.finite(null_rejection_rate)
        ) {
          sqrt(
            null_rejection_rate * (1 - null_rejection_rate) /
              nsim_success
          )
        } else {
          NA_real_
        }
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
        theta_true = theta_true
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

    type1_index <- type1_index + 1L
  }

  failure_list[[s]] <- dat %>%
    summarise(
      method = first(method),
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
      )
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
      theta_true = theta_true
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

  runtime_list[[s]] <- dat %>%
    summarise(
      method = first(method),
      mean_runtime_seconds = mean(elapsed_seconds, na.rm = TRUE),
      median_runtime_seconds = median(elapsed_seconds, na.rm = TRUE),
      p95_runtime_seconds = unname(
        quantile(elapsed_seconds, 0.95, na.rm = TRUE)
      ),
      total_cpu_hours = sum(elapsed_seconds, na.rm = TRUE) / 3600
    ) %>%
    mutate(
      sensitivity_scenario_id = sensitivity_scenario_id,
      primary_scenario_id = scenario_meta$primary_scenario_id[[1]],
      scenario_label = scenario_meta$scenario_label[[1]],
      scenario_key = scenario_meta$scenario_key[[1]],
      N = scenario_meta$N[[1]],
      phase2_fraction = scenario_meta$phase2_fraction[[1]],
      theta_true = theta_true
    ) %>%
    relocate(
      sensitivity_scenario_id,
      primary_scenario_id,
      scenario_label,
      scenario_key,
      N,
      phase2_fraction,
      theta_true,
      method
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
      mean_realized_phase2_fraction = mean(realized_phase2_fraction),
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
      theta_true = theta_true
    ) %>%
    relocate(
      sensitivity_scenario_id,
      primary_scenario_id,
      scenario_label,
      scenario_key,
      N,
      phase2_fraction,
      theta_true
    )
}


# ------------------------------------------------------------------------------
# 3. Save FCS-PMM performance and diagnostic summaries
# ------------------------------------------------------------------------------

rsimsum_tidy <- bind_rows(rsimsum_tidy_list)

fcs_pmm_atomic_write_csv(
  rsimsum_tidy,
  file.path(fcs_pmm_summary_dir, "fcs_pmm_rsimsum_tidy.csv")
)

fcs_pmm_atomic_save_rds(
  rsimsum_tidy,
  file.path(fcs_pmm_summary_dir, "fcs_pmm_rsimsum_tidy.rds"),
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
    mice_method,
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

fcs_pmm_atomic_write_csv(
  performance_wide,
  file.path(fcs_pmm_summary_dir, "fcs_pmm_performance_wide.csv")
)

fcs_pmm_atomic_write_csv(
  type1_summary,
  file.path(fcs_pmm_summary_dir, "fcs_pmm_type1_error.csv")
)

fcs_pmm_atomic_write_csv(
  failure_summary,
  file.path(fcs_pmm_summary_dir, "fcs_pmm_failure_rates.csv")
)

fcs_pmm_atomic_write_csv(
  runtime_summary,
  file.path(fcs_pmm_summary_dir, "fcs_pmm_runtime.csv")
)

fcs_pmm_atomic_write_csv(
  phase2_summary,
  file.path(fcs_pmm_summary_dir, "fcs_pmm_phase2_diagnostics.csv")
)


# ------------------------------------------------------------------------------
# 4. Aggregate comparison with primary FCS-Norm
# ------------------------------------------------------------------------------

primary_performance_file <- file.path(
  summary_dir,
  "primary_performance_wide.csv"
)

primary_type1_file <- file.path(summary_dir, "primary_type1_error.csv")
primary_failure_file <- file.path(summary_dir, "primary_failure_rates.csv")

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

comparison_metrics <- c(
  "est__thetamean",
  "est__bias",
  "est__empse",
  "est__modelse",
  "se_ratio",
  "est__cover",
  "est__becover",
  "est__mse"
)

primary_performance <- read_csv(
  primary_performance_file,
  show_col_types = FALSE
) %>%
  filter(
    scenario_id %in% fcs_pmm_sensitivity_grid$primary_scenario_id,
    method == "FCS-MI"
  )

if (
  nrow(primary_performance) != nrow(fcs_pmm_sensitivity_grid) ||
    anyDuplicated(primary_performance$scenario_id)
) {
  stop("Primary FCS-Norm performance rows are missing or duplicated.")
}

if (
  nrow(performance_wide) != nrow(fcs_pmm_sensitivity_grid) ||
    anyDuplicated(performance_wide$primary_scenario_id)
) {
  stop("FCS-PMM performance rows are missing or duplicated.")
}

missing_metrics <- setdiff(
  comparison_metrics,
  intersect(names(primary_performance), names(performance_wide))
)

if (length(missing_metrics) > 0L) {
  stop(
    "Required comparison metrics are missing: ",
    paste(missing_metrics, collapse = ", ")
  )
}

pmm_comparison_input <- performance_wide %>%
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
    all_of(comparison_metrics)
  ) %>%
  rename_with(
    ~ paste0("pmm_", .x),
    all_of(comparison_metrics)
  )

norm_comparison_input <- primary_performance %>%
  select(scenario_id, all_of(comparison_metrics)) %>%
  rename(primary_scenario_id = scenario_id) %>%
  rename_with(
    ~ paste0("norm_", .x),
    all_of(comparison_metrics)
  )

performance_comparison <- pmm_comparison_input %>%
  inner_join(norm_comparison_input, by = "primary_scenario_id") %>%
  mutate(
    delta_mean_estimate = pmm_est__thetamean - norm_est__thetamean,
    delta_bias = pmm_est__bias - norm_est__bias,
    delta_absolute_bias = abs(pmm_est__bias) - abs(norm_est__bias),
    delta_empirical_se = pmm_est__empse - norm_est__empse,
    delta_model_se = pmm_est__modelse - norm_est__modelse,
    delta_se_ratio = pmm_se_ratio - norm_se_ratio,
    delta_coverage = pmm_est__cover - norm_est__cover,
    delta_bias_eliminated_coverage =
      pmm_est__becover - norm_est__becover,
    delta_mse = pmm_est__mse - norm_est__mse
  )

fcs_pmm_atomic_write_csv(
  performance_comparison,
  file.path(fcs_pmm_summary_dir, "fcs_pmm_vs_norm_performance.csv")
)

primary_type1 <- read_csv(primary_type1_file, show_col_types = FALSE) %>%
  filter(
    scenario_id %in% fcs_pmm_sensitivity_grid$primary_scenario_id,
    method == "FCS-MI"
  ) %>%
  select(
    scenario_id,
    norm_null_rejection_rate = null_rejection_rate,
    norm_mcse_null_rejection = mcse_null_rejection
  ) %>%
  rename(primary_scenario_id = scenario_id)

if (
  nrow(primary_type1) != 2L ||
    anyDuplicated(primary_type1$primary_scenario_id)
) {
  stop("Primary FCS-Norm Type I error rows are missing or duplicated.")
}

if (
  nrow(type1_summary) != 2L ||
    anyDuplicated(type1_summary$primary_scenario_id)
) {
  stop("FCS-PMM Type I error rows are missing or duplicated.")
}

type1_comparison <- type1_summary %>%
  select(
    sensitivity_scenario_id,
    primary_scenario_id,
    scenario_label,
    theta_true,
    pmm_null_rejection_rate = null_rejection_rate,
    pmm_mcse_null_rejection = mcse_null_rejection
  ) %>%
  inner_join(primary_type1, by = "primary_scenario_id") %>%
  mutate(
    delta_null_rejection_rate =
      pmm_null_rejection_rate - norm_null_rejection_rate
  )

fcs_pmm_atomic_write_csv(
  type1_comparison,
  file.path(fcs_pmm_summary_dir, "fcs_pmm_vs_norm_type1.csv")
)

primary_failures <- read_csv(
  primary_failure_file,
  show_col_types = FALSE
) %>%
  filter(
    scenario_id %in% fcs_pmm_sensitivity_grid$primary_scenario_id,
    method == "FCS-MI"
  ) %>%
  select(
    scenario_id,
    norm_n_failed = n_failed,
    norm_failure_rate = failure_rate,
    norm_mcse_failure = mcse_failure
  ) %>%
  rename(primary_scenario_id = scenario_id)

if (
  nrow(primary_failures) != nrow(fcs_pmm_sensitivity_grid) ||
    anyDuplicated(primary_failures$primary_scenario_id)
) {
  stop("Primary FCS-Norm failure rows are missing or duplicated.")
}

if (
  nrow(failure_summary) != nrow(fcs_pmm_sensitivity_grid) ||
    anyDuplicated(failure_summary$primary_scenario_id)
) {
  stop("FCS-PMM failure rows are missing or duplicated.")
}

failure_comparison <- failure_summary %>%
  select(
    sensitivity_scenario_id,
    primary_scenario_id,
    scenario_label,
    theta_true,
    pmm_n_failed = n_failed,
    pmm_failure_rate = failure_rate,
    pmm_mcse_failure = mcse_failure
  ) %>%
  inner_join(primary_failures, by = "primary_scenario_id") %>%
  mutate(delta_failure_rate = pmm_failure_rate - norm_failure_rate)

fcs_pmm_atomic_write_csv(
  failure_comparison,
  file.path(fcs_pmm_summary_dir, "fcs_pmm_vs_norm_failures.csv")
)


# ------------------------------------------------------------------------------
# 5. Paired repetition-level FCS-PMM versus primary FCS-Norm comparison
# ------------------------------------------------------------------------------

paired_repetition_list <- vector("list", nrow(fcs_pmm_sensitivity_grid))
paired_summary_list <- vector("list", nrow(fcs_pmm_sensitivity_grid))

for (s in seq_len(nrow(fcs_pmm_sensitivity_grid))) {

  scenario_meta <- fcs_pmm_sensitivity_grid[s, , drop = FALSE]
  # Keep the scalar scenario truth separate from the theta_true column created
  # inside dplyr verbs below, preventing data-mask shadowing in scalar if().
  scenario_theta <- scenario_meta$theta[[1]]

  primary_file <- file.path(
    combined_dir,
    "scenarios",
    sprintf(
      "scenario_%03d_estimates.rds",
      scenario_meta$primary_scenario_id
    )
  )

  if (!file.exists(primary_file)) {
    stop("Matched primary scenario file missing: ", primary_file)
  }

  norm_repetitions <- readRDS(primary_file) %>%
    filter(method == "FCS-MI") %>%
    select(
      repetition,
      norm_status = status,
      norm_estimate = estimate,
      norm_se = se,
      norm_conf_low = conf_low,
      norm_conf_high = conf_high,
      norm_p_value = p_value
    )

  pmm_repetitions <- scenario_results_list[[s]] %>%
    select(
      repetition,
      pmm_status = status,
      pmm_estimate = estimate,
      pmm_se = se,
      pmm_conf_low = conf_low,
      pmm_conf_high = conf_high,
      pmm_p_value = p_value
    )

  if (
    nrow(norm_repetitions) != nsim_primary ||
      nrow(pmm_repetitions) != nsim_primary ||
      anyDuplicated(norm_repetitions$repetition) ||
      anyDuplicated(pmm_repetitions$repetition) ||
      !setequal(norm_repetitions$repetition, seq_len(nsim_primary)) ||
      !setequal(pmm_repetitions$repetition, seq_len(nsim_primary))
  ) {
    stop(
      "Paired repetition keys are incomplete for primary scenario ",
      scenario_meta$primary_scenario_id,
      "."
    )
  }

  paired_repetitions <- pmm_repetitions %>%
    inner_join(norm_repetitions, by = "repetition") %>%
    mutate(
      sensitivity_scenario_id =
        scenario_meta$sensitivity_scenario_id[[1]],
      primary_scenario_id = scenario_meta$primary_scenario_id[[1]],
      scenario_label = scenario_meta$scenario_label[[1]],
      theta_true = scenario_theta,
      pair_success =
        !is.na(pmm_status) &
        !is.na(norm_status) &
        pmm_status == "ok" &
        norm_status == "ok" &
        is.finite(pmm_estimate) &
        is.finite(norm_estimate) &
        is.finite(pmm_se) &
        pmm_se > 0 &
        is.finite(norm_se) &
        norm_se > 0 &
        is.finite(pmm_conf_low) &
        is.finite(pmm_conf_high) &
        is.finite(norm_conf_low) &
        is.finite(norm_conf_high) &
        if (scenario_theta == 0) {
          is.finite(pmm_p_value) & is.finite(norm_p_value)
        } else {
          TRUE
        },
      delta_estimate = pmm_estimate - norm_estimate,
      delta_model_se = pmm_se - norm_se,
      delta_absolute_error =
        abs(pmm_estimate - scenario_theta) -
        abs(norm_estimate - scenario_theta),
      delta_squared_error =
        (pmm_estimate - scenario_theta)^2 -
        (norm_estimate - scenario_theta)^2,
      pmm_covered =
        pmm_conf_low <= scenario_theta &
        pmm_conf_high >= scenario_theta,
      norm_covered =
        norm_conf_low <= scenario_theta &
        norm_conf_high >= scenario_theta,
      delta_coverage_indicator =
        as.numeric(pmm_covered) - as.numeric(norm_covered),
      delta_rejection_indicator = if (scenario_theta == 0) {
        as.numeric(pmm_p_value < alpha_level) -
          as.numeric(norm_p_value < alpha_level)
      } else {
        NA_real_
      }
    ) %>%
    relocate(
      sensitivity_scenario_id,
      primary_scenario_id,
      scenario_label,
      theta_true,
      repetition
    )

  paired_repetition_list[[s]] <- paired_repetitions

  paired_summary_list[[s]] <- paired_repetitions %>%
    summarise(
      n_pairs_total = n(),
      n_pairs_success = sum(pair_success),
      mean_delta_estimate = mean(delta_estimate[pair_success]),
      mcse_mean_delta_estimate = if (sum(pair_success) > 1L) {
        sd(delta_estimate[pair_success]) / sqrt(sum(pair_success))
      } else {
        NA_real_
      },
      mean_delta_model_se = mean(delta_model_se[pair_success]),
      mean_delta_absolute_error = mean(delta_absolute_error[pair_success]),
      mean_delta_squared_error = mean(delta_squared_error[pair_success]),
      mean_delta_coverage = mean(delta_coverage_indicator[pair_success]),
      mean_delta_null_rejection = if (scenario_theta == 0) {
        mean(delta_rejection_indicator[pair_success])
      } else {
        NA_real_
      }
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
      theta_true
    )
}

paired_repetition_differences <- bind_rows(paired_repetition_list)
paired_difference_summary <- bind_rows(paired_summary_list)

fcs_pmm_atomic_save_rds(
  paired_repetition_differences,
  file.path(
    fcs_pmm_summary_dir,
    "fcs_pmm_vs_norm_paired_repetition_differences.rds"
  ),
  compress = "xz"
)

fcs_pmm_atomic_write_csv(
  paired_difference_summary,
  file.path(fcs_pmm_summary_dir, "fcs_pmm_vs_norm_paired_summary.csv")
)


# ------------------------------------------------------------------------------
# 6. Final summary QC and compact archive
# ------------------------------------------------------------------------------

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

finite_key_performance <- fcs_pmm_all_true(
  vapply(
    performance_wide[key_performance_columns],
    function(x) fcs_pmm_all_true(is.finite(x)),
    logical(1)
  )
)

summary_qc <- data.frame(
  check = c(
    "Four FCS-PMM scenarios were summarized",
    "All four scenario performance rows are present",
    "Both null-scenario Type I error rows are present",
    "All key performance summaries are finite",
    "FCS-PMM versus FCS-Norm aggregate comparison contains four rows",
    "Type I error comparison contains two rows",
    "Failure comparison contains four rows",
    "Paired comparison contains four rows",
    "Each paired comparison joins all 2,000 repetitions",
    "Canary reproduction and primary FCS-Norm reproduction remain verified"
  ),
  passed = c(
    length(scenario_results_list) == 4L,
    nrow(performance_wide) == 4L,
    nrow(type1_summary) == 2L,
    finite_key_performance,
    nrow(performance_comparison) == 4L,
    nrow(type1_comparison) == 2L,
    nrow(failure_comparison) == 4L,
    nrow(paired_difference_summary) == 4L,
    fcs_pmm_all_true(
      paired_difference_summary$n_pairs_total == nsim_primary
    ) &
      fcs_pmm_all_true(
        paired_difference_summary$n_pairs_success > 0L
      ),
    nrow(canary_reproduction_qc) == 2L &
      nrow(norm_reproduction_qc) == 2L &
      fcs_pmm_all_true(
        canary_reproduction_qc$reproduced_within_tolerance
      ) &
      fcs_pmm_all_true(
        norm_reproduction_qc$reproduced_within_tolerance
      )
  ),
  stringsAsFactors = FALSE
)

summary_qc$overall_summary_pass <- fcs_pmm_all_true(summary_qc$passed)

fcs_pmm_atomic_write_csv(
  summary_qc,
  file.path(fcs_pmm_summary_dir, "fcs_pmm_summary_qc.csv")
)

review_flags <- failure_summary %>%
  mutate(
    exceeds_prespecified_failure_threshold =
      failure_rate > deterioration_thresholds$failure_rate
  ) %>%
  filter(exceeds_prespecified_failure_threshold)

fcs_pmm_atomic_write_csv(
  review_flags,
  file.path(fcs_pmm_summary_dir, "fcs_pmm_failure_review_flags.csv")
)

fcs_pmm_atomic_save_rds(
  list(
    scenarios = fcs_pmm_sensitivity_grid,
    rsimsum_tidy = rsimsum_tidy,
    performance_wide = performance_wide,
    type1_error = type1_summary,
    failure_rates = failure_summary,
    runtime = runtime_summary,
    phase2_diagnostics = phase2_summary,
    performance_comparison = performance_comparison,
    type1_comparison = type1_comparison,
    failure_comparison = failure_comparison,
    paired_difference_summary = paired_difference_summary,
    summary_qc = summary_qc,
    failure_review_flags = review_flags
  ),
  file.path(
    fcs_pmm_summary_dir,
    "fcs_pmm_production_summary_complete.rds"
  ),
  compress = "xz"
)

if (!fcs_pmm_all_true(summary_qc$passed)) {

  failed_checks <- summary_qc$check[!summary_qc$passed]

  stop(
    "FCS-PMM production summary FAILED. Failed checks:\n- ",
    paste(failed_checks, collapse = "\n- "),
    "\nSee outputs in: ",
    fcs_pmm_summary_dir
  )
}

fcs_pmm_atomic_write_lines(
  c(
    "FCS-PMM PRODUCTION SUMMARY: PASS",
    paste("Completed:", Sys.time()),
    paste("Production scenarios:", nrow(fcs_pmm_sensitivity_grid)),
    paste("Scenario performance rows:", nrow(performance_wide)),
    paste("Null-scenario rows:", nrow(type1_summary)),
    paste("Failure-rate review flags:", nrow(review_flags))
  ),
  summary_pass_file
)

cat(
  "\nFCS-PMM production summary PASSED.\n",
  "Main comparison file:\n",
  file.path(fcs_pmm_summary_dir, "fcs_pmm_vs_norm_performance.csv"),
  "\nSummary directory: ",
  fcs_pmm_summary_dir,
  "\n",
  sep = ""
)
