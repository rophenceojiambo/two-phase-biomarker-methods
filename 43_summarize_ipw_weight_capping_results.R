################################################################################
# 43_summarize_ipw_weight_capping_results.R
#
# ADEMP summary of capped IPW in four targeted primary MAR/MVN scenarios.
# rsimsum is the authoritative summary engine. Aggregate and paired comparisons
# are made against the completed primary untruncated IPW analysis.
################################################################################

source("00_config.R")
source("39_ipw_weight_capping_helpers.R")

library(dplyr)
library(tidyr)
library(readr)
library(rsimsum)


# ------------------------------------------------------------------------------
# 1. Paths and completed-input gate
# ------------------------------------------------------------------------------

ipw_cap_results_dir <- file.path(results_dir, "ipw_weight_capping_sensitivity")
ipw_cap_production_dir <- file.path(ipw_cap_results_dir, "production")
ipw_cap_combined_dir <- file.path(ipw_cap_production_dir, "combined")
ipw_cap_scenario_dir <- file.path(ipw_cap_combined_dir, "scenarios")
ipw_cap_summary_dir <- file.path(ipw_cap_production_dir, "summary")
ipw_cap_rsimsum_dir <- file.path(ipw_cap_summary_dir, "rsimsum_objects")

dir.create(ipw_cap_summary_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ipw_cap_rsimsum_dir, showWarnings = FALSE, recursive = TRUE)

summary_pass_file <- file.path(
  ipw_cap_summary_dir,
  "IPW_WEIGHT_CAPPING_PRODUCTION_SUMMARY_PASS.txt"
)

if (file.exists(summary_pass_file)) {
  unlink(summary_pass_file)
}

combine_pass_file <- file.path(
  ipw_cap_combined_dir,
  "IPW_WEIGHT_CAPPING_PRODUCTION_COMBINE_PASS.txt"
)

production_qc_file <- file.path(
  ipw_cap_combined_dir,
  "ipw_weight_capping_production_qc.csv"
)

canary_reproduction_file <- file.path(
  ipw_cap_combined_dir,
  "ipw_cap_canary_reproduction_qc.csv"
)

primary_reproduction_file <- file.path(
  ipw_cap_results_dir,
  "canary",
  "ipw_primary_reproduction_qc.csv"
)

validated_function_reproduction_file <- file.path(
  ipw_cap_results_dir,
  "canary",
  "ipw_uncapped_validated_function_reproduction_qc.csv"
)

required_files <- c(
  combine_pass_file,
  production_qc_file,
  canary_reproduction_file,
  primary_reproduction_file,
  validated_function_reproduction_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required IPW-capping combined files are missing:\n- ",
    paste(missing_files, collapse = "\n- "),
    "\nRun 42_combine_ipw_weight_capping_results.R first."
  )
}

if (!any(grepl(
  "IPW WEIGHT CAPPING PRODUCTION COMBINATION: PASS",
  readLines(combine_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The IPW-capping production combination PASS marker is invalid.")
}

production_qc <- read_csv(production_qc_file, show_col_types = FALSE)
canary_reproduction_qc <- read_csv(
  canary_reproduction_file,
  show_col_types = FALSE
)
primary_reproduction_qc <- read_csv(
  primary_reproduction_file,
  show_col_types = FALSE
)
validated_function_reproduction_qc <- read_csv(
  validated_function_reproduction_file,
  show_col_types = FALSE
)

if (
  nrow(production_qc) == 0L ||
    !ipw_cap_all_true(production_qc$passed) ||
    !ipw_cap_all_true(production_qc$overall_production_content_pass) ||
    !ipw_cap_all_true(
      canary_reproduction_qc$reproduced_within_tolerance
    ) ||
    !ipw_cap_all_true(
      primary_reproduction_qc$reproduced_within_tolerance
    ) ||
    !ipw_cap_all_true(
      validated_function_reproduction_qc$reproduced_within_tolerance
    )
) {
  stop("Combined IPW-capping outputs do not record a full pass.")
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

rsimsum_tidy_list <- vector("list", nrow(ipw_cap_sensitivity_grid))
type1_list <- list()
type1_index <- 1L
failure_list <- vector("list", nrow(ipw_cap_sensitivity_grid))
runtime_list <- vector("list", nrow(ipw_cap_sensitivity_grid))
phase2_list <- vector("list", nrow(ipw_cap_sensitivity_grid))
weight_list <- vector("list", nrow(ipw_cap_sensitivity_grid))
scenario_results_list <- vector("list", nrow(ipw_cap_sensitivity_grid))

for (s in seq_len(nrow(ipw_cap_sensitivity_grid))) {

  scenario_meta <- ipw_cap_sensitivity_grid[s, , drop = FALSE]
  sensitivity_scenario_id <- scenario_meta$sensitivity_scenario_id
  scenario_theta <- scenario_meta$theta[[1]]

  cat(
    "Summarising capped-IPW scenario ",
    sensitivity_scenario_id,
    " / ",
    nrow(ipw_cap_sensitivity_grid),
    " with rsimsum\n",
    sep = ""
  )

  scenario_file <- file.path(
    ipw_cap_scenario_dir,
    sprintf(
      "ipw_cap_scenario_%02d_estimates.rds",
      sensitivity_scenario_id
    )
  )

  if (!file.exists(scenario_file)) {
    stop("Capped-IPW scenario estimate file missing:\n", scenario_file)
  }

  dat <- readRDS(scenario_file)

  required_scenario_columns <- c(
    "sensitivity_scenario_id", "primary_scenario_id", "repetition", "method",
    "status", "estimate", "se", "p_value", "conf_low", "conf_high",
    "elapsed_seconds", "realized_phase2_fraction", "n_phase2",
    "min_pi_true", "p01_pi_true", "median_pi_true", "p99_pi_true",
    "max_pi_true", "cap_lower_probability", "cap_upper_probability",
    "cap_quantile_type", "cap_lower_threshold", "cap_upper_threshold",
    "n_capped_lower", "n_capped_upper", "proportion_capped_lower",
    "proportion_capped_upper", "raw_weight_min", "raw_weight_p01",
    "raw_weight_p99", "raw_weight_max", "raw_weight_cv", "raw_weight_ess",
    "capped_weight_min", "capped_weight_p01", "capped_weight_p99",
    "capped_weight_max", "capped_weight_cv", "capped_weight_ess"
  )
  missing_scenario_columns <- setdiff(
    required_scenario_columns,
    names(dat)
  )

  if (length(missing_scenario_columns) > 0L) {
    stop(
      "Capped-IPW scenario file is missing required columns: ",
      paste(missing_scenario_columns, collapse = ", ")
    )
  }

  scenario_keys_valid <-
    nrow(dat) == nsim_primary &&
    !anyDuplicated(dat[c("sensitivity_scenario_id", "repetition", "method")]) &&
    setequal(dat$repetition, seq_len(nsim_primary)) &&
    identical(unique(dat$method), ipw_capped_label) &&
    ipw_cap_all_true(
      dat$sensitivity_scenario_id == sensitivity_scenario_id
    ) &&
    ipw_cap_all_true(
      dat$primary_scenario_id == scenario_meta$primary_scenario_id[[1]]
    )

  if (!scenario_keys_valid) {
    stop("Capped-IPW scenario file has incomplete or unexpected keys.")
  }

  scenario_results_list[[s]] <- dat

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
    stop("No successful capped-IPW fits are available for this scenario.")
  }

  rs_obj <- rsimsum::simsum(
    data = analysis_dat,
    estvarname = "estimate",
    true = scenario_theta,
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

  ipw_cap_atomic_save_rds(
    rs_obj,
    file.path(
      ipw_cap_rsimsum_dir,
      sprintf("ipw_cap_scenario_%02d_rsimsum.rds", sensitivity_scenario_id)
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
      theta_true = scenario_theta,
      selection = "mar",
      marker_error = "mvn",
      cap_lower_probability = ipw_cap_lower_probability,
      cap_upper_probability = ipw_cap_upper_probability,
      cap_quantile_type = ipw_cap_quantile_type
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
      cap_lower_probability,
      cap_upper_probability,
      cap_quantile_type
    )

  rsimsum_tidy_list[[s]] <- rs_tidy

  if (isTRUE(all.equal(scenario_theta, 0))) {

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
      theta_true = scenario_theta
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
      theta_true = scenario_theta
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

  weight_list[[s]] <- analysis_dat %>%
    summarise(
      nsim_success = n(),
      mean_n_phase2 = mean(n_phase2),
      mean_cap_lower_threshold = mean(cap_lower_threshold),
      mean_cap_upper_threshold = mean(cap_upper_threshold),
      mean_n_capped_lower = mean(n_capped_lower),
      mean_n_capped_upper = mean(n_capped_upper),
      mean_proportion_capped_lower = mean(proportion_capped_lower),
      mean_proportion_capped_upper = mean(proportion_capped_upper),
      mean_raw_weight_min = mean(raw_weight_min),
      mean_raw_weight_p01 = mean(raw_weight_p01),
      mean_raw_weight_p99 = mean(raw_weight_p99),
      mean_raw_weight_max = mean(raw_weight_max),
      mean_raw_weight_cv = mean(raw_weight_cv),
      p95_raw_weight_cv = unname(quantile(raw_weight_cv, 0.95)),
      mean_raw_weight_ess = mean(raw_weight_ess),
      min_raw_weight_ess = min(raw_weight_ess),
      mean_capped_weight_min = mean(capped_weight_min),
      mean_capped_weight_p01 = mean(capped_weight_p01),
      mean_capped_weight_p99 = mean(capped_weight_p99),
      mean_capped_weight_max = mean(capped_weight_max),
      mean_capped_weight_cv = mean(capped_weight_cv),
      p95_capped_weight_cv = unname(quantile(capped_weight_cv, 0.95)),
      mean_capped_weight_ess = mean(capped_weight_ess),
      min_capped_weight_ess = min(capped_weight_ess),
      mean_delta_weight_cv = mean(capped_weight_cv - raw_weight_cv),
      mean_delta_weight_ess = mean(capped_weight_ess - raw_weight_ess)
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
      theta_true = scenario_theta,
      cap_lower_probability = ipw_cap_lower_probability,
      cap_upper_probability = ipw_cap_upper_probability,
      cap_quantile_type = ipw_cap_quantile_type
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
      cap_lower_probability,
      cap_upper_probability,
      cap_quantile_type
    )
}


# ------------------------------------------------------------------------------
# 3. Save capped-IPW performance and diagnostic summaries
# ------------------------------------------------------------------------------

rsimsum_tidy <- bind_rows(rsimsum_tidy_list)

ipw_cap_atomic_write_csv(
  rsimsum_tidy,
  file.path(ipw_cap_summary_dir, "ipw_weight_capping_rsimsum_tidy.csv")
)

ipw_cap_atomic_save_rds(
  rsimsum_tidy,
  file.path(ipw_cap_summary_dir, "ipw_weight_capping_rsimsum_tidy.rds"),
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
    cap_lower_probability,
    cap_upper_probability,
    cap_quantile_type,
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

ipw_cap_atomic_write_csv(
  performance_wide,
  file.path(ipw_cap_summary_dir, "ipw_weight_capping_performance_wide.csv")
)

ipw_cap_atomic_write_csv(
  type1_summary,
  file.path(ipw_cap_summary_dir, "ipw_weight_capping_type1_error.csv")
)

ipw_cap_atomic_write_csv(
  failure_summary,
  file.path(ipw_cap_summary_dir, "ipw_weight_capping_failure_rates.csv")
)

ipw_cap_atomic_write_csv(
  runtime_summary,
  file.path(ipw_cap_summary_dir, "ipw_weight_capping_runtime.csv")
)

ipw_cap_atomic_write_csv(
  phase2_summary,
  file.path(ipw_cap_summary_dir, "ipw_weight_capping_phase2_diagnostics.csv")
)

ipw_cap_atomic_write_csv(
  weight_summary,
  file.path(ipw_cap_summary_dir, "ipw_weight_capping_weight_diagnostics.csv")
)


# ------------------------------------------------------------------------------
# 4. Aggregate comparison with primary untruncated IPW
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
    scenario_id %in% ipw_cap_sensitivity_grid$primary_scenario_id,
    method == "IPW"
  )

if (
  nrow(primary_performance) != nrow(ipw_cap_sensitivity_grid) ||
    anyDuplicated(primary_performance$scenario_id) ||
    nrow(performance_wide) != nrow(ipw_cap_sensitivity_grid) ||
    anyDuplicated(performance_wide$primary_scenario_id)
) {
  stop("Primary or capped performance rows are incomplete or duplicated.")
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

capped_comparison_input <- performance_wide %>%
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
    ~ paste0("capped_", .x),
    all_of(comparison_metrics)
  )

uncapped_comparison_input <- primary_performance %>%
  select(scenario_id, all_of(comparison_metrics)) %>%
  rename(primary_scenario_id = scenario_id) %>%
  rename_with(
    ~ paste0("uncapped_", .x),
    all_of(comparison_metrics)
  )

performance_comparison <- capped_comparison_input %>%
  inner_join(uncapped_comparison_input, by = "primary_scenario_id") %>%
  mutate(
    delta_mean_estimate =
      capped_est__thetamean - uncapped_est__thetamean,
    delta_bias = capped_est__bias - uncapped_est__bias,
    delta_absolute_bias =
      abs(capped_est__bias) - abs(uncapped_est__bias),
    delta_empirical_se = capped_est__empse - uncapped_est__empse,
    delta_model_se = capped_est__modelse - uncapped_est__modelse,
    delta_se_ratio = capped_se_ratio - uncapped_se_ratio,
    delta_coverage = capped_est__cover - uncapped_est__cover,
    delta_bias_eliminated_coverage =
      capped_est__becover - uncapped_est__becover,
    delta_mse = capped_est__mse - uncapped_est__mse
  )

ipw_cap_atomic_write_csv(
  performance_comparison,
  file.path(
    ipw_cap_summary_dir,
    "ipw_capped_vs_uncapped_performance.csv"
  )
)

primary_type1 <- read_csv(primary_type1_file, show_col_types = FALSE) %>%
  filter(
    scenario_id %in% ipw_cap_sensitivity_grid$primary_scenario_id,
    method == "IPW"
  ) %>%
  select(
    scenario_id,
    uncapped_null_rejection_rate = null_rejection_rate,
    uncapped_mcse_null_rejection = mcse_null_rejection
  ) %>%
  rename(primary_scenario_id = scenario_id)

if (
  nrow(primary_type1) != 2L ||
    anyDuplicated(primary_type1$primary_scenario_id) ||
    nrow(type1_summary) != 2L ||
    anyDuplicated(type1_summary$primary_scenario_id)
) {
  stop("Primary or capped Type I error rows are incomplete or duplicated.")
}

type1_comparison <- type1_summary %>%
  select(
    sensitivity_scenario_id,
    primary_scenario_id,
    scenario_label,
    theta_true,
    capped_null_rejection_rate = null_rejection_rate,
    capped_mcse_null_rejection = mcse_null_rejection
  ) %>%
  inner_join(primary_type1, by = "primary_scenario_id") %>%
  mutate(
    delta_null_rejection_rate =
      capped_null_rejection_rate - uncapped_null_rejection_rate
  )

ipw_cap_atomic_write_csv(
  type1_comparison,
  file.path(ipw_cap_summary_dir, "ipw_capped_vs_uncapped_type1.csv")
)

primary_failures <- read_csv(
  primary_failure_file,
  show_col_types = FALSE
) %>%
  filter(
    scenario_id %in% ipw_cap_sensitivity_grid$primary_scenario_id,
    method == "IPW"
  ) %>%
  select(
    scenario_id,
    uncapped_n_failed = n_failed,
    uncapped_failure_rate = failure_rate,
    uncapped_mcse_failure = mcse_failure
  ) %>%
  rename(primary_scenario_id = scenario_id)

if (
  nrow(primary_failures) != nrow(ipw_cap_sensitivity_grid) ||
    anyDuplicated(primary_failures$primary_scenario_id) ||
    nrow(failure_summary) != nrow(ipw_cap_sensitivity_grid) ||
    anyDuplicated(failure_summary$primary_scenario_id)
) {
  stop("Primary or capped failure rows are incomplete or duplicated.")
}

failure_comparison <- failure_summary %>%
  select(
    sensitivity_scenario_id,
    primary_scenario_id,
    scenario_label,
    theta_true,
    capped_n_failed = n_failed,
    capped_failure_rate = failure_rate,
    capped_mcse_failure = mcse_failure
  ) %>%
  inner_join(primary_failures, by = "primary_scenario_id") %>%
  mutate(
    delta_failure_rate = capped_failure_rate - uncapped_failure_rate
  )

ipw_cap_atomic_write_csv(
  failure_comparison,
  file.path(ipw_cap_summary_dir, "ipw_capped_vs_uncapped_failures.csv")
)


# ------------------------------------------------------------------------------
# 5. Paired repetition-level capped versus untruncated IPW comparison
# ------------------------------------------------------------------------------

paired_repetition_list <- vector("list", nrow(ipw_cap_sensitivity_grid))
paired_summary_list <- vector("list", nrow(ipw_cap_sensitivity_grid))

for (s in seq_len(nrow(ipw_cap_sensitivity_grid))) {

  scenario_meta <- ipw_cap_sensitivity_grid[s, , drop = FALSE]
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

  uncapped_repetitions <- readRDS(primary_file) %>%
    filter(method == "IPW") %>%
    select(
      repetition,
      uncapped_status = status,
      uncapped_estimate = estimate,
      uncapped_se = se,
      uncapped_conf_low = conf_low,
      uncapped_conf_high = conf_high,
      uncapped_p_value = p_value,
      uncapped_weight_min = weight_min,
      uncapped_weight_p99 = weight_p99,
      uncapped_weight_max = weight_max,
      uncapped_weight_cv = weight_cv,
      uncapped_weight_ess = weight_ess
    )

  if (
    nrow(uncapped_repetitions) != nsim_primary ||
      anyDuplicated(uncapped_repetitions$repetition) ||
      !setequal(uncapped_repetitions$repetition, seq_len(nsim_primary))
  ) {
    stop("Matched primary IPW repetitions are incomplete or duplicated.")
  }

  capped_repetitions <- scenario_results_list[[s]] %>%
    select(
      repetition,
      capped_status = status,
      capped_estimate = estimate,
      capped_se = se,
      capped_conf_low = conf_low,
      capped_conf_high = conf_high,
      capped_p_value = p_value,
      cap_lower_threshold,
      cap_upper_threshold,
      n_capped_lower,
      n_capped_upper,
      raw_weight_min,
      raw_weight_p99,
      raw_weight_max,
      raw_weight_cv,
      raw_weight_ess,
      capped_weight_min,
      capped_weight_p99,
      capped_weight_max,
      capped_weight_cv,
      capped_weight_ess
    )

  if (
    nrow(capped_repetitions) != nsim_primary ||
      anyDuplicated(capped_repetitions$repetition) ||
      !setequal(capped_repetitions$repetition, seq_len(nsim_primary))
  ) {
    stop("Capped-IPW repetitions are incomplete or duplicated.")
  }

  paired_repetitions <- capped_repetitions %>%
    inner_join(uncapped_repetitions, by = "repetition") %>%
    mutate(
      sensitivity_scenario_id =
        scenario_meta$sensitivity_scenario_id[[1]],
      primary_scenario_id = scenario_meta$primary_scenario_id[[1]],
      scenario_label = scenario_meta$scenario_label[[1]],
      theta_true = scenario_theta,
      pair_success =
        !is.na(capped_status) &
        !is.na(uncapped_status) &
        capped_status == "ok" &
        uncapped_status == "ok" &
        is.finite(capped_estimate) &
        is.finite(uncapped_estimate) &
        is.finite(capped_se) &
        is.finite(uncapped_se) &
        capped_se > 0 &
        uncapped_se > 0 &
        is.finite(capped_conf_low) &
        is.finite(capped_conf_high) &
        is.finite(uncapped_conf_low) &
        is.finite(uncapped_conf_high) &
        (scenario_theta != 0 |
          (is.finite(capped_p_value) & is.finite(uncapped_p_value))),
      delta_estimate = capped_estimate - uncapped_estimate,
      delta_model_se = capped_se - uncapped_se,
      delta_absolute_error =
        abs(capped_estimate - scenario_theta) -
          abs(uncapped_estimate - scenario_theta),
      delta_squared_error =
        (capped_estimate - scenario_theta)^2 -
          (uncapped_estimate - scenario_theta)^2,
      capped_covered =
        capped_conf_low <= scenario_theta &
          capped_conf_high >= scenario_theta,
      uncapped_covered =
        uncapped_conf_low <= scenario_theta &
          uncapped_conf_high >= scenario_theta,
      delta_coverage_indicator =
        as.numeric(capped_covered) - as.numeric(uncapped_covered),
      delta_rejection_indicator = if (scenario_theta == 0) {
        as.numeric(capped_p_value < alpha_level) -
          as.numeric(uncapped_p_value < alpha_level)
      } else {
        NA_real_
      },
      raw_weight_min_difference =
        raw_weight_min - uncapped_weight_min,
      raw_weight_p99_difference =
        raw_weight_p99 - uncapped_weight_p99,
      raw_weight_max_difference =
        raw_weight_max - uncapped_weight_max,
      raw_weight_cv_difference =
        raw_weight_cv - uncapped_weight_cv,
      raw_weight_ess_difference =
        raw_weight_ess - uncapped_weight_ess,
      raw_diagnostics_match =
        is.finite(raw_weight_min_difference) &
        is.finite(raw_weight_p99_difference) &
        is.finite(raw_weight_max_difference) &
        is.finite(raw_weight_cv_difference) &
        is.finite(raw_weight_ess_difference) &
        abs(raw_weight_min_difference) < 1e-12 &
        abs(raw_weight_p99_difference) < 1e-12 &
        abs(raw_weight_max_difference) < 1e-12 &
        abs(raw_weight_cv_difference) < 1e-12 &
        abs(raw_weight_ess_difference) < 1e-12,
      delta_weight_cv = capped_weight_cv - raw_weight_cv,
      delta_weight_ess = capped_weight_ess - raw_weight_ess
    ) %>%
    relocate(
      sensitivity_scenario_id,
      primary_scenario_id,
      scenario_label,
      theta_true,
      repetition
    )

  paired_repetition_list[[s]] <- paired_repetitions

  if (
    nrow(paired_repetitions) != nsim_primary ||
      anyDuplicated(paired_repetitions$repetition) ||
      !setequal(paired_repetitions$repetition, seq_len(nsim_primary))
  ) {
    stop("Paired capped and primary IPW repetitions did not join one-to-one.")
  }

  paired_summary_list[[s]] <- paired_repetitions %>%
    summarise(
      n_pairs_total = n(),
      n_pairs_success = sum(pair_success),
      n_raw_diagnostic_matches = sum(raw_diagnostics_match),
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
      },
      mean_delta_weight_cv = mean(delta_weight_cv[pair_success]),
      mean_delta_weight_ess = mean(delta_weight_ess[pair_success]),
      max_abs_raw_weight_diagnostic_difference = max(
        abs(c(
          raw_weight_min_difference,
          raw_weight_p99_difference,
          raw_weight_max_difference,
          raw_weight_cv_difference,
          raw_weight_ess_difference
        )),
        na.rm = TRUE
      )
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

ipw_cap_atomic_save_rds(
  paired_repetition_differences,
  file.path(
    ipw_cap_summary_dir,
    "ipw_capped_vs_uncapped_paired_repetition_differences.rds"
  ),
  compress = "xz"
)

ipw_cap_atomic_write_csv(
  paired_difference_summary,
  file.path(
    ipw_cap_summary_dir,
    "ipw_capped_vs_uncapped_paired_summary.csv"
  )
)


# ------------------------------------------------------------------------------
# 6. Final summary QC, review flags, and compact archive
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

missing_key_performance_columns <- setdiff(
  key_performance_columns,
  names(performance_wide)
)

if (length(missing_key_performance_columns) > 0L) {
  stop(
    "Key performance columns are missing: ",
    paste(missing_key_performance_columns, collapse = ", ")
  )
}

finite_key_performance <- ipw_cap_all_true(
  vapply(
    performance_wide[key_performance_columns],
    function(x) ipw_cap_all_true(is.finite(x)),
    logical(1)
  )
)

summary_qc <- data.frame(
  check = c(
    "Four capped-IPW scenarios were summarized",
    "All four scenario performance rows are present",
    "Both null-scenario Type I error rows are present",
    "All key performance summaries are finite",
    "Capped versus untruncated performance comparison contains four rows",
    "Type I error comparison contains two rows",
    "Failure comparison contains four rows",
    "Weight diagnostics contain four rows",
    "Paired comparison contains four rows",
    "Each paired comparison joins all 2,000 repetitions and has valid pairs",
    "Raw-weight diagnostics reproduce primary IPW for all repetitions",
    "Canary, primary-IPW, and validated-function reproduction remain verified"
  ),
  passed = c(
    length(scenario_results_list) == 4L,
    nrow(performance_wide) == 4L,
    nrow(type1_summary) == 2L,
    finite_key_performance,
    nrow(performance_comparison) == 4L,
    nrow(type1_comparison) == 2L,
    nrow(failure_comparison) == 4L,
    nrow(weight_summary) == 4L,
    nrow(paired_difference_summary) == 4L,
    ipw_cap_all_true(
      paired_difference_summary$n_pairs_total == nsim_primary
    ) &
      ipw_cap_all_true(paired_difference_summary$n_pairs_success > 0L),
    ipw_cap_all_true(
      paired_difference_summary$n_raw_diagnostic_matches == nsim_primary
    ) &
      ipw_cap_all_true(is.finite(
        paired_difference_summary$max_abs_raw_weight_diagnostic_difference
      )) &
      ipw_cap_all_true(
      paired_difference_summary$max_abs_raw_weight_diagnostic_difference <
        1e-12
      ),
    nrow(canary_reproduction_qc) == nrow(ipw_cap_canary_grid) &
      nrow(primary_reproduction_qc) == nrow(ipw_cap_canary_grid) &
      nrow(validated_function_reproduction_qc) ==
        nrow(ipw_cap_canary_grid) * 2L &
      ipw_cap_all_true(
        canary_reproduction_qc$reproduced_within_tolerance
      ) &
      ipw_cap_all_true(
        primary_reproduction_qc$reproduced_within_tolerance
      ) &
      ipw_cap_all_true(
        validated_function_reproduction_qc$reproduced_within_tolerance
      )
  ),
  stringsAsFactors = FALSE
)

summary_qc$overall_summary_pass <- ipw_cap_all_true(summary_qc$passed)

ipw_cap_atomic_write_csv(
  summary_qc,
  file.path(ipw_cap_summary_dir, "ipw_weight_capping_summary_qc.csv")
)

review_flags <- performance_comparison %>%
  transmute(
    sensitivity_scenario_id,
    primary_scenario_id,
    scenario_label,
    theta_true,
    capped_bias = capped_est__bias,
    capped_se_ratio = capped_se_ratio,
    capped_coverage = capped_est__cover,
    delta_absolute_bias,
    delta_empirical_se,
    delta_coverage,
    delta_mse,
    absolute_bias_flag = abs(capped_est__bias) >
      deterioration_thresholds$abs_bias,
    se_ratio_flag = capped_se_ratio < deterioration_thresholds$se_ratio_low |
      capped_se_ratio > deterioration_thresholds$se_ratio_high,
    coverage_flag = capped_est__cover <
      deterioration_thresholds$coverage_low |
      capped_est__cover > deterioration_thresholds$coverage_high
  ) %>%
  left_join(
    type1_comparison %>%
      select(
        primary_scenario_id,
        capped_null_rejection_rate,
        delta_null_rejection_rate
      ),
    by = "primary_scenario_id"
  ) %>%
  mutate(
    type1_flag = !is.na(capped_null_rejection_rate) &
      capped_null_rejection_rate > deterioration_thresholds$type1_high,
    any_review_flag =
      absolute_bias_flag |
      se_ratio_flag |
      coverage_flag |
      type1_flag
  )

ipw_cap_atomic_write_csv(
  review_flags,
  file.path(ipw_cap_summary_dir, "ipw_weight_capping_review_flags.csv")
)

ipw_cap_atomic_save_rds(
  list(
    scenarios = ipw_cap_sensitivity_grid,
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
    review_flags = review_flags
  ),
  file.path(
    ipw_cap_summary_dir,
    "ipw_weight_capping_production_summary_complete.rds"
  ),
  compress = "xz"
)

if (!ipw_cap_all_true(summary_qc$passed)) {

  failed_checks <- summary_qc$check[!summary_qc$passed]

  stop(
    "IPW weight-capping production summary FAILED. Failed checks:\n- ",
    paste(failed_checks, collapse = "\n- "),
    "\nSee outputs in: ",
    ipw_cap_summary_dir
  )
}

ipw_cap_atomic_write_lines(
  c(
    "IPW WEIGHT CAPPING PRODUCTION SUMMARY: PASS",
    paste("Completed:", Sys.time()),
    paste("Production scenarios:", nrow(ipw_cap_sensitivity_grid)),
    paste("Scenario performance rows:", nrow(performance_wide)),
    paste("Null-scenario rows:", nrow(type1_summary)),
    paste("Scenarios with review flags:", sum(review_flags$any_review_flag))
  ),
  summary_pass_file
)

cat(
  "\nIPW weight-capping production summary PASSED.\n",
  "Main comparison file:\n",
  file.path(
    ipw_cap_summary_dir,
    "ipw_capped_vs_uncapped_performance.csv"
  ),
  "\nSummary directory: ",
  ipw_cap_summary_dir,
  "\n",
  sep = ""
)
