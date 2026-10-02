################################################################################
# 40_run_ipw_weight_capping_canary.R
#
# Canary for the IPW 1st/99th-percentile weight-capping sensitivity.
#
# The new uncapped helper must reproduce both the validated fit_ipw() function
# and the stored primary IPW results. The capped fit must satisfy weight-bound,
# derivative, structural, and finite-result checks.
################################################################################

source("00_config.R")
source("02_dgm.R")
source("03_methods.R")
source("39_ipw_weight_capping_helpers.R")

library(dplyr)
library(readr)


# ------------------------------------------------------------------------------
# 1. Settings, paths, and sequential gate
# ------------------------------------------------------------------------------

nsim_canary <- read_ipw_cap_positive_integer_env(
  "SIM_IPW_CAP_CANARY_NSIM",
  2L
)

if (nsim_canary != 2L) {
  stop("The prespecified IPW-capping canary requires 2 repetitions per scenario.")
}

ipw_cap_results_dir <- file.path(results_dir, "ipw_weight_capping_sensitivity")
ipw_cap_canary_dir <- file.path(ipw_cap_results_dir, "canary")

dir.create(ipw_cap_canary_dir, showWarnings = FALSE, recursive = TRUE)

canary_pass_file <- file.path(
  ipw_cap_canary_dir,
  "IPW_WEIGHT_CAPPING_CANARY_PASS.txt"
)

if (file.exists(canary_pass_file)) {
  unlink(canary_pass_file)
}

fcs_pmm_summary_pass_file <- file.path(
  results_dir,
  "fcs_pmm_sensitivity",
  "production",
  "summary",
  "FCS_PMM_PRODUCTION_SUMMARY_PASS.txt"
)

fcs_pmm_summary_qc_file <- file.path(
  results_dir,
  "fcs_pmm_sensitivity",
  "production",
  "summary",
  "fcs_pmm_summary_qc.csv"
)

if (
  !file.exists(fcs_pmm_summary_pass_file) ||
    !file.exists(fcs_pmm_summary_qc_file)
) {
  stop(
    "The FCS-PMM production summary is not complete. Required files:\n- ",
    paste(
      c(fcs_pmm_summary_pass_file, fcs_pmm_summary_qc_file),
      collapse = "\n- "
    ),
    "\nDo not run the IPW-capping canary yet."
  )
}

if (!any(grepl(
  "FCS-PMM PRODUCTION SUMMARY: PASS",
  readLines(fcs_pmm_summary_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The FCS-PMM summary PASS file has an unexpected marker.")
}

fcs_pmm_summary_qc <- read_csv(
  fcs_pmm_summary_qc_file,
  show_col_types = FALSE
)

if (
  nrow(fcs_pmm_summary_qc) == 0L ||
    !ipw_cap_all_true(fcs_pmm_summary_qc$passed) ||
    !ipw_cap_all_true(fcs_pmm_summary_qc$overall_summary_pass)
) {
  stop("The FCS-PMM summary QC file does not record a full pass.")
}

if (!file.exists(calibration_rds)) {
  stop("Calibration object not found: ", calibration_rds)
}

calibration <- readRDS(calibration_rds)

ipw_cap_atomic_write_csv(
  ipw_cap_sensitivity_grid,
  file.path(ipw_cap_canary_dir, "ipw_weight_capping_all_scenarios.csv")
)

ipw_cap_atomic_write_csv(
  ipw_cap_canary_grid,
  file.path(ipw_cap_canary_dir, "ipw_weight_capping_canary_scenarios.csv")
)

cap_specification <- data.frame(
  lower_probability = ipw_cap_lower_probability,
  upper_probability = ipw_cap_upper_probability,
  quantile_type = ipw_cap_quantile_type,
  threshold_population = "estimated raw weights among Phase-2 observations",
  operation = "two-sided winsorization",
  stringsAsFactors = FALSE
)

ipw_cap_atomic_write_csv(
  cap_specification,
  file.path(ipw_cap_canary_dir, "ipw_weight_capping_specification.csv")
)


# ------------------------------------------------------------------------------
# 2. Verify the unchanged primary DGM caches
# ------------------------------------------------------------------------------

dgm_list <- vector("list", nrow(ipw_cap_sensitivity_grid))
dgm_qc_list <- vector("list", nrow(ipw_cap_sensitivity_grid))

for (i in seq_len(nrow(ipw_cap_sensitivity_grid))) {

  row <- ipw_cap_sensitivity_grid[i, , drop = FALSE]
  dgm_file <- ipw_cap_primary_cache_filename(row)

  if (!file.exists(dgm_file)) {
    stop("Matched primary DGM cache not found: ", dgm_file)
  }

  dgm <- readRDS(dgm_file)
  dgm_list[[i]] <- dgm
  dgm_qc_list[[i]] <- validate_ipw_cap_primary_dgm(dgm, row)
}

dgm_qc <- bind_rows(dgm_qc_list)

ipw_cap_atomic_write_csv(
  dgm_qc,
  file.path(ipw_cap_canary_dir, "ipw_weight_capping_primary_dgm_qc.csv")
)


# ------------------------------------------------------------------------------
# 3. Run uncapped and capped IPW in both canary regimes
# ------------------------------------------------------------------------------

canary_results_list <- list()
dataset_qc_list <- list()
validated_function_comparison_list <- list()
jacobian_check_list <- list()
original_result_list <- list()
result_index <- 1L

comparison_columns <- c(
  "estimate",
  "se",
  "df",
  "p_value",
  "conf_low",
  "conf_high",
  "weight_min",
  "weight_p99",
  "weight_max",
  "weight_cv",
  "weight_ess"
)

for (i in seq_len(nrow(ipw_cap_canary_grid))) {

  row <- ipw_cap_canary_grid[i, , drop = FALSE]
  dgm <- dgm_list[[row$sensitivity_scenario_id]]

  for (repetition in seq_len(nsim_canary)) {

    states <- ipw_cap_get_repetition_states(
      primary_scenario_id = row$primary_scenario_id,
      repetition = repetition
    )

    .Random.seed <- states$data_state

    simulated <- generate_two_phase_data(
      N = row$N,
      calibration = calibration,
      dgm_parameters = dgm
    )

    dat <- simulated$observed

    .Random.seed <- states$data_state

    repeated <- generate_two_phase_data(
      N = row$N,
      calibration = calibration,
      dgm_parameters = dgm
    )

    structural_qc <- ipw_cap_dataset_qc(
      dat = dat,
      marker_names = calibration$marker_names,
      row = row
    )

    original_result <- tryCatch(
      fit_ipw(dat, calibration$marker_names),
      error = function(e) failed_result("IPW", conditionMessage(e))
    )

    uncapped_result <- fit_ipw_uncapped_reference(
      dat = dat,
      marker_names = calibration$marker_names
    )

    capped_result <- fit_ipw_capped_sensitivity(
      dat = dat,
      marker_names = calibration$marker_names
    )

    original_result_list[[result_index]] <- original_result %>%
      mutate(
        primary_scenario_id = row$primary_scenario_id,
        repetition = repetition
      )

    missing_comparison_columns <- setdiff(
      comparison_columns,
      intersect(names(original_result), names(uncapped_result))
    )

    if (length(missing_comparison_columns) > 0L) {
      stop(
        "Cannot compare the validated and uncapped IPW fits; missing: ",
        paste(missing_comparison_columns, collapse = ", ")
      )
    }

    common_columns <- comparison_columns

    function_comparison <- all.equal(
      uncapped_result[, common_columns, drop = FALSE],
      original_result[, common_columns, drop = FALSE],
      tolerance = 1e-12,
      check.attributes = FALSE
    )

    numeric_differences <- unlist(
      lapply(
        common_columns,
        function(column_name) {
          difference <- abs(
            uncapped_result[[column_name]] - original_result[[column_name]]
          )
          difference[is.finite(difference)]
        }
      ),
      use.names = FALSE
    )

    validated_function_comparison_list[[result_index]] <- data.frame(
      primary_scenario_id = row$primary_scenario_id,
      repetition = repetition,
      compared_columns = paste(common_columns, collapse = ","),
      max_abs_numeric_difference = if (
        length(numeric_differences) == 0L
      ) {
        0
      } else {
        max(numeric_differences)
      },
      reproduced_within_tolerance = isTRUE(function_comparison),
      comparison_message = if (isTRUE(function_comparison)) {
        ""
      } else {
        paste(function_comparison, collapse = "; ")
      },
      stringsAsFactors = FALSE
    )

    jacobian_check_list[[result_index]] <- bind_cols(
      data.frame(
        primary_scenario_id = row$primary_scenario_id,
        repetition = repetition,
        stringsAsFactors = FALSE
      ),
      ipw_cap_jacobian_check(dat, calibration$marker_names)
    )

    dataset_qc_list[[result_index]] <- bind_cols(
      row,
      data.frame(
        repetition = repetition,
        data_rng_state = ipw_cap_state_to_string(states$data_state),
        method_rng_state = ipw_cap_state_to_string(states$method_state),
        generation_reproducible = identical(dat, repeated$observed),
        stringsAsFactors = FALSE
      ),
      structural_qc
    )

    canary_results_list[[result_index]] <- bind_rows(
      uncapped_result,
      capped_result
    ) %>%
      mutate(
        sensitivity_scenario_id = row$sensitivity_scenario_id,
        primary_scenario_id = row$primary_scenario_id,
        scenario_label = row$scenario_label,
        scenario_key = row$scenario_key,
        repetition = repetition,
        N = row$N,
        target_phase2_fraction = row$phase2_fraction,
        realized_phase2_fraction = mean(dat$phase2),
        n_phase2 = sum(dat$phase2),
        r2_a_marker = row$r2_a_marker,
        r2_y_marker = row$r2_y_marker,
        theta_true = row$theta,
        selection = dgm$selection,
        marker_error = dgm$marker_error
      )

    result_index <- result_index + 1L
  }
}

canary_results <- bind_rows(canary_results_list)
dataset_qc <- bind_rows(dataset_qc_list)
validated_function_comparison_qc <- bind_rows(
  validated_function_comparison_list
)
jacobian_qc <- bind_rows(jacobian_check_list)
original_results <- bind_rows(original_result_list)

ipw_cap_atomic_write_csv(
  canary_results,
  file.path(ipw_cap_canary_dir, "ipw_weight_capping_canary_results.csv")
)

ipw_cap_atomic_write_csv(
  dataset_qc,
  file.path(ipw_cap_canary_dir, "ipw_weight_capping_canary_dataset_qc.csv")
)

ipw_cap_atomic_write_csv(
  validated_function_comparison_qc,
  file.path(
    ipw_cap_canary_dir,
    "ipw_uncapped_validated_function_reproduction_qc.csv"
  )
)

ipw_cap_atomic_write_csv(
  jacobian_qc,
  file.path(ipw_cap_canary_dir, "ipw_weight_capping_jacobian_qc.csv")
)


# ------------------------------------------------------------------------------
# 4. Exact reproduction of the stored primary IPW results
# ------------------------------------------------------------------------------

primary_reproduction_list <- vector("list", nrow(ipw_cap_canary_grid))

for (i in seq_len(nrow(ipw_cap_canary_grid))) {

  row <- ipw_cap_canary_grid[i, , drop = FALSE]

  primary_result_file <- file.path(
    combined_dir,
    "scenarios",
    sprintf("scenario_%03d_estimates.rds", row$primary_scenario_id)
  )

  if (!file.exists(primary_result_file)) {
    stop("Matched primary result file not found: ", primary_result_file)
  }

  stored_primary <- readRDS(primary_result_file) %>%
    filter(
      repetition %in% seq_len(nsim_canary),
      method == "IPW"
    ) %>%
    arrange(repetition)

  regenerated_primary <- original_results %>%
    filter(primary_scenario_id == row$primary_scenario_id) %>%
    arrange(repetition)

  stored_columns <- c("repetition", comparison_columns)
  missing_stored_columns <- setdiff(
    stored_columns,
    intersect(names(stored_primary), names(regenerated_primary))
  )

  if (length(missing_stored_columns) > 0L) {
    stop(
      "Cannot compare regenerated and stored primary IPW results; missing: ",
      paste(missing_stored_columns, collapse = ", ")
    )
  }

  stored_keys_complete <-
    nrow(stored_primary) == nsim_canary &&
    nrow(regenerated_primary) == nsim_canary &&
    !anyDuplicated(stored_primary$repetition) &&
    !anyDuplicated(regenerated_primary$repetition) &&
    identical(stored_primary$repetition, regenerated_primary$repetition)

  comparison <- if (stored_keys_complete) {
    all.equal(
      regenerated_primary[, stored_columns, drop = FALSE],
      stored_primary[, stored_columns, drop = FALSE],
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  } else {
    "Stored and regenerated repetition keys are incomplete or unequal."
  }

  numeric_columns <- setdiff(stored_columns, "repetition")
  numeric_differences <- if (stored_keys_complete) {
    unlist(
      lapply(
        numeric_columns,
        function(column_name) {
          difference <- abs(
            regenerated_primary[[column_name]] -
              stored_primary[[column_name]]
          )
          difference[is.finite(difference)]
        }
      ),
      use.names = FALSE
    )
  } else {
    numeric()
  }

  primary_reproduction_list[[i]] <- data.frame(
    primary_scenario_id = row$primary_scenario_id,
    expected_rows = nsim_canary,
    stored_primary_rows = nrow(stored_primary),
    regenerated_rows = nrow(regenerated_primary),
    max_abs_numeric_difference = if (!stored_keys_complete) {
      NA_real_
    } else if (length(numeric_differences) == 0L) {
      0
    } else {
      max(numeric_differences)
    },
    reproduced_within_tolerance = isTRUE(comparison),
    comparison_message = if (isTRUE(comparison)) {
      ""
    } else {
      paste(comparison, collapse = "; ")
    },
    stringsAsFactors = FALSE
  )
}

primary_reproduction_qc <- bind_rows(primary_reproduction_list)

ipw_cap_atomic_write_csv(
  primary_reproduction_qc,
  file.path(ipw_cap_canary_dir, "ipw_primary_reproduction_qc.csv")
)


# ------------------------------------------------------------------------------
# 5. Final canary QC and compact archive
# ------------------------------------------------------------------------------

expected_methods <- c(ipw_uncapped_label, ipw_capped_label)
expected_result_rows <-
  nrow(ipw_cap_canary_grid) * nsim_canary * length(expected_methods)

method_keys <- canary_results %>%
  count(
    sensitivity_scenario_id,
    repetition,
    method,
    name = "n_rows"
  )

method_sets <- canary_results %>%
  group_by(sensitivity_scenario_id, repetition) %>%
  summarise(
    methods = paste(sort(unique(method)), collapse = " | "),
    n_methods = n_distinct(method),
    .groups = "drop"
  )

expected_method_set <- paste(sort(expected_methods), collapse = " | ")

capped_results <- canary_results %>%
  filter(method == ipw_capped_label)

uncapped_results <- canary_results %>%
  filter(method == ipw_uncapped_label)

raw_diagnostic_match <- capped_results %>%
  select(
    primary_scenario_id,
    repetition,
    raw_weight_min,
    raw_weight_p99,
    raw_weight_max,
    raw_weight_cv,
    raw_weight_ess
  ) %>%
  inner_join(
    uncapped_results %>%
      select(
        primary_scenario_id,
        repetition,
        weight_min,
        weight_p99,
        weight_max,
        weight_cv,
        weight_ess
      ),
    by = c("primary_scenario_id", "repetition")
  )

raw_diagnostic_keys_complete <-
  nrow(raw_diagnostic_match) ==
    nrow(ipw_cap_canary_grid) * nsim_canary &&
  !anyDuplicated(
    raw_diagnostic_match[c("primary_scenario_id", "repetition")]
  )

raw_diagnostics_equal <- raw_diagnostic_keys_complete &&
  isTRUE(all.equal(
    unname(as.matrix(
      raw_diagnostic_match %>%
        select(starts_with("raw_weight"))
    )),
    unname(as.matrix(
      raw_diagnostic_match %>%
        select(starts_with("weight"))
    )),
    tolerance = 1e-12,
    check.attributes = FALSE
  ))

cap_bounds_valid <- ipw_cap_all_true(
  capped_results$capped_weight_min >=
    capped_results$cap_lower_threshold - 1e-12 &
    capped_results$capped_weight_max <=
      capped_results$cap_upper_threshold + 1e-12 &
    abs(
      capped_results$raw_weight_p01 -
        capped_results$cap_lower_threshold
    ) <= 1e-12 &
    abs(
      capped_results$raw_weight_p99 -
        capped_results$cap_upper_threshold
    ) <= 1e-12
)

qc_checks <- data.frame(
  check = c(
    "FCS-PMM production summary passed first",
    "Four matched primary MAR/MVN scenarios are defined",
    "All primary DGM caches match the intended scenarios exactly",
    "All canary datasets are reproducible",
    "All marker missingness is blockwise",
    "phase2 agrees with marker observation",
    "All MAR probabilities are finite and strictly between zero and one",
    "Every canary dataset has enough Phase-2 observations",
    "Expected uncapped and capped IPW result rows are present",
    "Every scenario-repetition-method key appears once",
    "Both expected method labels are present",
    "All uncapped and capped IPW fits return status ok",
    "All estimates are finite",
    "All standard errors are finite and positive",
    "New uncapped helper reproduces validated fit_ipw",
    "Regenerated IPW exactly reproduces stored primary IPW",
    "Capped fit uses the prespecified 1st/99th type-7 thresholds",
    "At least one weight is capped in each tail of every canary dataset",
    "Capped weights obey their dataset-specific bounds",
    "Raw-weight diagnostics match the uncapped fit",
    "Analytic capped-weight Jacobian agrees with central differences"
  ),
  passed = c(
    ipw_cap_all_true(fcs_pmm_summary_qc$passed) &
      ipw_cap_all_true(fcs_pmm_summary_qc$overall_summary_pass),
    nrow(ipw_cap_sensitivity_grid) == 4L,
    ipw_cap_all_true(dgm_qc$dgm_matches_exactly),
    ipw_cap_all_true(dataset_qc$generation_reproducible),
    ipw_cap_all_true(dataset_qc$marker_blockwise),
    ipw_cap_all_true(dataset_qc$phase2_matches_markers),
    ipw_cap_all_true(dataset_qc$finite_valid_probabilities),
    ipw_cap_all_true(dataset_qc$sufficient_phase2_n),
    nrow(canary_results) == expected_result_rows,
    nrow(method_keys) == expected_result_rows &
      ipw_cap_all_true(method_keys$n_rows == 1L),
    nrow(method_sets) == nrow(ipw_cap_canary_grid) * nsim_canary &
      ipw_cap_all_true(method_sets$n_methods == length(expected_methods)) &
      ipw_cap_all_true(method_sets$methods == expected_method_set),
    ipw_cap_all_true(canary_results$status == "ok"),
    ipw_cap_all_true(is.finite(canary_results$estimate)),
    ipw_cap_all_true(is.finite(canary_results$se)) &
      ipw_cap_all_true(canary_results$se > 0),
    ipw_cap_all_true(
      validated_function_comparison_qc$reproduced_within_tolerance
    ),
    nrow(primary_reproduction_qc) == nrow(ipw_cap_canary_grid) &
      ipw_cap_all_true(
        primary_reproduction_qc$stored_primary_rows == nsim_canary
      ) &
      ipw_cap_all_true(
        primary_reproduction_qc$regenerated_rows == nsim_canary
      ) &
      ipw_cap_all_true(
        primary_reproduction_qc$reproduced_within_tolerance
      ),
    ipw_cap_all_true(
      capped_results$cap_lower_probability == ipw_cap_lower_probability
    ) &
      ipw_cap_all_true(
        capped_results$cap_upper_probability == ipw_cap_upper_probability
      ) &
      ipw_cap_all_true(
        capped_results$cap_quantile_type == ipw_cap_quantile_type
      ),
    ipw_cap_all_true(capped_results$n_capped_lower > 0L) &
      ipw_cap_all_true(capped_results$n_capped_upper > 0L),
    cap_bounds_valid,
    raw_diagnostics_equal,
    ipw_cap_all_true(jacobian_qc$max_scaled_difference < 1e-5)
  ),
  stringsAsFactors = FALSE
)

qc_checks$overall_canary_pass <- ipw_cap_all_true(qc_checks$passed)

ipw_cap_atomic_write_csv(
  qc_checks,
  file.path(ipw_cap_canary_dir, "ipw_weight_capping_canary_qc.csv")
)

ipw_cap_atomic_save_rds(
  list(
    settings = list(
      nsim_canary = nsim_canary,
      lower_probability = ipw_cap_lower_probability,
      upper_probability = ipw_cap_upper_probability,
      quantile_type = ipw_cap_quantile_type,
      primary_rng_seed_master = rng_seed_master
    ),
    scenarios = ipw_cap_sensitivity_grid,
    dgm_qc = dgm_qc,
    dataset_qc = dataset_qc,
    validated_function_comparison_qc = validated_function_comparison_qc,
    primary_reproduction_qc = primary_reproduction_qc,
    jacobian_qc = jacobian_qc,
    results = canary_results,
    qc = qc_checks
  ),
  file.path(
    ipw_cap_canary_dir,
    "ipw_weight_capping_canary_complete.rds"
  ),
  compress = "xz"
)

if (!ipw_cap_all_true(qc_checks$passed)) {

  failed_checks <- qc_checks$check[!qc_checks$passed]

  stop(
    "IPW weight-capping canary FAILED. Failed checks:\n- ",
    paste(failed_checks, collapse = "\n- "),
    "\nSee outputs in: ",
    ipw_cap_canary_dir
  )
}

ipw_cap_atomic_write_lines(
  c(
    "IPW WEIGHT CAPPING CANARY: PASS",
    paste("Completed:", Sys.time()),
    paste("Computational scenarios:", nrow(ipw_cap_canary_grid)),
    paste("Repetitions per scenario:", nsim_canary),
    "Cap: empirical 1st/99th percentiles among Phase-2 raw weights",
    "Quantile algorithm: R type 7",
    "Primary uncapped IPW reproduction: PASS",
    "Capped-weight Jacobian validation: PASS"
  ),
  canary_pass_file
)

cat(
  "\nIPW weight-capping canary PASSED.\n",
  "Results directory: ",
  ipw_cap_canary_dir,
  "\n",
  sep = ""
)
