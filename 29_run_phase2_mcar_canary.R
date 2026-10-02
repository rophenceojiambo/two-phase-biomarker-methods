################################################################################
# 29_run_phase2_mcar_canary.R
#
# Canary for the Phase-2 MCAR sensitivity.
#
# The full-data DGM and all six analysis methods remain unchanged. Only the
# Phase-2 selection mechanism changes from MAR to MCAR.
################################################################################

source("00_config.R")
source("02_dgm.R")
source("03_methods.R")
source("28_phase2_mcar_sensitivity_helpers.R")

library(dplyr)
library(readr)


# ------------------------------------------------------------------------------
# 1. Settings, paths, and sequential gate
# ------------------------------------------------------------------------------

nsim_canary <- read_mcar_positive_integer_env(
  "SIM_MCAR_CANARY_NSIM",
  2L
)

validation_n <- read_mcar_positive_integer_env(
  "SIM_MCAR_CANARY_VALIDATION_N",
  100000L
)

if (nsim_canary != 2L) {
  stop("The prespecified MCAR canary requires exactly 2 repetitions per computational scenario.")
}

if (validation_n < 100000L) {
  stop("SIM_MCAR_CANARY_VALIDATION_N must be at least 100000.")
}

mcar_results_dir <- file.path(results_dir, "phase2_mcar_sensitivity")
mcar_canary_dir <- file.path(mcar_results_dir, "canary")
mcar_dgm_dir <- file.path(mcar_canary_dir, "dgm")

invisible(
  lapply(
    c(mcar_results_dir, mcar_canary_dir, mcar_dgm_dir),
    dir.create,
    showWarnings = FALSE,
    recursive = TRUE
  )
)

canary_pass_file <- file.path(
  mcar_canary_dir,
  "PHASE2_MCAR_CANARY_PASS.txt"
)

if (file.exists(canary_pass_file)) {
  unlink(canary_pass_file)
}

# This gate enforces the prespecified sensitivity order. It prevents accidental
# execution while the empirical-residual production analysis is unfinished.
empirical_summary_pass_file <- file.path(
  results_dir,
  "empirical_residual_sensitivity",
  "production",
  "summary",
  "EMPIRICAL_RESIDUAL_PRODUCTION_SUMMARY_PASS.txt"
)

if (!file.exists(empirical_summary_pass_file)) {
  stop(
    "The empirical-residual production summary is not complete. Missing:\n",
    empirical_summary_pass_file,
    "\nDo not run the Phase-2 MCAR canary yet."
  )
}

if (!any(grepl(
  "PRODUCTION SUMMARY: PASS",
  readLines(empirical_summary_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The empirical-residual summary PASS file has an unexpected marker.")
}

if (!file.exists(calibration_rds)) {
  stop("Calibration object not found: ", calibration_rds)
}

calibration <- readRDS(calibration_rds)

mcar_atomic_write_csv(
  mcar_sensitivity_grid,
  file.path(mcar_canary_dir, "phase2_mcar_all_scenarios.csv")
)

mcar_atomic_write_csv(
  mcar_canary_grid,
  file.path(mcar_canary_dir, "phase2_mcar_canary_scenarios.csv")
)


# ------------------------------------------------------------------------------
# 2. Build and validate all four MCAR DGMs
# ------------------------------------------------------------------------------

dgm_list <- vector("list", nrow(mcar_sensitivity_grid))
parameter_check_list <- vector("list", nrow(mcar_sensitivity_grid))
dgm_validation_list <- vector("list", nrow(mcar_sensitivity_grid))
independence_validation_list <- vector("list", nrow(mcar_sensitivity_grid))

for (i in seq_len(nrow(mcar_sensitivity_grid))) {

  row <- mcar_sensitivity_grid[i, , drop = FALSE]
  dgm <- build_mcar_dgm(calibration, row)

  primary_cache_file <- mcar_primary_cache_filename(row)

  if (!file.exists(primary_cache_file)) {
    stop("Matching primary DGM cache not found: ", primary_cache_file)
  }

  primary_dgm <- readRDS(primary_cache_file)

  mcar_dgm_file <- mcar_dgm_cache_filename(row, mcar_dgm_dir)
  mcar_atomic_save_rds(dgm, mcar_dgm_file, compress = FALSE)

  parameter_check_list[[i]] <- compare_mcar_with_primary_dgm(
    row = row,
    mcar_dgm = dgm,
    primary_dgm = primary_dgm,
    primary_cache_file = primary_cache_file,
    mcar_cache_file = mcar_dgm_file
  )

  validation_seed <- as.integer(
    20260821L + 20000L + 100L * row$primary_scenario_id
  )

  validation <- validate_mcar_dgm(
    calibration = calibration,
    dgm = dgm,
    row = row,
    validation_n = validation_n,
    validation_seed = validation_seed
  )

  dgm_validation_list[[i]] <- validation$dgm_validation
  independence_validation_list[[i]] <-
    validation$independence_validation
  dgm_list[[i]] <- dgm
}

dgm_parameter_check <- bind_rows(parameter_check_list)
dgm_validation <- bind_rows(dgm_validation_list)
independence_validation <- bind_rows(independence_validation_list)

mcar_atomic_write_csv(
  dgm_parameter_check,
  file.path(mcar_canary_dir, "mcar_vs_primary_dgm_parameter_check.csv")
)

mcar_atomic_write_csv(
  dgm_validation,
  file.path(mcar_canary_dir, "phase2_mcar_dgm_validation.csv")
)

mcar_atomic_write_csv(
  independence_validation,
  file.path(
    mcar_canary_dir,
    "phase2_mcar_independence_validation.csv"
  )
)


# ------------------------------------------------------------------------------
# 3. Run all six methods in both canary regimes
# ------------------------------------------------------------------------------

expected_methods <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

canary_results_list <- list()
dataset_qc_list <- list()
result_index <- 1L

for (i in seq_len(nrow(mcar_canary_grid))) {

  row <- mcar_canary_grid[i, , drop = FALSE]
  dgm <- dgm_list[[row$sensitivity_scenario_id]]

  for (repetition in seq_len(nsim_canary)) {

    states <- mcar_get_repetition_states(
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

    n_marker_observed <- rowSums(
      !is.na(dat[, calibration$marker_names, drop = FALSE])
    )

    dataset_qc_list[[result_index]] <- bind_cols(
      row,
      data.frame(
        repetition = repetition,
        data_rng_state = mcar_state_to_string(states$data_state),
        method_rng_state = mcar_state_to_string(states$method_state),
        generation_reproducible = identical(dat, repeated$observed),
        marker_blockwise = all(
          n_marker_observed %in% c(0L, length(calibration$marker_names))
        ),
        phase2_matches_markers = identical(
          as.integer(
            n_marker_observed == length(calibration$marker_names)
          ),
          as.integer(dat$phase2)
        ),
        n_phase2 = sum(dat$phase2),
        realized_phase2_fraction = mean(dat$phase2),
        min_pi_true = min(dat$pi_true),
        max_pi_true = max(dat$pi_true),
        max_abs_pi_deviation = max(
          abs(dat$pi_true - row$phase2_fraction)
        ),
        sufficient_phase2_n =
          sum(dat$phase2) > length(calibration$marker_names) + 10L,
        stringsAsFactors = FALSE
      )
    )

    .Random.seed <- states$method_state
    start_time <- proc.time()[["elapsed"]]

    ans <- run_all_mcar_methods_timed(
      dat = dat,
      marker_names = calibration$marker_names
    )

    total_elapsed_seconds <- proc.time()[["elapsed"]] - start_time

    canary_results_list[[result_index]] <- ans %>%
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
        marker_error = dgm$marker_error,
        min_pi_true = min(dat$pi_true),
        max_pi_true = max(dat$pi_true),
        nimp = nimp_primary,
        mice_method = mice_method_primary,
        mice_maxit = mice_maxit_primary,
        jomo_nburn = jomo_nburn_primary,
        jomo_nbetween = jomo_nbetween_primary,
        total_all_methods_elapsed_seconds = total_elapsed_seconds
      )

    result_index <- result_index + 1L
  }
}

canary_results <- bind_rows(canary_results_list)
dataset_qc <- bind_rows(dataset_qc_list)

mcar_atomic_write_csv(
  canary_results,
  file.path(mcar_canary_dir, "phase2_mcar_canary_results.csv")
)

mcar_atomic_write_csv(
  dataset_qc,
  file.path(mcar_canary_dir, "phase2_mcar_canary_dataset_qc.csv")
)


# ------------------------------------------------------------------------------
# 4. Verify that paired full-data Naive results reproduce the primary run
# ------------------------------------------------------------------------------

naive_reproduction_list <- vector("list", nrow(mcar_canary_grid))

for (i in seq_len(nrow(mcar_canary_grid))) {

  row <- mcar_canary_grid[i, , drop = FALSE]

  primary_result_file <- file.path(
    combined_dir,
    "scenarios",
    sprintf(
      "scenario_%03d_estimates.rds",
      row$primary_scenario_id
    )
  )

  if (!file.exists(primary_result_file)) {
    stop("Matched primary result file not found: ", primary_result_file)
  }

  primary_naive <- readRDS(primary_result_file) %>%
    filter(
      repetition %in% seq_len(nsim_canary),
      method == "Naive"
    ) %>%
    arrange(repetition)

  mcar_naive <- canary_results %>%
    filter(
      primary_scenario_id == row$primary_scenario_id,
      repetition %in% seq_len(nsim_canary),
      method == "Naive"
    ) %>%
    arrange(repetition)

  comparison_columns <- intersect(
    c("repetition", "estimate", "se", "p_value", "conf_low", "conf_high"),
    intersect(names(primary_naive), names(mcar_naive))
  )

  comparison <- all.equal(
    mcar_naive[, comparison_columns, drop = FALSE],
    primary_naive[, comparison_columns, drop = FALSE],
    tolerance = 1e-12,
    check.attributes = FALSE
  )

  numeric_columns <- setdiff(comparison_columns, "repetition")

  rows_complete <-
    nrow(primary_naive) == nsim_canary &&
    nrow(mcar_naive) == nsim_canary

  keys_complete <-
    rows_complete &&
    !anyDuplicated(primary_naive$repetition) &&
    !anyDuplicated(mcar_naive$repetition) &&
    identical(primary_naive$repetition, mcar_naive$repetition)

  numeric_differences <- if (keys_complete) unlist(
    lapply(
      numeric_columns,
      function(column_name) {
        difference <- abs(
          mcar_naive[[column_name]] -
            primary_naive[[column_name]]
        )
        difference[is.finite(difference)]
      }
    ),
    use.names = FALSE
  ) else numeric()

  max_difference <- if (length(numeric_differences) == 0L) {
    NA_real_
  } else {
    max(numeric_differences)
  }

  reproduced <- rows_complete && keys_complete && isTRUE(comparison)

  naive_reproduction_list[[i]] <- data.frame(
    primary_scenario_id = row$primary_scenario_id,
    expected_rows = nsim_canary,
    primary_rows = nrow(primary_naive),
    mcar_rows = nrow(mcar_naive),
    max_abs_numeric_difference = max_difference,
    reproduced_within_tolerance = reproduced,
    comparison_message = if (reproduced) {
      ""
    } else {
      paste(comparison, collapse = "; ")
    },
    stringsAsFactors = FALSE
  )
}

naive_reproduction_qc <- bind_rows(naive_reproduction_list)

mcar_atomic_write_csv(
  naive_reproduction_qc,
  file.path(mcar_canary_dir, "phase2_mcar_primary_naive_reproduction_qc.csv")
)


# ------------------------------------------------------------------------------
# 5. Final canary QC and compact archive
# ------------------------------------------------------------------------------

method_count_qc <- canary_results %>%
  count(
    sensitivity_scenario_id,
    repetition,
    method,
    name = "n_rows"
  )

expected_result_rows <-
  nrow(mcar_canary_grid) * nsim_canary * length(expected_methods)

dataset_key_qc <- dataset_qc %>%
  count(
    sensitivity_scenario_id,
    repetition,
    name = "n_rows"
  )

method_groups <- split(
  canary_results$method,
  interaction(
    canary_results$sensitivity_scenario_id,
    canary_results$repetition,
    drop = TRUE
  )
)

all_method_sets_complete <-
  length(method_groups) == nrow(mcar_canary_grid) * nsim_canary &&
  mcar_all_true(
    vapply(
      method_groups,
      setequal,
      logical(1),
      y = expected_methods
    )
  )

qc_checks <- data.frame(
  check = c(
    "Empirical-residual sensitivity summary passed first",
    "Four intended MCAR scenarios are defined",
    "Exposure parameters exactly match the primary DGMs",
    "Outcome parameters exactly match the primary DGMs",
    "Common design values exactly match the primary DGMs",
    "Only Phase-2 selection changes from MAR to MCAR",
    "Marker-error generation remains multivariate normal",
    "MCAR needs no calibrated selection intercept",
    "All large-sample DGM checks pass",
    "Constant true Phase-2 probabilities and independence checks pass",
    "All generated datasets are reproducible",
    "Marker missingness is blockwise",
    "phase2 agrees with marker observation",
    "Every canary dataset has enough Phase-2 observations",
    "All canary true probabilities equal the target fraction",
    "Expected number of method-result rows is present",
    "Every scenario-repetition-method key appears once",
    "All six expected methods are present",
    "All methods return status ok",
    "All estimates are finite",
    "All standard errors are finite and positive",
    "All canary results record MCAR selection and MVN marker errors",
    "All canary analyses use production MI settings",
    "Paired Naive results reproduce the matched primary run"
  ),
  passed = c(
    TRUE,
    nrow(mcar_sensitivity_grid) == 4L,
    mcar_all_true(dgm_parameter_check$exposure_parameters_identical),
    mcar_all_true(dgm_parameter_check$outcome_parameters_identical),
    mcar_all_true(dgm_parameter_check$common_design_values_identical),
    mcar_all_true(dgm_parameter_check$selection_changed_as_intended),
    mcar_all_true(dgm_parameter_check$marker_error_unchanged),
    mcar_all_true(dgm_parameter_check$mcar_intercept_not_required) &
      mcar_all_true(dgm_parameter_check$mcar_expected_fraction_exact),
    mcar_all_true(dgm_validation$pass_all),
    mcar_all_true(independence_validation$pass_all),
    nrow(dataset_key_qc) == nrow(mcar_canary_grid) * nsim_canary &
      mcar_all_true(dataset_key_qc$n_rows == 1L) &
      mcar_all_true(dataset_qc$generation_reproducible),
    mcar_all_true(dataset_qc$marker_blockwise),
    mcar_all_true(dataset_qc$phase2_matches_markers),
    mcar_all_true(dataset_qc$sufficient_phase2_n),
    mcar_all_true(dataset_qc$max_abs_pi_deviation <= 1e-15),
    nrow(canary_results) == expected_result_rows,
    nrow(method_count_qc) == expected_result_rows &
      mcar_all_true(method_count_qc$n_rows == 1L),
    all_method_sets_complete,
    mcar_all_true(canary_results$status == "ok"),
    mcar_all_true(is.finite(canary_results$estimate)),
    mcar_all_true(is.finite(canary_results$se)) &
      mcar_all_true(canary_results$se > 0),
    mcar_all_true(canary_results$selection == "mcar") &
      mcar_all_true(canary_results$marker_error == "mvn"),
    mcar_all_true(canary_results$nimp == nimp_primary) &
      mcar_all_true(canary_results$mice_method == mice_method_primary) &
      mcar_all_true(canary_results$mice_maxit == mice_maxit_primary) &
      mcar_all_true(canary_results$jomo_nburn == jomo_nburn_primary) &
      mcar_all_true(canary_results$jomo_nbetween == jomo_nbetween_primary),
    nrow(naive_reproduction_qc) == nrow(mcar_canary_grid) &
      mcar_all_true(naive_reproduction_qc$primary_rows == nsim_canary) &
      mcar_all_true(naive_reproduction_qc$mcar_rows == nsim_canary) &
      mcar_all_true(naive_reproduction_qc$reproduced_within_tolerance)
  ),
  stringsAsFactors = FALSE
)

qc_checks$overall_canary_pass <- mcar_all_true(qc_checks$passed)

mcar_atomic_write_csv(
  qc_checks,
  file.path(mcar_canary_dir, "phase2_mcar_canary_qc.csv")
)

mcar_atomic_save_rds(
  list(
    settings = list(
      nsim_canary = nsim_canary,
      validation_n = validation_n,
      nimp = nimp_primary,
      mice_method = mice_method_primary,
      mice_maxit = mice_maxit_primary,
      jomo_nburn = jomo_nburn_primary,
      jomo_nbetween = jomo_nbetween_primary,
      primary_rng_seed_master = rng_seed_master
    ),
    scenarios = mcar_sensitivity_grid,
    dgm_parameter_check = dgm_parameter_check,
    dgm_validation = dgm_validation,
    independence_validation = independence_validation,
    dataset_qc = dataset_qc,
    naive_reproduction_qc = naive_reproduction_qc,
    results = canary_results,
    qc = qc_checks
  ),
  file.path(mcar_canary_dir, "phase2_mcar_canary_complete.rds"),
  compress = "xz"
)

if (!mcar_all_true(qc_checks$passed)) {

  failed_checks <- qc_checks$check[!qc_checks$passed]

  stop(
    "Phase-2 MCAR canary FAILED. Failed checks:\n- ",
    paste(failed_checks, collapse = "\n- "),
    "\nSee outputs in: ",
    mcar_canary_dir
  )
}

mcar_atomic_write_lines(
  c(
    "PHASE-2 MCAR CANARY: PASS",
    paste("Completed:", Sys.time()),
    paste("Validated DGM scenarios:", nrow(mcar_sensitivity_grid)),
    paste("Computational canary scenarios:", nrow(mcar_canary_grid)),
    paste("Canary repetitions per computational scenario:", nsim_canary),
    paste("Validation sample size per DGM:", validation_n),
    paste("Production nimp:", nimp_primary),
    paste("Production JOMO nburn:", jomo_nburn_primary),
    paste("Production JOMO nbetween:", jomo_nbetween_primary)
  ),
  canary_pass_file
)

cat(
  "\nPhase-2 MCAR canary PASSED.\n",
  "Results directory: ",
  mcar_canary_dir,
  "\n",
  sep = ""
)
