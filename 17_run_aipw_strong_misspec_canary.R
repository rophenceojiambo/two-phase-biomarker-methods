################################################################################
# 17_run_aipw_strong_misspec_canary.R
#
# Runs a small canary for the strong AIPW misspecification sensitivity.
#
# By default, the canary uses 20 exact primary datasets from each of the two
# generating settings and fits four AIPW specifications to every dataset.
#
# Run 18_check_aipw_strong_misspec_canary.R before production.
################################################################################

source("16_aipw_strong_misspec_helpers.R")

library(dplyr)
library(readr)

# Set the canary output directory.
sensitivity_root <- file.path(
  results_dir,
  "sensitivity",
  "aipw_strong_misspecification"
)

canary_dir <- file.path(
  sensitivity_root,
  "canary"
)

dir.create(
  canary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Allow a smaller local canary through an environment variable.
n_canary <- suppressWarnings(
  as.integer(
    Sys.getenv(
      "AIPW_STRONG_CANARY_N",
      unset = "20"
    )
  )
)

if (
  length(n_canary) != 1L ||
  is.na(n_canary) ||
  n_canary < 1L ||
  n_canary > nsim_primary
) {
  stop(
    "AIPW_STRONG_CANARY_N must be between 1 and ",
    nsim_primary,
    "."
  )
}

if (!file.exists(calibration_rds)) {
  stop("Calibration object not found: ", calibration_rds)
}

calibration <- readRDS(calibration_rds)

n_settings <- nrow(aipw_strong_settings)
n_datasets <- n_settings * n_canary

result_list <- vector("list", n_datasets)
diagnostic_list <- vector("list", n_datasets)

dataset_index <- 0L

# Run the canary separately for each generating setting.
for (setting_index in seq_len(n_settings)) {
  setting <- aipw_strong_settings[
    setting_index,
    ,
    drop = FALSE
  ]

  cache_file <- strong_cache_filename(
    phase2_fraction = setting$phase2_fraction,
    r2_a_marker = setting$r2_a_marker,
    r2_y_marker = setting$r2_y_marker,
    theta = setting$theta
  )

  if (!file.exists(cache_file)) {
    stop(
      "Required DGM cache not found:\n",
      cache_file,
      "\nRun 04_build_dgm_cache.R first."
    )
  }

  dgm <- readRDS(cache_file)

  # Read the validated primary AIPW results for exact comparison.
  primary_estimate_file <- file.path(
    combined_dir,
    "scenarios",
    sprintf(
      "scenario_%03d_estimates.rds",
      setting$primary_scenario_id
    )
  )

  if (!file.exists(primary_estimate_file)) {
    stop(
      "Primary estimates not found:\n",
      primary_estimate_file,
      "\nRun 06_combine_primary_results.R first."
    )
  }

  primary_results <- readRDS(primary_estimate_file)

  required_primary_columns <- c(
    "repetition",
    "method",
    "estimate",
    "se",
    "status"
  )

  missing_primary_columns <- setdiff(
    required_primary_columns,
    names(primary_results)
  )

  if (length(missing_primary_columns) > 0L) {
    stop(
      "Primary results are missing: ",
      paste(missing_primary_columns, collapse = ", "),
      "."
    )
  }

  primary_aipw <- primary_results %>%
    filter(
      .data$method == "AIPW",
      .data$repetition <= .env$n_canary
    ) %>%
    transmute(
      repetition,
      primary_estimate = estimate,
      primary_se = se,
      primary_status = status
    )

  duplicate_primary_rows <- primary_aipw %>%
    count(repetition, name = "n") %>%
    filter(.data$n != 1L)

  if (
    nrow(primary_aipw) != n_canary ||
    !setequal(
      primary_aipw$repetition,
      seq_len(n_canary)
    ) ||
    nrow(duplicate_primary_rows) > 0L
  ) {
    stop(
      "Primary AIPW results must contain exactly one row for each ",
      "canary repetition."
    )
  }

  cat(
    "\nStrong AIPW misspecification canary\n",
    "Setting: ", setting$setting, "\n",
    "Primary scenario: ", setting$primary_scenario_id, "\n",
    "N: ", setting$N, "\n",
    "Target Phase-2 fraction: ",
    setting$phase2_fraction,
    "\n",
    "Partial R2 for A and M: ",
    setting$r2_a_marker,
    "\n",
    "Partial R2 for Y and M: ",
    setting$r2_y_marker,
    "\n",
    "True exposure effect: ", setting$theta, "\n",
    "Canary repetitions: ", n_canary, "\n",
    "Misspecified nuisance models use X only.\n\n",
    sep = ""
  )

  for (repetition in seq_len(n_canary)) {
    dataset_index <- dataset_index + 1L

    rng_states <- get_strong_primary_repetition_states(
      scenario_id = setting$primary_scenario_id,
      repetition = repetition
    )

    # Regenerate the exact primary dataset.
    .Random.seed <- rng_states$data_state

    simulated_data <- generate_two_phase_data(
      N = setting$N,
      calibration = calibration,
      dgm_parameters = dgm
    )

    analysis_data <- simulated_data$observed

    # Fit all four specifications to the same dataset.
    .Random.seed <- rng_states$method_state

    robustness_results <- run_strong_aipw_specs(
      dat = analysis_data,
      marker_names = calibration$marker_names
    ) %>%
      mutate(
        setting_id = .env$setting$setting_id,
        setting = .env$setting$setting,
        primary_scenario_id =
          .env$setting$primary_scenario_id,
        repetition = .env$repetition,
        n_canary = .env$n_canary,
        N = .env$setting$N,
        target_phase2_fraction =
          .env$setting$phase2_fraction,
        realized_phase2_fraction =
          mean(.env$analysis_data$phase2),
        n_phase2 = sum(.env$analysis_data$phase2),
        r2_a_marker = .env$setting$r2_a_marker,
        r2_y_marker = .env$setting$r2_y_marker,
        theta_true = .env$setting$theta,
        min_pi_true = min(.env$analysis_data$pi_true),
        p01_pi_true = unname(
          quantile(.env$analysis_data$pi_true, 0.01)
        ),
        median_pi_true = median(
          .env$analysis_data$pi_true
        ),
        p99_pi_true = unname(
          quantile(.env$analysis_data$pi_true, 0.99)
        ),
        max_pi_true = max(.env$analysis_data$pi_true),
        data_rng_state = strong_state_to_string(
          .env$rng_states$data_state
        )
      )

    primary_row <- primary_aipw %>%
      filter(
        .data$repetition == .env$repetition
      )

    if (nrow(primary_row) != 1L) {
      stop(
        "Primary AIPW repetition not uniquely found: ",
        repetition,
        "."
      )
    }

    # Store the primary comparison only for "Both correct".
    robustness_results <- robustness_results %>%
      mutate(
        primary_estimate = if_else(
          specification == "Both correct",
          primary_row$primary_estimate[[1]],
          NA_real_
        ),
        primary_se = if_else(
          specification == "Both correct",
          primary_row$primary_se[[1]],
          NA_real_
        ),
        primary_status = if_else(
          specification == "Both correct",
          as.character(
            primary_row$primary_status[[1]]
          ),
          NA_character_
        ),
        estimate_difference_from_primary = if_else(
          specification == "Both correct",
          estimate -
            primary_row$primary_estimate[[1]],
          NA_real_
        ),
        se_difference_from_primary = if_else(
          specification == "Both correct",
          se - primary_row$primary_se[[1]],
          NA_real_
        )
      )

    result_list[[dataset_index]] <- robustness_results

    # Measure the strength of nuisance-model misspecification.
    diagnostic_result <- tryCatch(
      compute_strong_misspec_diagnostics(
        dat = analysis_data,
        marker_names = calibration$marker_names
      ),
      error = function(e) {
        data.frame(
          mean_abs_pi_difference = NA_real_,
          max_abs_pi_difference = NA_real_,
          rmse_pi_difference = NA_real_,
          cor_pi_correct_misspecified = NA_real_,
          mean_abs_marker_mean_difference = NA_real_,
          max_abs_marker_mean_difference = NA_real_,
          rmse_marker_mean_difference = NA_real_,
          relative_frobenius_sigma_difference = NA_real_,
          pi_correct_min = NA_real_,
          pi_correct_max = NA_real_,
          pi_misspecified_min = NA_real_,
          pi_misspecified_max = NA_real_,
          diagnostic_message = conditionMessage(e)
        )
      }
    )

    if (
      !("diagnostic_message" %in%
        names(diagnostic_result))
    ) {
      diagnostic_result$diagnostic_message <-
        NA_character_
    }

    diagnostic_result <- diagnostic_result %>%
      mutate(
        setting_id = .env$setting$setting_id,
        setting = .env$setting$setting,
        primary_scenario_id =
          .env$setting$primary_scenario_id,
        repetition = .env$repetition,
        n_canary = .env$n_canary,
        N = .env$setting$N,
        target_phase2_fraction =
          .env$setting$phase2_fraction,
        realized_phase2_fraction =
          mean(.env$analysis_data$phase2),
        n_phase2 = sum(.env$analysis_data$phase2),
        r2_a_marker = .env$setting$r2_a_marker,
        r2_y_marker = .env$setting$r2_y_marker,
        theta_true = .env$setting$theta
      ) %>%
      relocate(
        setting_id,
        setting,
        primary_scenario_id,
        repetition,
        n_canary,
        N,
        target_phase2_fraction,
        realized_phase2_fraction,
        n_phase2,
        r2_a_marker,
        r2_y_marker,
        theta_true
      )

    diagnostic_list[[dataset_index]] <-
      diagnostic_result

    cat(
      setting$setting,
      ": completed repetition ",
      repetition,
      " of ",
      n_canary,
      ".\n",
      sep = ""
    )
  }
}

canary_results <- bind_rows(result_list)
canary_diagnostics <- bind_rows(diagnostic_list)

# Confirm the expected number of output rows.
expected_result_rows <- (
  n_datasets * nrow(aipw_strong_specs)
)

if (nrow(canary_results) != expected_result_rows) {
  stop(
    "Expected ",
    expected_result_rows,
    " result rows; found ",
    nrow(canary_results),
    "."
  )
}

if (nrow(canary_diagnostics) != n_datasets) {
  stop(
    "Expected ",
    n_datasets,
    " diagnostic rows; found ",
    nrow(canary_diagnostics),
    "."
  )
}

# Save the canary estimates and diagnostics.
results_csv <- file.path(
  canary_dir,
  "aipw_strong_canary_results.csv"
)

results_rds <- file.path(
  canary_dir,
  "aipw_strong_canary_results.rds"
)

diagnostics_csv <- file.path(
  canary_dir,
  "aipw_strong_canary_misspecification_diagnostics.csv"
)

write_csv(canary_results, results_csv)

saveRDS(
  canary_results,
  results_rds,
  compress = "xz"
)

write_csv(
  canary_diagnostics,
  diagnostics_csv
)

cat(
  "\nStrong AIPW misspecification canary complete.\n",
  "Results: ", results_csv, "\n",
  "Diagnostics: ", diagnostics_csv, "\n",
  sep = ""
)
