################################################################################
# 19_run_aipw_strong_misspec_production.R
#
# Runs the full strong-misspecification AIPW sensitivity analysis.
#
# Each array task handles one generating setting and 2,000 repetitions:
#   task 1: MIDUS-like setting
#   task 2: stress setting
#
# Every repetition regenerates the exact primary dataset and fits four AIPW
# specifications. The misspecified nuisance models use X only, omitting A and Y.
################################################################################

source("16_aipw_strong_misspec_helpers.R")

library(dplyr)

# Identify the setting assigned to this array task.
task_id_text <- Sys.getenv(
  "SLURM_ARRAY_TASK_ID",
  unset = Sys.getenv("SIM_ARRAY_TASK_ID", unset = "")
)

if (!nzchar(task_id_text)) {
  stop(
    "No array task ID found. Use SLURM_ARRAY_TASK_ID on Torch or ",
    "SIM_ARRAY_TASK_ID for a local test."
  )
}

task_id <- suppressWarnings(as.integer(task_id_text))
valid_task_ids <- aipw_strong_settings$setting_id

if (is.na(task_id) || !(task_id %in% valid_task_ids)) {
  stop(
    "The array task ID must be one of: ",
    paste(valid_task_ids, collapse = ", "),
    "."
  )
}

setting <- aipw_strong_settings %>%
  filter(.data$setting_id == .env$task_id)

if (nrow(setting) != 1L) {
  stop(
    "Could not uniquely identify the strong-misspecification setting ",
    "for task ",
    task_id,
    "."
  )
}

# Set the production and checkpoint paths.
sensitivity_root <- file.path(
  results_dir,
  "sensitivity",
  "aipw_strong_misspecification"
)

production_dir <- file.path(
  sensitivity_root,
  "production"
)

setting_slug <- gsub(
  "[^A-Za-z0-9]+",
  "_",
  setting$setting
)

setting_dir <- file.path(
  production_dir,
  sprintf(
    "setting_%02d_%s",
    setting$setting_id,
    setting_slug
  )
)

dir.create(
  setting_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

final_file <- file.path(
  setting_dir,
  sprintf(
    "aipw_strong_setting_%02d_results.rds",
    setting$setting_id
  )
)

checkpoint_file <- file.path(
  setting_dir,
  sprintf(
    "aipw_strong_setting_%02d_checkpoint.rds",
    setting$setting_id
  )
)

if (file.exists(final_file)) {
  cat(
    "Final production file already exists; exiting:\n",
    final_file,
    "\n"
  )

  quit(save = "no", status = 0L)
}

# Read the MIDUS calibration object.
if (!file.exists(calibration_rds)) {
  stop("Calibration object not found: ", calibration_rds)
}

calibration <- readRDS(calibration_rds)

# Read the cached DGM for this setting.
cache_file <- strong_cache_filename(
  phase2_fraction = setting$phase2_fraction,
  r2_a_marker = setting$r2_a_marker,
  r2_y_marker = setting$r2_y_marker,
  theta = setting$theta
)

if (!file.exists(cache_file)) {
  stop(
    "Required primary DGM cache is missing:\n",
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
    "Validated primary scenario estimates not found:\n",
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
  filter(.data$method == "AIPW") %>%
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
  nrow(primary_aipw) != nsim_primary ||
  !setequal(
    primary_aipw$repetition,
    seq_len(nsim_primary)
  ) ||
  nrow(duplicate_primary_rows) > 0L
) {
  stop(
    "Primary AIPW comparison data must contain exactly one row for each ",
    "of the ",
    nsim_primary,
    " repetitions."
  )
}

# Initialize storage or resume from a checkpoint.
results_list <- vector("list", nsim_primary)
diagnostics_list <- vector("list", nsim_primary)
next_rep <- 1L
checkpoint_interval <- 50L

if (file.exists(checkpoint_file)) {
  checkpoint <- readRDS(checkpoint_file)

  required_checkpoint_objects <- c(
    "task_id",
    "next_rep",
    "results_list",
    "diagnostics_list"
  )

  missing_checkpoint_objects <- setdiff(
    required_checkpoint_objects,
    names(checkpoint)
  )

  if (length(missing_checkpoint_objects) > 0L) {
    stop(
      "Checkpoint is missing: ",
      paste(missing_checkpoint_objects, collapse = ", "),
      "."
    )
  }

  if (!identical(as.integer(checkpoint$task_id), task_id)) {
    stop("Checkpoint task ID does not match the current task.")
  }

  next_rep <- as.integer(checkpoint$next_rep)

  if (
    length(next_rep) != 1L ||
    is.na(next_rep) ||
    next_rep < 1L ||
    next_rep > nsim_primary + 1L
  ) {
    stop("Checkpoint contains an invalid next repetition.")
  }

  results_list <- checkpoint$results_list
  diagnostics_list <- checkpoint$diagnostics_list

  length(results_list) <- nsim_primary
  length(diagnostics_list) <- nsim_primary

  cat(
    "Resuming ",
    setting$setting,
    " at repetition ",
    next_rep,
    ".\n",
    sep = ""
  )
}

cat(
  "\nStrong AIPW misspecification production\n",
  "Task ID: ", task_id, "\n",
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
  "Repetitions: ", nsim_primary, "\n",
  "Misspecified nuisance models use X only.\n\n",
  sep = ""
)

# Run all remaining repetitions.
if (next_rep <= nsim_primary) {
  for (repetition in seq.int(next_rep, nsim_primary)) {
    rng_states <- get_strong_primary_repetition_states(
      scenario_id = setting$primary_scenario_id,
      repetition = repetition
    )

    # Regenerate the exact dataset from the primary simulation.
    .Random.seed <- rng_states$data_state

    simulated_data <- generate_two_phase_data(
      N = setting$N,
      calibration = calibration,
      dgm_parameters = dgm
    )

    analysis_data <- simulated_data$observed

    # Fit all four nuisance-model specifications to the same dataset.
    .Random.seed <- rng_states$method_state

    sensitivity_results <- run_strong_aipw_specs(
      dat = analysis_data,
      marker_names = calibration$marker_names
    ) %>%
      mutate(
        task_id = .env$task_id,
        setting_id = .env$setting$setting_id,
        setting = .env$setting$setting,
        primary_scenario_id =
          .env$setting$primary_scenario_id,
        repetition = .env$repetition,
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
        "Primary AIPW repetition was not uniquely found: ",
        repetition,
        "."
      )
    }

    # Store the primary comparison only for "Both correct".
    sensitivity_results <- sensitivity_results %>%
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

    results_list[[repetition]] <- sensitivity_results

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
        N,
        target_phase2_fraction,
        realized_phase2_fraction,
        n_phase2,
        r2_a_marker,
        r2_y_marker,
        theta_true
      )

    diagnostics_list[[repetition]] <-
      diagnostic_result

    # Save progress so an interrupted task can resume.
    if (
      repetition %% checkpoint_interval == 0L ||
      repetition == nsim_primary
    ) {
      saveRDS(
        list(
          task_id = task_id,
          setting = setting$setting,
          next_rep = repetition + 1L,
          results_list = results_list,
          diagnostics_list = diagnostics_list
        ),
        checkpoint_file
      )

      cat(
        setting$setting,
        ": completed repetition ",
        repetition,
        " of ",
        nsim_primary,
        ".\n",
        sep = ""
      )
    }
  }
}

# Combine the stored repetition-level results.
results <- bind_rows(results_list)
diagnostics <- bind_rows(diagnostics_list)

expected_repetitions <- seq_len(nsim_primary)

expected_specifications <- sort(
  aipw_strong_specs$specification
)

expected_specification_string <- paste(
  expected_specifications,
  collapse = " | "
)

expected_rows <- (
  nsim_primary * length(expected_specifications)
)

if (nrow(results) != expected_rows) {
  stop(
    "Expected ",
    expected_rows,
    " result rows; found ",
    nrow(results),
    "."
  )
}

if (!setequal(results$repetition, expected_repetitions)) {
  stop("The production repetition set is incomplete.")
}

# Confirm one row per repetition and specification.
duplicate_result_keys <- results %>%
  count(
    setting_id,
    repetition,
    specification,
    name = "n"
  ) %>%
  filter(.data$n != 1L)

if (nrow(duplicate_result_keys) > 0L) {
  stop("Duplicate result keys were detected.")
}

specification_check <- results %>%
  group_by(repetition) %>%
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
    .data$n_specifications !=
      length(expected_specifications) |
      .data$specifications !=
      expected_specification_string
  )

if (nrow(incorrect_specification_sets) > 0L) {
  stop(
    "One or more repetitions do not contain the four expected ",
    "AIPW specifications."
  )
}

# Confirm one diagnostic row per repetition.
if (nrow(diagnostics) != nsim_primary) {
  stop(
    "Expected ",
    nsim_primary,
    " diagnostic rows; found ",
    nrow(diagnostics),
    "."
  )
}

duplicate_diagnostic_keys <- diagnostics %>%
  count(setting_id, repetition, name = "n") %>%
  filter(.data$n != 1L)

if (
  !setequal(
    diagnostics$repetition,
    expected_repetitions
  ) ||
  nrow(duplicate_diagnostic_keys) > 0L
) {
  stop(
    "Diagnostic repetition keys are incomplete or duplicated."
  )
}

# Confirm exact reproduction of the primary AIPW results.
both_correct <- results %>%
  filter(.data$specification == "Both correct")

exact_tolerance <- 1e-10

primary_exact <- (
  nrow(both_correct) == nsim_primary &&
    all(
      is.finite(
        both_correct$estimate_difference_from_primary
      )
    ) &&
    all(
      is.finite(
        both_correct$se_difference_from_primary
      )
    ) &&
    max(
      abs(
        both_correct$estimate_difference_from_primary
      )
    ) < exact_tolerance &&
    max(
      abs(
        both_correct$se_difference_from_primary
      )
    ) < exact_tolerance
)

if (!primary_exact) {
  stop(
    "The both-correct specification did not exactly reproduce ",
    "the primary AIPW results."
  )
}

# Save the validated production object.
metadata <- list(
  task_id = task_id,
  setting_id = setting$setting_id,
  setting = setting$setting,
  primary_scenario_id = setting$primary_scenario_id,
  N = setting$N,
  phase2_fraction = setting$phase2_fraction,
  r2_a_marker = setting$r2_a_marker,
  r2_y_marker = setting$r2_y_marker,
  theta = setting$theta,
  nsim = nsim_primary,
  misspecification =
    "Omit both A and Y; nuisance model uses X only",
  specifications = aipw_strong_specs,
  created = as.character(Sys.time())
)

saveRDS(
  list(
    metadata = metadata,
    results = results,
    diagnostics = diagnostics
  ),
  final_file,
  compress = "xz"
)

# Remove the checkpoint only after the final file is saved.
if (file.exists(checkpoint_file)) {
  unlink(checkpoint_file)
}

n_failures <- sum(
  is.na(results$status) |
    results$status != "ok" |
    !is.finite(results$estimate) |
    !is.finite(results$se)
)

n_diagnostic_failures <- sum(
  !is.na(diagnostics$diagnostic_message) &
    nzchar(diagnostics$diagnostic_message)
)

cat(
  "\nStrong-misspecification setting complete.\n",
  "File: ", final_file, "\n",
  "Result rows: ", nrow(results), "\n",
  "AIPW failures: ", n_failures, "\n",
  "Diagnostic failures: ",
  n_diagnostic_failures,
  "\n",
  sep = ""
)
