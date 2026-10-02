################################################################################
# 35_run_fcs_pmm_canary.R
#
# Canary for the FCS predictive-mean-matching sensitivity.
#
# Runs the validated FCS-Norm analysis and the FCS-PMM sensitivity from the
# same method RNG state. FCS-Norm must reproduce the stored primary result.
################################################################################

source("00_config.R")
source("02_dgm.R")
source("03_methods.R")
source("34_fcs_pmm_sensitivity_helpers.R")

library(dplyr)
library(readr)


# ------------------------------------------------------------------------------
# 1. Settings, output paths, and sequential gate
# ------------------------------------------------------------------------------

nsim_canary <- read_fcs_pmm_positive_integer_env(
  "SIM_FCS_PMM_CANARY_NSIM",
  2L
)

if (nsim_canary != 2L) {
  stop("The prespecified FCS-PMM canary requires exactly 2 repetitions per scenario.")
}

fcs_pmm_results_dir <- file.path(results_dir, "fcs_pmm_sensitivity")
fcs_pmm_canary_dir <- file.path(fcs_pmm_results_dir, "canary")

dir.create(fcs_pmm_canary_dir, showWarnings = FALSE, recursive = TRUE)

canary_pass_file <- file.path(
  fcs_pmm_canary_dir,
  "FCS_PMM_CANARY_PASS.txt"
)

if (file.exists(canary_pass_file)) {
  unlink(canary_pass_file)
}

mcar_summary_pass_file <- file.path(
  results_dir,
  "phase2_mcar_sensitivity",
  "production",
  "summary",
  "PHASE2_MCAR_PRODUCTION_SUMMARY_PASS.txt"
)

mcar_summary_qc_file <- file.path(
  results_dir,
  "phase2_mcar_sensitivity",
  "production",
  "summary",
  "phase2_mcar_summary_qc.csv"
)

if (
  !file.exists(mcar_summary_pass_file) ||
    !file.exists(mcar_summary_qc_file)
) {
  stop(
    "The Phase-2 MCAR production summary is not complete. Required files:\n",
    mcar_summary_pass_file,
    "\n",
    mcar_summary_qc_file,
    "\nDo not run the FCS-PMM canary yet."
  )
}

mcar_summary_gate_pass <- any(grepl(
  "MCAR PRODUCTION SUMMARY: PASS",
  readLines(mcar_summary_pass_file, warn = FALSE),
  fixed = TRUE
))

if (!mcar_summary_gate_pass) {
  stop("The Phase-2 MCAR summary PASS file has an unexpected marker.")
}

mcar_summary_qc <- read_csv(
  mcar_summary_qc_file,
  show_col_types = FALSE
)

if (
  nrow(mcar_summary_qc) == 0L ||
    !fcs_pmm_all_true(mcar_summary_qc$passed) ||
    !fcs_pmm_all_true(mcar_summary_qc$overall_summary_pass)
) {
  stop("The Phase-2 MCAR summary QC file does not record a full pass.")
}

if (!file.exists(calibration_rds)) {
  stop("Calibration object not found: ", calibration_rds)
}

calibration <- readRDS(calibration_rds)
pmm_defaults <- fcs_pmm_default_settings()

fcs_pmm_atomic_write_csv(
  fcs_pmm_sensitivity_grid,
  file.path(fcs_pmm_canary_dir, "fcs_pmm_all_scenarios.csv")
)

fcs_pmm_atomic_write_csv(
  fcs_pmm_canary_grid,
  file.path(fcs_pmm_canary_dir, "fcs_pmm_canary_scenarios.csv")
)

fcs_pmm_atomic_write_csv(
  pmm_defaults,
  file.path(fcs_pmm_canary_dir, "fcs_pmm_default_settings.csv")
)


# ------------------------------------------------------------------------------
# 2. Verify the unchanged primary DGM caches
# ------------------------------------------------------------------------------

dgm_list <- vector("list", nrow(fcs_pmm_sensitivity_grid))
dgm_qc_list <- vector("list", nrow(fcs_pmm_sensitivity_grid))

for (i in seq_len(nrow(fcs_pmm_sensitivity_grid))) {

  row <- fcs_pmm_sensitivity_grid[i, , drop = FALSE]
  dgm_file <- fcs_pmm_primary_cache_filename(row)

  if (!file.exists(dgm_file)) {
    stop("Matched primary DGM cache not found: ", dgm_file)
  }

  dgm <- readRDS(dgm_file)
  dgm_list[[i]] <- dgm
  dgm_qc_list[[i]] <- validate_fcs_pmm_primary_dgm(dgm, row)
}

dgm_qc <- bind_rows(dgm_qc_list)

fcs_pmm_atomic_write_csv(
  dgm_qc,
  file.path(fcs_pmm_canary_dir, "fcs_pmm_primary_dgm_qc.csv")
)


# ------------------------------------------------------------------------------
# 3. Run FCS-Norm and FCS-PMM in both canary regimes
# ------------------------------------------------------------------------------

canary_results_list <- list()
dataset_qc_list <- list()
result_index <- 1L

for (i in seq_len(nrow(fcs_pmm_canary_grid))) {

  row <- fcs_pmm_canary_grid[i, , drop = FALSE]
  dgm <- dgm_list[[row$sensitivity_scenario_id]]

  for (repetition in seq_len(nsim_canary)) {

    states <- fcs_pmm_get_repetition_states(
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

    structural_qc <- fcs_pmm_dataset_qc(
      dat = dat,
      marker_names = calibration$marker_names,
      row = row
    )

    dataset_qc_list[[result_index]] <- bind_cols(
      row,
      data.frame(
        repetition = repetition,
        data_rng_state = fcs_pmm_state_to_string(states$data_state),
        method_rng_state = fcs_pmm_state_to_string(states$method_state),
        generation_reproducible = identical(dat, repeated$observed),
        stringsAsFactors = FALSE
      ),
      structural_qc
    )

    .Random.seed <- states$method_state
    norm_result <- fit_fcs_norm_reference(
      dat = dat,
      marker_names = calibration$marker_names
    )

    .Random.seed <- states$method_state
    pmm_result <- fit_fcs_pmm_sensitivity(
      dat = dat,
      marker_names = calibration$marker_names
    )

    canary_results_list[[result_index]] <- bind_rows(
      norm_result,
      pmm_result
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
        marker_error = dgm$marker_error,
        mice_method = if_else(method == fcs_norm_label, "norm", "pmm"),
        nimp = nimp_primary,
        mice_maxit = mice_maxit_primary,
        pmm_donors = pmm_defaults$donors,
        pmm_matchtype = pmm_defaults$matchtype,
        pmm_ridge = pmm_defaults$ridge,
        mice_version = pmm_defaults$mice_version
      )

    result_index <- result_index + 1L
  }
}

canary_results <- bind_rows(canary_results_list)
dataset_qc <- bind_rows(dataset_qc_list)

fcs_pmm_atomic_write_csv(
  canary_results,
  file.path(fcs_pmm_canary_dir, "fcs_pmm_canary_results.csv")
)

fcs_pmm_atomic_write_csv(
  dataset_qc,
  file.path(fcs_pmm_canary_dir, "fcs_pmm_canary_dataset_qc.csv")
)


# ------------------------------------------------------------------------------
# 4. Exact FCS-Norm reproduction of the primary simulation
# ------------------------------------------------------------------------------

norm_reproduction_list <- vector("list", nrow(fcs_pmm_canary_grid))

for (i in seq_len(nrow(fcs_pmm_canary_grid))) {

  row <- fcs_pmm_canary_grid[i, , drop = FALSE]

  primary_result_file <- file.path(
    combined_dir,
    "scenarios",
    sprintf("scenario_%03d_estimates.rds", row$primary_scenario_id)
  )

  if (!file.exists(primary_result_file)) {
    stop("Matched primary result file not found: ", primary_result_file)
  }

  primary_norm <- readRDS(primary_result_file) %>%
    filter(
      repetition %in% seq_len(nsim_canary),
      method == "FCS-MI"
    ) %>%
    arrange(repetition)

  canary_norm <- canary_results %>%
    filter(
      primary_scenario_id == row$primary_scenario_id,
      repetition %in% seq_len(nsim_canary),
      method == fcs_norm_label
    ) %>%
    arrange(repetition)

  comparison_columns <- c(
    "repetition",
    "estimate",
    "se",
    "df",
    "p_value",
    "conf_low",
    "conf_high"
  )

  missing_comparison_columns <- setdiff(
    comparison_columns,
    intersect(names(primary_norm), names(canary_norm))
  )

  if (length(missing_comparison_columns) > 0L) {
    stop(
      "FCS-Norm comparison is missing columns: ",
      paste(missing_comparison_columns, collapse = ", ")
    )
  }

  primary_values <- as.data.frame(
    primary_norm[, comparison_columns, drop = FALSE]
  )

  canary_values <- as.data.frame(
    canary_norm[, comparison_columns, drop = FALSE]
  )

  rownames(primary_values) <- NULL
  rownames(canary_values) <- NULL

  rows_complete <-
    nrow(primary_values) == nsim_canary &&
    nrow(canary_values) == nsim_canary

  keys_complete <-
    rows_complete &&
    !anyDuplicated(primary_values$repetition) &&
    !anyDuplicated(canary_values$repetition) &&
    identical(primary_values$repetition, canary_values$repetition)

  comparison <- if (keys_complete) {
    all.equal(
      canary_values,
      primary_values,
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  } else {
    "Expected row counts or repetition keys do not match."
  }

  numeric_columns <- setdiff(comparison_columns, "repetition")
  numeric_differences <- if (keys_complete) unlist(
    lapply(
      numeric_columns,
      function(column_name) {
        difference <- abs(
          canary_values[[column_name]] - primary_values[[column_name]]
        )
        difference[is.finite(difference)]
      }
    ),
    use.names = FALSE
  ) else numeric()

  reproduced <- rows_complete && keys_complete && isTRUE(comparison)

  norm_reproduction_list[[i]] <- data.frame(
    primary_scenario_id = row$primary_scenario_id,
    expected_rows = nsim_canary,
    primary_rows = nrow(primary_values),
    canary_rows = nrow(canary_values),
    max_abs_numeric_difference = if (
      length(numeric_differences) == 0L
    ) {
      NA_real_
    } else {
      max(numeric_differences)
    },
    reproduced_within_tolerance = reproduced,
    comparison_message = if (reproduced) {
      ""
    } else {
      paste(comparison, collapse = "; ")
    },
    stringsAsFactors = FALSE
  )
}

norm_reproduction_qc <- bind_rows(norm_reproduction_list)

fcs_pmm_atomic_write_csv(
  norm_reproduction_qc,
  file.path(fcs_pmm_canary_dir, "fcs_norm_primary_reproduction_qc.csv")
)


# ------------------------------------------------------------------------------
# 5. Final canary QC and compact archive
# ------------------------------------------------------------------------------

expected_methods <- c(fcs_norm_label, fcs_pmm_label)
expected_result_rows <-
  nrow(fcs_pmm_canary_grid) * nsim_canary * length(expected_methods)

method_keys <- canary_results %>%
  count(
    sensitivity_scenario_id,
    repetition,
    method,
    name = "n_rows"
  )

dataset_keys <- dataset_qc %>%
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
  length(method_groups) == nrow(fcs_pmm_canary_grid) * nsim_canary &&
  fcs_pmm_all_true(
    vapply(
      method_groups,
      setequal,
      logical(1),
      y = expected_methods
    )
  )

method_assignment_correct <- fcs_pmm_all_true(
  canary_results$mice_method ==
    ifelse(canary_results$method == fcs_norm_label, "norm", "pmm")
)

qc_checks <- data.frame(
  check = c(
    "Phase-2 MCAR production summary passed first",
    "Four matched primary MAR/MVN scenarios are defined",
    "All primary DGM caches match the intended scenarios exactly",
    "mice PMM defaults are donors 5, matchtype 1, and ridge 1e-5",
    "All canary datasets are reproducible",
    "All marker missingness is blockwise",
    "phase2 agrees with marker observation",
    "All MAR probabilities are finite and strictly between zero and one",
    "Every canary dataset has enough Phase-2 observations",
    "Expected FCS-Norm and FCS-PMM result rows are present",
    "Every scenario-repetition-method key appears once",
    "Both expected method labels are present",
    "All FCS-Norm and FCS-PMM fits return status ok",
    "All estimates are finite",
    "All standard errors are finite and positive",
    "Both analyses use 20 imputations and 10 MICE iterations",
    "Only the FCS imputation method differs",
    "FCS-Norm exactly reproduces the matched primary FCS result"
  ),
  passed = c(
    mcar_summary_gate_pass &
      fcs_pmm_all_true(mcar_summary_qc$passed) &
      fcs_pmm_all_true(mcar_summary_qc$overall_summary_pass),
    nrow(fcs_pmm_sensitivity_grid) == 4L,
    fcs_pmm_all_true(dgm_qc$dgm_matches_exactly),
    fcs_pmm_all_true(pmm_defaults$defaults_match_prespecified),
    nrow(dataset_keys) == nrow(fcs_pmm_canary_grid) * nsim_canary &
      fcs_pmm_all_true(dataset_keys$n_rows == 1L) &
      fcs_pmm_all_true(dataset_qc$generation_reproducible),
    fcs_pmm_all_true(dataset_qc$marker_blockwise),
    fcs_pmm_all_true(dataset_qc$phase2_matches_markers),
    fcs_pmm_all_true(dataset_qc$finite_valid_probabilities),
    fcs_pmm_all_true(dataset_qc$sufficient_phase2_n),
    nrow(canary_results) == expected_result_rows,
    nrow(method_keys) == expected_result_rows &
      fcs_pmm_all_true(method_keys$n_rows == 1L),
    all_method_sets_complete,
    fcs_pmm_all_true(canary_results$status == "ok"),
    fcs_pmm_all_true(is.finite(canary_results$estimate)),
    fcs_pmm_all_true(
      is.finite(canary_results$se) & canary_results$se > 0
    ),
    fcs_pmm_all_true(canary_results$nimp == nimp_primary) &
      fcs_pmm_all_true(
        canary_results$mice_maxit == mice_maxit_primary
      ),
    identical(mice_method_primary, "norm") & method_assignment_correct,
    nrow(norm_reproduction_qc) == nrow(fcs_pmm_canary_grid) &
      fcs_pmm_all_true(norm_reproduction_qc$primary_rows == nsim_canary) &
      fcs_pmm_all_true(norm_reproduction_qc$canary_rows == nsim_canary) &
      fcs_pmm_all_true(norm_reproduction_qc$reproduced_within_tolerance)
  ),
  stringsAsFactors = FALSE
)

qc_checks$overall_canary_pass <- fcs_pmm_all_true(qc_checks$passed)

fcs_pmm_atomic_write_csv(
  qc_checks,
  file.path(fcs_pmm_canary_dir, "fcs_pmm_canary_qc.csv")
)

fcs_pmm_atomic_save_rds(
  list(
    settings = list(
      nsim_canary = nsim_canary,
      nimp = nimp_primary,
      mice_maxit = mice_maxit_primary,
      primary_method = "norm",
      sensitivity_method = "pmm",
      primary_rng_seed_master = rng_seed_master
    ),
    scenarios = fcs_pmm_sensitivity_grid,
    pmm_defaults = pmm_defaults,
    dgm_qc = dgm_qc,
    dataset_qc = dataset_qc,
    norm_reproduction_qc = norm_reproduction_qc,
    results = canary_results,
    qc = qc_checks
  ),
  file.path(fcs_pmm_canary_dir, "fcs_pmm_canary_complete.rds"),
  compress = "xz"
)

if (!fcs_pmm_all_true(qc_checks$passed)) {

  failed_checks <- qc_checks$check[!qc_checks$passed]

  stop(
    "FCS-PMM canary FAILED. Failed checks:\n- ",
    paste(failed_checks, collapse = "\n- "),
    "\nSee outputs in: ",
    fcs_pmm_canary_dir
  )
}

fcs_pmm_atomic_write_lines(
  c(
    "FCS-PMM CANARY: PASS",
    paste("Completed:", Sys.time()),
    paste("Computational scenarios:", nrow(fcs_pmm_canary_grid)),
    paste("Repetitions per scenario:", nsim_canary),
    paste("Imputations:", nimp_primary),
    paste("MICE iterations:", mice_maxit_primary),
    paste("PMM donors:", pmm_defaults$donors),
    "Primary FCS-Norm reproduction: PASS"
  ),
  canary_pass_file
)

cat(
  "\nFCS-PMM canary PASSED.\n",
  "Results directory: ",
  fcs_pmm_canary_dir,
  "\n",
  sep = ""
)
