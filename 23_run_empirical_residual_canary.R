################################################################################
# 23_run_empirical_residual_canary.R
#
# Validates the empirical marker-residual resampling sensitivity before
# production.
#
# The canary:
#   1. covariance-matches the empirical residual pool to Sigma_marker
#   2. builds MIDUS-like and stress empirical-residual DGMs
#   3. compares their parameters with the matched primary DGMs
#   4. performs large-sample DGM checks
#   5. runs all six validated analysis methods
#   6. creates a pass file only when every QC check succeeds
################################################################################

source("00_config.R")
source("02_dgm.R")
source("03_methods.R")
source("22_empirical_residual_sensitivity_helpers.R")

library(dplyr)
library(readr)

# Read a positive integer from an environment variable.
read_positive_integer_env <- function(name, default) {
  value_text <- Sys.getenv(
    name,
    unset = as.character(default)
  )
  
  value <- suppressWarnings(
    as.integer(value_text)
  )
  
  if (
    length(value) != 1L ||
    is.na(value) ||
    value < 1L
  ) {
    stop(
      name,
      " must be a positive integer; received '",
      value_text,
      "'."
    )
  }
  
  value
}

# Convert a logical vector into one nonmissing pass/fail value.
all_true <- function(x) {
  isTRUE(all(x))
}

# Set canary sizes. Production MI settings are always used.
nsim_canary <- read_positive_integer_env(
  "SIM_EMPIRICAL_CANARY_NSIM",
  2L
)

validation_n <- read_positive_integer_env(
  "SIM_EMPIRICAL_CANARY_VALIDATION_N",
  100000L
)

selection_calibration_n <- read_positive_integer_env(
  "SIM_EMPIRICAL_SELECTION_CALIBRATION_N",
  50000L
)

empirical_seed_master <- 20260820L

# Set the canary output directories.
empirical_results_dir <- file.path(
  results_dir,
  "empirical_residual_sensitivity"
)

empirical_canary_dir <- file.path(
  empirical_results_dir,
  "canary"
)

empirical_dgm_dir <- file.path(
  empirical_canary_dir,
  "dgm"
)

invisible(
  lapply(
    c(
      empirical_results_dir,
      empirical_canary_dir,
      empirical_dgm_dir
    ),
    dir.create,
    showWarnings = FALSE,
    recursive = TRUE
  )
)

pass_file <- file.path(
  empirical_canary_dir,
  "EMPIRICAL_RESIDUAL_CANARY_PASS.txt"
)

# Prevent an old pass file from surviving a failed rerun.
if (file.exists(pass_file)) {
  unlink(pass_file)
}

# Read and covariance-match the empirical residual pool.
if (!file.exists(calibration_rds)) {
  stop(
    "Calibration object not found:\n",
    calibration_rds,
    "\nRun 01_calibrate_midus.R first or set SIM_CALIBRATION_RDS."
  )
}

calibration <- readRDS(calibration_rds)

empirical_preparation <- prepare_empirical_residual_calibration(
  calibration
)

empirical_calibration <- empirical_preparation$calibration

write_csv(
  empirical_preparation$pool_summary,
  file.path(
    empirical_canary_dir,
    "empirical_residual_pool_diagnostics.csv"
  )
)

write_csv(
  empirical_preparation$marginal_diagnostics,
  file.path(
    empirical_canary_dir,
    "empirical_residual_marginal_diagnostics.csv"
  )
)

write_csv(
  empirical_canary_grid,
  file.path(
    empirical_canary_dir,
    "empirical_residual_canary_scenarios.csv"
  )
)

# Return the primary DGM cache filename for one scenario.
primary_cache_filename <- function(scenario) {
  file.path(
    dgm_cache_dir,
    sprintf(
      "dgm_f%03d_rA%03d_rY%03d_th%03d.rds",
      round(100 * scenario$phase2_fraction),
      round(100 * scenario$r2_a_marker),
      round(100 * scenario$r2_y_marker),
      round(100 * scenario$theta)
    )
  )
}

# Build and validate the two alternative-scenario DGMs.
n_canary_scenarios <- nrow(empirical_canary_grid)

dgm_list <- vector("list", n_canary_scenarios)
dgm_validation_list <- vector(
  "list",
  n_canary_scenarios
)

dgm_parameter_comparison_list <- vector(
  "list",
  n_canary_scenarios
)

required_validation_columns <- c(
  "observed_r2_A_M_given_X",
  "target_r2_A_M_given_X",
  "observed_r2_Y_M_given_AX",
  "target_r2_Y_M_given_AX",
  "observed_theta_full_model",
  "target_theta",
  "mean_A",
  "sd_A",
  "mean_Y",
  "sd_Y",
  "realized_phase2_fraction",
  "target_phase2_fraction",
  "min_pi",
  "max_pi"
)

for (i in seq_len(n_canary_scenarios)) {
  scenario <- empirical_canary_grid[
    i,
    ,
    drop = FALSE
  ]
  
  selection_seed <- as.integer(
    empirical_seed_master +
      10000L +
      100L * scenario$primary_scenario_id
  )
  
  empirical_dgm <- build_dgm_parameters(
    calibration = empirical_calibration,
    r2_a_marker = scenario$r2_a_marker,
    r2_y_marker = scenario$r2_y_marker,
    theta = scenario$theta,
    phase2_fraction = scenario$phase2_fraction,
    selection = "mar",
    selection_calibration_n =
      selection_calibration_n,
    selection_seed = selection_seed,
    marker_error = "empirical"
  )
  
  if (!identical(
    empirical_dgm$marker_error,
    "empirical"
  )) {
    stop(
      "Canary DGM was not built with empirical marker errors."
    )
  }
  
  # Compare with the matched primary multivariate-normal DGM.
  primary_cache_file <- primary_cache_filename(
    scenario
  )
  
  if (!file.exists(primary_cache_file)) {
    stop(
      "Matching primary DGM cache not found:\n",
      primary_cache_file,
      "\nRun 04_build_dgm_cache.R first."
    )
  }
  
  primary_dgm <- readRDS(primary_cache_file)
  
  exposure_parameters_identical <- isTRUE(
    all.equal(
      empirical_dgm$exposure,
      primary_dgm$exposure,
      tolerance = 0,
      check.attributes = TRUE
    )
  )
  
  outcome_parameters_identical <- isTRUE(
    all.equal(
      empirical_dgm$outcome,
      primary_dgm$outcome,
      tolerance = 0,
      check.attributes = TRUE
    )
  )
  
  design_values_identical <- all_true(
    c(
      isTRUE(all.equal(
        empirical_dgm$r2_a_marker,
        primary_dgm$r2_a_marker,
        tolerance = 0
      )),
      isTRUE(all.equal(
        empirical_dgm$r2_y_marker,
        primary_dgm$r2_y_marker,
        tolerance = 0
      )),
      isTRUE(all.equal(
        empirical_dgm$theta,
        primary_dgm$theta,
        tolerance = 0
      )),
      isTRUE(all.equal(
        empirical_dgm$phase2_fraction,
        primary_dgm$phase2_fraction,
        tolerance = 0
      )),
      identical(
        empirical_dgm$selection,
        primary_dgm$selection
      )
    )
  )
  
  marker_error_changed_as_intended <- (
    identical(primary_dgm$marker_error, "mvn") &&
      identical(
        empirical_dgm$marker_error,
        "empirical"
      )
  )
  
  primary_selection_intercept <- (
    primary_dgm$selection_parameters$eta0
  )
  
  empirical_selection_intercept <- (
    empirical_dgm$selection_parameters$eta0
  )
  
  selection_intercept_recalibrated <- !isTRUE(
    all.equal(
      primary_selection_intercept,
      empirical_selection_intercept,
      tolerance = 1e-12
    )
  )
  
  dgm_parameter_comparison_list[[i]] <- bind_cols(
    scenario,
    data.frame(
      primary_cache_file =
        basename(primary_cache_file),
      exposure_parameters_identical =
        exposure_parameters_identical,
      outcome_parameters_identical =
        outcome_parameters_identical,
      design_values_identical =
        design_values_identical,
      primary_marker_error =
        primary_dgm$marker_error,
      sensitivity_marker_error =
        empirical_dgm$marker_error,
      marker_error_changed_as_intended =
        marker_error_changed_as_intended,
      primary_selection_intercept =
        primary_selection_intercept,
      empirical_selection_intercept =
        empirical_selection_intercept,
      selection_intercept_recalibrated =
        selection_intercept_recalibrated,
      stringsAsFactors = FALSE
    )
  )
  
  dgm_list[[i]] <- empirical_dgm
  
  saveRDS(
    empirical_dgm,
    file.path(
      empirical_dgm_dir,
      paste0(
        "empirical_dgm_",
        scenario$scenario_key,
        ".rds"
      )
    )
  )
  
  # Check whether a large generated sample recovers the DGM targets.
  validation_seed <- as.integer(
    empirical_seed_master +
      20000L +
      100L * scenario$primary_scenario_id
  )
  
  dgm_validation <- check_dgm(
    calibration = empirical_calibration,
    dgm_parameters = empirical_dgm,
    N_check = validation_n,
    seed = validation_seed
  )
  
  missing_validation_columns <- setdiff(
    required_validation_columns,
    names(dgm_validation)
  )
  
  if (length(missing_validation_columns) > 0L) {
    stop(
      "check_dgm() output is missing: ",
      paste(
        missing_validation_columns,
        collapse = ", "
      ),
      "."
    )
  }
  
  dgm_validation$selection_target_check <- (
    empirical_dgm$selection_parameters$
      expected_fraction_check
  )
  
  dgm_validation$dgm_marker_error <- (
    empirical_dgm$marker_error
  )
  
  dgm_validation_list[[i]] <- bind_cols(
    scenario[
      rep(1L, nrow(dgm_validation)),
      ,
      drop = FALSE
    ],
    dgm_validation
  )
}

dgm_parameter_comparison <- bind_rows(
  dgm_parameter_comparison_list
)

write_csv(
  dgm_parameter_comparison,
  file.path(
    empirical_canary_dir,
    "empirical_vs_primary_dgm_parameter_check.csv"
  )
)

# Apply prespecified tolerances to the large-sample checks.
dgm_validation <- bind_rows(
  dgm_validation_list
) %>%
  mutate(
    abs_error_r2_A = abs(
      observed_r2_A_M_given_X -
        target_r2_A_M_given_X
    ),
    abs_error_r2_Y = abs(
      observed_r2_Y_M_given_AX -
        target_r2_Y_M_given_AX
    ),
    abs_error_theta = abs(
      observed_theta_full_model -
        target_theta
    ),
    abs_error_phase2_fraction = abs(
      realized_phase2_fraction -
        target_phase2_fraction
    ),
    abs_error_selection_target = abs(
      selection_target_check -
        target_phase2_fraction
    ),
    pass_r2_A = abs_error_r2_A <= 0.01,
    pass_r2_Y = abs_error_r2_Y <= 0.01,
    pass_theta = abs_error_theta <= 0.01,
    pass_mean_A = abs(mean_A) <= 0.02,
    pass_sd_A = abs(sd_A - 1) <= 0.02,
    pass_mean_Y = abs(mean_Y) <= 0.02,
    pass_sd_Y = abs(sd_Y - 1) <= 0.02,
    pass_phase2_fraction =
      abs_error_phase2_fraction <= 0.01,
    pass_selection_target =
      abs_error_selection_target <= 1e-8,
    pass_probability_range = (
      is.finite(min_pi) &
        is.finite(max_pi) &
        min_pi > 0 &
        max_pi < 1
    ),
    pass_marker_error = (
      dgm_marker_error == "empirical"
    ),
    pass_all = (
      pass_r2_A &
        pass_r2_Y &
        pass_theta &
        pass_mean_A &
        pass_sd_A &
        pass_mean_Y &
        pass_sd_Y &
        pass_phase2_fraction &
        pass_selection_target &
        pass_probability_range &
        pass_marker_error
    )
  )

write_csv(
  dgm_validation,
  file.path(
    empirical_canary_dir,
    "empirical_residual_dgm_validation.csv"
  )
)

# Run all six validated methods using production MI settings.
RNGkind("L'Ecuyer-CMRG")

expected_methods <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

n_canary_datasets <- (
  n_canary_scenarios * nsim_canary
)

canary_results_list <- vector(
  "list",
  n_canary_datasets
)

dataset_qc_list <- vector(
  "list",
  n_canary_datasets
)

dataset_index <- 0L

for (i in seq_len(n_canary_scenarios)) {
  scenario <- empirical_canary_grid[
    i,
    ,
    drop = FALSE
  ]
  
  empirical_dgm <- dgm_list[[i]]
  
  for (repetition in seq_len(nsim_canary)) {
    dataset_index <- dataset_index + 1L
    
    data_seed <- as.integer(
      empirical_seed_master +
        100000L +
        1000L * scenario$primary_scenario_id +
        2L * repetition
    )
    
    method_seed <- data_seed + 1L
    
    # Generate one empirical-residual two-phase dataset.
    set.seed(data_seed)
    
    simulated_data <- generate_two_phase_data(
      N = scenario$N,
      calibration = empirical_calibration,
      dgm_parameters = empirical_dgm
    )
    
    analysis_data <- simulated_data$observed
    
    # Regenerate from the same seed to check deterministic generation.
    set.seed(data_seed)
    
    reproduced_data <- generate_two_phase_data(
      N = scenario$N,
      calibration = empirical_calibration,
      dgm_parameters = empirical_dgm
    )
    
    generation_reproducible <- identical(
      analysis_data,
      reproduced_data$observed
    )
    
    n_markers_observed <- rowSums(
      !is.na(
        analysis_data[
          ,
          empirical_calibration$marker_names,
          drop = FALSE
        ]
      )
    )
    
    n_markers <- length(
      empirical_calibration$marker_names
    )
    
    marker_blockwise <- all(
      n_markers_observed %in% c(0L, n_markers)
    )
    
    phase2_matches_markers <- identical(
      as.integer(n_markers_observed == n_markers),
      as.integer(analysis_data$phase2)
    )
    
    dataset_qc_list[[dataset_index]] <- bind_cols(
      scenario,
      data.frame(
        repetition = repetition,
        data_seed = data_seed,
        method_seed = method_seed,
        generation_reproducible =
          generation_reproducible,
        marker_blockwise = marker_blockwise,
        phase2_matches_markers =
          phase2_matches_markers,
        n_phase2 = sum(analysis_data$phase2),
        realized_phase2_fraction =
          mean(analysis_data$phase2),
        min_pi_true = min(analysis_data$pi_true),
        max_pi_true = max(analysis_data$pi_true),
        sufficient_phase2_n = (
          sum(analysis_data$phase2) >
            n_markers + 10L
        ),
        finite_probability_range = (
          all(is.finite(analysis_data$pi_true)) &
            all(analysis_data$pi_true > 0) &
            all(analysis_data$pi_true < 1)
        ),
        stringsAsFactors = FALSE
      )
    )
    
    # Analyze the fixed dataset with all six methods.
    set.seed(method_seed)
    
    start_time <- proc.time()[["elapsed"]]
    
    method_results <- run_all_methods(
      dat = analysis_data,
      marker_names =
        empirical_calibration$marker_names,
      nimp = nimp_primary,
      mice_maxit = mice_maxit_primary,
      jomo_nburn = jomo_nburn_primary,
      jomo_nbetween = jomo_nbetween_primary
    )
    
    total_elapsed_seconds <- (
      proc.time()[["elapsed"]] -
        start_time
    )
    
    method_results <- method_results %>%
      mutate(
        sensitivity_scenario_id =
          scenario$sensitivity_scenario_id,
        primary_scenario_id =
          scenario$primary_scenario_id,
        scenario_label = scenario$scenario_label,
        scenario_key = scenario$scenario_key,
        repetition = repetition,
        data_seed = data_seed,
        method_seed = method_seed,
        N = scenario$N,
        target_phase2_fraction =
          scenario$phase2_fraction,
        realized_phase2_fraction =
          mean(analysis_data$phase2),
        n_phase2 = sum(analysis_data$phase2),
        r2_a_marker = scenario$r2_a_marker,
        r2_y_marker = scenario$r2_y_marker,
        theta_true = scenario$theta,
        marker_error = empirical_dgm$marker_error,
        nimp = nimp_primary,
        mice_maxit = mice_maxit_primary,
        jomo_nburn = jomo_nburn_primary,
        jomo_nbetween = jomo_nbetween_primary,
        total_all_methods_elapsed_seconds =
          total_elapsed_seconds
      )
    
    canary_results_list[[dataset_index]] <-
      method_results
  }
}

canary_results <- bind_rows(
  canary_results_list
)

dataset_qc <- bind_rows(
  dataset_qc_list
)

results_csv <- file.path(
  empirical_canary_dir,
  "empirical_residual_canary_results.csv"
)

dataset_qc_csv <- file.path(
  empirical_canary_dir,
  "empirical_residual_canary_dataset_qc.csv"
)

write_csv(canary_results, results_csv)
write_csv(dataset_qc, dataset_qc_csv)

# Save all canary components together for reproducibility.
saveRDS(
  list(
    settings = list(
      nsim_canary = nsim_canary,
      validation_n = validation_n,
      selection_calibration_n =
        selection_calibration_n,
      nimp = nimp_primary,
      mice_maxit = mice_maxit_primary,
      jomo_nburn = jomo_nburn_primary,
      jomo_nbetween = jomo_nbetween_primary,
      empirical_seed_master = empirical_seed_master
    ),
    scenarios = empirical_canary_grid,
    residual_pool_summary =
      empirical_preparation$pool_summary,
    residual_marginal_diagnostics =
      empirical_preparation$marginal_diagnostics,
    dgm_parameter_comparison =
      dgm_parameter_comparison,
    dgm_validation = dgm_validation,
    dataset_qc = dataset_qc,
    results = canary_results
  ),
  file.path(
    empirical_canary_dir,
    "empirical_residual_canary_complete.rds"
  )
)

# Make the final canary decision.
method_key_counts <- canary_results %>%
  count(
    sensitivity_scenario_id,
    repetition,
    method,
    name = "n_rows"
  )

expected_result_rows <- (
  n_canary_datasets * length(expected_methods)
)

qc_checks <- data.frame(
  check = c(
    "Residual-pool means match zero",
    "Residual-pool covariance matches Sigma_marker",
    "Exposure parameters exactly match the primary DGM",
    "Outcome parameters exactly match the primary DGM",
    "Design values exactly match the primary DGM",
    "Only marker-error generation changes as intended",
    "Phase-2 intercept is recalibrated for empirical residuals",
    "All large-sample DGM checks pass",
    "All generated datasets are reproducible",
    "Marker missingness is blockwise",
    "phase2 agrees with marker observation",
    "Every canary dataset has enough Phase-2 observations",
    "All true selection probabilities are finite and in (0,1)",
    "Expected number of result rows is present",
    "Every scenario-repetition-method combination appears once",
    "All six expected methods are present",
    "All methods return status ok",
    "All estimates are finite",
    "All standard errors are finite and positive",
    "All canary results record empirical marker errors",
    "All canary analyses use production MI settings"
  ),
  passed = c(
    isTRUE(
      empirical_preparation$pool_summary$
        max_abs_matched_mean <= 1e-10
    ),
    isTRUE(
      empirical_preparation$pool_summary$
        max_relative_covariance_error_after <= 1e-8
    ),
    all_true(
      dgm_parameter_comparison$
        exposure_parameters_identical
    ),
    all_true(
      dgm_parameter_comparison$
        outcome_parameters_identical
    ),
    all_true(
      dgm_parameter_comparison$
        design_values_identical
    ),
    all_true(
      dgm_parameter_comparison$
        marker_error_changed_as_intended
    ),
    all_true(
      dgm_parameter_comparison$
        selection_intercept_recalibrated
    ),
    all_true(dgm_validation$pass_all),
    all_true(dataset_qc$generation_reproducible),
    all_true(dataset_qc$marker_blockwise),
    all_true(dataset_qc$phase2_matches_markers),
    all_true(dataset_qc$sufficient_phase2_n),
    all_true(dataset_qc$finite_probability_range),
    nrow(canary_results) == expected_result_rows,
    (
      nrow(method_key_counts) ==
        expected_result_rows &&
        all_true(method_key_counts$n_rows == 1L)
    ),
    setequal(
      unique(canary_results$method),
      expected_methods
    ),
    all_true(
      !is.na(canary_results$status) &
        canary_results$status == "ok"
    ),
    all_true(is.finite(canary_results$estimate)),
    all_true(
      is.finite(canary_results$se) &
        canary_results$se > 0
    ),
    all_true(
      !is.na(canary_results$marker_error) &
        canary_results$marker_error == "empirical"
    ),
    all_true(
      canary_results$nimp == nimp_primary &
        canary_results$mice_maxit ==
        mice_maxit_primary &
        canary_results$jomo_nburn ==
        jomo_nburn_primary &
        canary_results$jomo_nbetween ==
        jomo_nbetween_primary
    )
  ),
  stringsAsFactors = FALSE
)

qc_checks$overall_canary_pass <- all_true(
  qc_checks$passed
)

write_csv(
  qc_checks,
  file.path(
    empirical_canary_dir,
    "empirical_residual_canary_qc.csv"
  )
)

if (!all_true(qc_checks$passed)) {
  failed_checks <- qc_checks$check[
    !qc_checks$passed
  ]
  
  stop(
    "Empirical residual-resampling canary failed. Failed checks:\n- ",
    paste(failed_checks, collapse = "\n- "),
    "\nSee outputs in: ",
    empirical_canary_dir
  )
}

# Create the production authorization file only after every check passes.
writeLines(
  c(
    "EMPIRICAL RESIDUAL-RESAMPLING CANARY: PASS",
    paste("Completed:", Sys.time()),
    paste(
      "Canary repetitions per scenario:",
      nsim_canary
    ),
    paste(
      "Validation sample size per DGM:",
      validation_n
    ),
    paste("Production nimp:", nimp_primary),
    paste(
      "Production jomo nburn:",
      jomo_nburn_primary
    ),
    paste(
      "Production jomo nbetween:",
      jomo_nbetween_primary
    )
  ),
  pass_file
)

cat(
  "\nEmpirical residual-resampling canary passed.\n",
  "Results: ", results_csv, "\n",
  "Dataset QC: ", dataset_qc_csv, "\n",
  "Pass file: ", pass_file, "\n",
  sep = ""
)

