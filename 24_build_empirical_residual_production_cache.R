################################################################################
# 24_build_empirical_residual_production_cache.R
#
# Builds and validates the four empirical-residual production DGMs.
#
# The script:
#   1. requires a completed empirical-residual canary
#   2. saves the covariance-matched empirical calibration
#   3. builds the four null and alternative DGMs
#   4. compares them with their matched primary DGMs
#   5. performs large-sample DGM checks
#   6. creates a production-cache pass file only after all checks succeed
################################################################################

source("00_config.R")
source("02_dgm.R")
source("22_empirical_residual_sensitivity_helpers.R")

library(dplyr)
library(readr)

read_positive_integer_env <- function(name, default) {
  value_text <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.integer(value_text))
  
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

all_true <- function(x) {
  isTRUE(all(x))
}

# Write cache files through temporary files to avoid partial outputs.
atomic_save_rds <- function(
    object,
    path,
    compress = FALSE
) {
  temporary_path <- paste0(
    path,
    ".tmp_",
    Sys.getpid()
  )
  
  on.exit(
    unlink(temporary_path),
    add = TRUE
  )
  
  saveRDS(
    object,
    temporary_path,
    compress = compress
  )
  
  if (!file.rename(temporary_path, path)) {
    stop("Could not publish RDS file: ", path)
  }
  
  invisible(path)
}

atomic_write_csv <- function(object, path) {
  temporary_path <- paste0(
    path,
    ".tmp_",
    Sys.getpid()
  )
  
  on.exit(
    unlink(temporary_path),
    add = TRUE
  )
  
  write_csv(object, temporary_path)
  
  if (!file.rename(temporary_path, path)) {
    stop("Could not publish CSV file: ", path)
  }
  
  invisible(path)
}

atomic_write_lines <- function(text, path) {
  temporary_path <- paste0(
    path,
    ".tmp_",
    Sys.getpid()
  )
  
  on.exit(
    unlink(temporary_path),
    add = TRUE
  )
  
  writeLines(text, temporary_path)
  
  if (!file.rename(temporary_path, path)) {
    stop("Could not publish text file: ", path)
  }
  
  invisible(path)
}

validation_n <- read_positive_integer_env(
  "SIM_EMPIRICAL_PRODUCTION_VALIDATION_N",
  100000L
)

selection_calibration_n <- read_positive_integer_env(
  "SIM_EMPIRICAL_SELECTION_CALIBRATION_N",
  50000L
)

empirical_seed_master <- 20260820L

# Set canary and production-cache directories.
empirical_results_dir <- file.path(
  results_dir,
  "empirical_residual_sensitivity"
)

empirical_canary_dir <- file.path(
  empirical_results_dir,
  "canary"
)

empirical_production_dir <- file.path(
  empirical_results_dir,
  "production"
)

empirical_cache_dir <- file.path(
  empirical_production_dir,
  "cache"
)

empirical_dgm_dir <- file.path(
  empirical_cache_dir,
  "dgm"
)

invisible(
  lapply(
    c(
      empirical_results_dir,
      empirical_production_dir,
      empirical_cache_dir,
      empirical_dgm_dir
    ),
    dir.create,
    showWarnings = FALSE,
    recursive = TRUE
  )
)

canary_pass_file <- file.path(
  empirical_canary_dir,
  "EMPIRICAL_RESIDUAL_CANARY_PASS.txt"
)

canary_qc_file <- file.path(
  empirical_canary_dir,
  "empirical_residual_canary_qc.csv"
)

cache_pass_file <- file.path(
  empirical_cache_dir,
  "EMPIRICAL_RESIDUAL_PRODUCTION_CACHE_PASS.txt"
)

# Prevent a previous pass file from surviving a failed rebuild.
if (file.exists(cache_pass_file)) {
  unlink(cache_pass_file)
}

# Require a completed and fully passing canary.
required_canary_files <- c(
  canary_pass_file,
  canary_qc_file
)

missing_canary_files <- required_canary_files[
  !file.exists(required_canary_files)
]

if (length(missing_canary_files) > 0L) {
  stop(
    "Required canary files are missing:\n- ",
    paste(missing_canary_files, collapse = "\n- "),
    "\nRun 23_run_empirical_residual_canary.R first."
  )
}

canary_qc <- read_csv(
  canary_qc_file,
  show_col_types = FALSE
)

required_canary_columns <- c(
  "check",
  "passed",
  "overall_canary_pass"
)

missing_canary_columns <- setdiff(
  required_canary_columns,
  names(canary_qc)
)

if (length(missing_canary_columns) > 0L) {
  stop(
    "Canary QC file is missing: ",
    paste(missing_canary_columns, collapse = ", "),
    "."
  )
}

if (
  nrow(canary_qc) == 0L ||
  !all_true(canary_qc$passed) ||
  !all_true(canary_qc$overall_canary_pass)
) {
  stop(
    "The empirical-residual canary QC file does not record a full pass."
  )
}

# Prepare and save the covariance-matched empirical calibration.
if (!file.exists(calibration_rds)) {
  stop(
    "Calibration object not found:\n",
    calibration_rds,
    "\nRun 01_calibrate_midus.R first."
  )
}

calibration <- readRDS(calibration_rds)

empirical_preparation <- prepare_empirical_residual_calibration(
  calibration
)

empirical_calibration <- empirical_preparation$calibration

empirical_calibration_file <- file.path(
  empirical_cache_dir,
  "empirical_residual_calibration.rds"
)

atomic_save_rds(
  empirical_calibration,
  empirical_calibration_file,
  compress = FALSE
)

atomic_write_csv(
  empirical_preparation$pool_summary,
  file.path(
    empirical_cache_dir,
    "empirical_residual_production_pool_diagnostics.csv"
  )
)

atomic_write_csv(
  empirical_preparation$marginal_diagnostics,
  file.path(
    empirical_cache_dir,
    "empirical_residual_production_marginal_diagnostics.csv"
  )
)

atomic_write_csv(
  empirical_sensitivity_grid,
  file.path(
    empirical_cache_dir,
    "empirical_residual_production_scenarios.csv"
  )
)

# Return the matched primary DGM filename.
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

# Return the empirical-residual DGM filename.
empirical_cache_filename <- function(scenario) {
  file.path(
    empirical_dgm_dir,
    sprintf(
      "empirical_dgm_%02d_%s.rds",
      scenario$sensitivity_scenario_id,
      scenario$scenario_key
    )
  )
}

n_scenarios <- nrow(empirical_sensitivity_grid)

dgm_parameter_comparison_list <- vector(
  "list",
  n_scenarios
)

dgm_validation_list <- vector(
  "list",
  n_scenarios
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

# Build and validate each empirical-residual DGM.
for (i in seq_len(n_scenarios)) {
  scenario <- empirical_sensitivity_grid[
    i,
    ,
    drop = FALSE
  ]
  
  cat(
    "Building empirical production DGM ",
    i,
    " of ",
    n_scenarios,
    ": ",
    scenario$scenario_label,
    "\n",
    sep = ""
  )
  
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
      "Production DGM was not built with empirical marker errors."
    )
  }
  
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
  
  empirical_dgm_file <- empirical_cache_filename(
    scenario
  )
  
  atomic_save_rds(
    empirical_dgm,
    empirical_dgm_file,
    compress = FALSE
  )
  
  dgm_parameter_comparison_list[[i]] <- bind_cols(
    scenario,
    data.frame(
      primary_cache_file =
        basename(primary_cache_file),
      empirical_cache_file =
        basename(empirical_dgm_file),
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
      selection_seed = selection_seed,
      stringsAsFactors = FALSE
    )
  )
  
  validation_seed <- as.integer(
    empirical_seed_master +
      20000L +
      100L * scenario$primary_scenario_id
  )
  
  validation <- check_dgm(
    calibration = empirical_calibration,
    dgm_parameters = empirical_dgm,
    N_check = validation_n,
    seed = validation_seed
  )
  
  missing_validation_columns <- setdiff(
    required_validation_columns,
    names(validation)
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
  
  validation$selection_target_check <- (
    empirical_dgm$selection_parameters$
      expected_fraction_check
  )
  
  validation$dgm_marker_error <- (
    empirical_dgm$marker_error
  )
  
  validation$validation_seed <- validation_seed
  
  dgm_validation_list[[i]] <- bind_cols(
    scenario[
      rep(1L, nrow(validation)),
      ,
      drop = FALSE
    ],
    validation
  )
}

dgm_parameter_comparison <- bind_rows(
  dgm_parameter_comparison_list
)

atomic_write_csv(
  dgm_parameter_comparison,
  file.path(
    empirical_cache_dir,
    "empirical_vs_primary_production_dgm_parameter_check.csv"
  )
)

# Apply the same DGM tolerances used by the canary.
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

atomic_write_csv(
  dgm_validation,
  file.path(
    empirical_cache_dir,
    "empirical_residual_production_dgm_validation.csv"
  )
)

# Make the final production-cache decision.
expected_primary_ids <- c(
  60L,
  168L,
  100L,
  208L
)

expected_dgm_files <- vapply(
  seq_len(n_scenarios),
  function(i) {
    empirical_cache_filename(
      empirical_sensitivity_grid[
        i,
        ,
        drop = FALSE
      ]
    )
  },
  character(1)
)

cache_qc <- data.frame(
  check = c(
    "Completed canary records a full pass",
    "Four intended production scenarios are present",
    "Primary scenario IDs match the planned scenarios",
    "Residual-pool means match zero",
    "Residual-pool covariance matches Sigma_marker",
    "Exposure parameters match the primary DGMs",
    "Outcome parameters match the primary DGMs",
    "Design values match the primary DGMs",
    "Only marker-error generation changes",
    "Phase-2 intercepts are recalibrated",
    "All large-sample DGM checks pass",
    "All empirical DGM cache files exist",
    "Empirical calibration cache exists"
  ),
  passed = c(
    all_true(
      canary_qc$passed &
        canary_qc$overall_canary_pass
    ),
    n_scenarios == 4L,
    identical(
      empirical_sensitivity_grid$primary_scenario_id,
      expected_primary_ids
    ),
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
    all(file.exists(expected_dgm_files)),
    file.exists(empirical_calibration_file)
  ),
  stringsAsFactors = FALSE
)

cache_qc$overall_cache_pass <- all_true(
  cache_qc$passed
)

atomic_write_csv(
  cache_qc,
  file.path(
    empirical_cache_dir,
    "empirical_residual_production_cache_qc.csv"
  )
)

if (!all_true(cache_qc$passed)) {
  failed_checks <- cache_qc$check[
    !cache_qc$passed
  ]
  
  stop(
    "Empirical residual production cache failed. Failed checks:\n- ",
    paste(failed_checks, collapse = "\n- "),
    "\nSee outputs in: ",
    empirical_cache_dir
  )
}

atomic_write_lines(
  c(
    "EMPIRICAL RESIDUAL PRODUCTION CACHE: PASS",
    paste("Completed:", Sys.time()),
    paste("Production scenarios:", n_scenarios),
    paste(
      "Validation sample size per DGM:",
      validation_n
    ),
    paste(
      "Selection calibration size per DGM:",
      selection_calibration_n
    ),
    paste(
      "Empirical seed master:",
      empirical_seed_master
    )
  ),
  cache_pass_file
)

cat(
  "\nEmpirical residual production cache passed.\n",
  "Cache directory: ",
  empirical_cache_dir,
  "\n",
  sep = ""
)
