################################################################################
# 30_build_phase2_mcar_production_cache.R
#
# Builds and validates the four Phase-2 MCAR production DGMs after the canary.
################################################################################

source("00_config.R")
source("02_dgm.R")
source("03_methods.R")
source("28_phase2_mcar_sensitivity_helpers.R")

library(dplyr)
library(readr)


# ------------------------------------------------------------------------------
# 1. Settings, paths, and canary gate
# ------------------------------------------------------------------------------

validation_n <- read_mcar_positive_integer_env(
  "SIM_MCAR_PRODUCTION_VALIDATION_N",
  100000L
)

if (validation_n < 100000L) {
  stop("SIM_MCAR_PRODUCTION_VALIDATION_N must be at least 100000.")
}

mcar_results_dir <- file.path(results_dir, "phase2_mcar_sensitivity")
mcar_canary_dir <- file.path(mcar_results_dir, "canary")
mcar_production_dir <- file.path(mcar_results_dir, "production")
mcar_cache_dir <- file.path(mcar_production_dir, "cache")
mcar_dgm_dir <- file.path(mcar_cache_dir, "dgm")

invisible(
  lapply(
    c(mcar_results_dir, mcar_production_dir, mcar_cache_dir, mcar_dgm_dir),
    dir.create,
    showWarnings = FALSE,
    recursive = TRUE
  )
)

canary_pass_file <- file.path(
  mcar_canary_dir,
  "PHASE2_MCAR_CANARY_PASS.txt"
)

canary_qc_file <- file.path(
  mcar_canary_dir,
  "phase2_mcar_canary_qc.csv"
)

cache_pass_file <- file.path(
  mcar_cache_dir,
  "PHASE2_MCAR_PRODUCTION_CACHE_PASS.txt"
)

if (file.exists(cache_pass_file)) {
  unlink(cache_pass_file)
}

if (!file.exists(canary_pass_file) || !file.exists(canary_qc_file)) {
  stop("Completed Phase-2 MCAR canary outputs were not found.")
}

if (!any(grepl(
  "MCAR CANARY: PASS",
  readLines(canary_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The Phase-2 MCAR canary PASS marker is invalid.")
}

canary_qc <- read_csv(canary_qc_file, show_col_types = FALSE)

if (
  nrow(canary_qc) == 0L ||
    !mcar_all_true(canary_qc$passed) ||
    !mcar_all_true(canary_qc$overall_canary_pass)
) {
  stop("The Phase-2 MCAR canary QC file does not record a full pass.")
}

if (!file.exists(calibration_rds)) {
  stop("Calibration object not found: ", calibration_rds)
}

calibration <- readRDS(calibration_rds)

mcar_atomic_write_csv(
  mcar_sensitivity_grid,
  file.path(mcar_cache_dir, "phase2_mcar_production_scenarios.csv")
)


# ------------------------------------------------------------------------------
# 2. Build, compare, and validate all four production DGMs
# ------------------------------------------------------------------------------

parameter_check_list <- vector("list", nrow(mcar_sensitivity_grid))
dgm_validation_list <- vector("list", nrow(mcar_sensitivity_grid))
independence_validation_list <- vector("list", nrow(mcar_sensitivity_grid))

for (i in seq_len(nrow(mcar_sensitivity_grid))) {

  row <- mcar_sensitivity_grid[i, , drop = FALSE]

  cat(
    "Building MCAR production DGM ",
    i,
    " / ",
    nrow(mcar_sensitivity_grid),
    ": ",
    row$scenario_label,
    "\n",
    sep = ""
  )

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
    20260821L + 30000L + 100L * row$primary_scenario_id
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
}

dgm_parameter_check <- bind_rows(parameter_check_list)
dgm_validation <- bind_rows(dgm_validation_list)
independence_validation <- bind_rows(independence_validation_list)

mcar_atomic_write_csv(
  dgm_parameter_check,
  file.path(
    mcar_cache_dir,
    "mcar_vs_primary_production_dgm_parameter_check.csv"
  )
)

mcar_atomic_write_csv(
  dgm_validation,
  file.path(mcar_cache_dir, "phase2_mcar_production_dgm_validation.csv")
)

mcar_atomic_write_csv(
  independence_validation,
  file.path(
    mcar_cache_dir,
    "phase2_mcar_production_independence_validation.csv"
  )
)


# ------------------------------------------------------------------------------
# 3. Final production-cache gate
# ------------------------------------------------------------------------------

expected_primary_ids <- c(60L, 168L, 100L, 208L)

expected_mcar_files <- vapply(
  seq_len(nrow(mcar_sensitivity_grid)),
  function(i) {
    mcar_dgm_cache_filename(
      mcar_sensitivity_grid[i, , drop = FALSE],
      mcar_dgm_dir
    )
  },
  character(1)
)

cache_qc <- data.frame(
  check = c(
    "Completed Phase-2 MCAR canary records a full pass",
    "Four intended production scenarios are present",
    "Primary scenario IDs are 60, 168, 100, and 208",
    "Exposure parameters exactly match the primary DGMs",
    "Outcome parameters exactly match the primary DGMs",
    "Common design values exactly match the primary DGMs",
    "Only Phase-2 selection changes from MAR to MCAR",
    "Marker-error generation remains multivariate normal",
    "MCAR selection target is exact without intercept calibration",
    "All large-sample DGM checks pass",
    "All constant-probability and independence checks pass",
    "All four MCAR DGM cache files exist"
  ),
  passed = c(
    mcar_all_true(canary_qc$passed) &
      mcar_all_true(canary_qc$overall_canary_pass),
    nrow(mcar_sensitivity_grid) == 4L,
    identical(
      mcar_sensitivity_grid$primary_scenario_id,
      expected_primary_ids
    ),
    mcar_all_true(dgm_parameter_check$exposure_parameters_identical),
    mcar_all_true(dgm_parameter_check$outcome_parameters_identical),
    mcar_all_true(dgm_parameter_check$common_design_values_identical),
    mcar_all_true(dgm_parameter_check$selection_changed_as_intended),
    mcar_all_true(dgm_parameter_check$marker_error_unchanged),
    mcar_all_true(dgm_parameter_check$mcar_intercept_not_required) &
      mcar_all_true(dgm_parameter_check$mcar_expected_fraction_exact),
    mcar_all_true(dgm_validation$pass_all),
    mcar_all_true(independence_validation$pass_all),
    mcar_all_true(file.exists(expected_mcar_files))
  ),
  stringsAsFactors = FALSE
)

cache_qc$overall_cache_pass <- mcar_all_true(cache_qc$passed)

mcar_atomic_write_csv(
  cache_qc,
  file.path(mcar_cache_dir, "phase2_mcar_production_cache_qc.csv")
)

if (!mcar_all_true(cache_qc$passed)) {

  failed_checks <- cache_qc$check[!cache_qc$passed]

  stop(
    "Phase-2 MCAR production cache FAILED. Failed checks:\n- ",
    paste(failed_checks, collapse = "\n- "),
    "\nSee outputs in: ",
    mcar_cache_dir
  )
}

mcar_atomic_write_lines(
  c(
    "PHASE-2 MCAR PRODUCTION CACHE: PASS",
    paste("Completed:", Sys.time()),
    paste("Production scenarios:", nrow(mcar_sensitivity_grid)),
    paste("Validation sample size per DGM:", validation_n),
    paste("Primary RNG seed master retained:", rng_seed_master)
  ),
  cache_pass_file
)

cat(
  "\nPhase-2 MCAR production cache PASSED.\n",
  "Cache directory: ",
  mcar_cache_dir,
  "\n",
  sep = ""
)
