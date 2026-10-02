# 04_build_dgm_cache.R
#
# Builds and saves the fixed DGM parameters used by the production simulation.
# Run this script once on NYU Torch after completing MIDUS calibration and
# before submitting the Slurm production array.
#
# The primary design contains 324 scenarios, but phase-1 sample size does not
# affect the DGM parameters. The three sample-size levels therefore share the
# same 108 combinations of:
# - phase-2 sampling fraction;
# - marker contribution to the exposure;
# - marker contribution to the outcome; and
# - true exposure effect.
#
# Each unique parameter object is calibrated once and saved as an RDS file.
# Existing cache files are reused so an interrupted run can be restarted.

# Load shared settings and DGM functions
source("00_config.R")
source("02_dgm.R")

library(dplyr)
library(readr)

# Confirm that the MIDUS calibration object is available
if (!file.exists(calibration_rds)) {
  stop(
    "Calibration object not found at:\n",
    calibration_rds
  )
}

calibration <- readRDS(calibration_rds)

# Save the complete 324-scenario production manifest
write_csv(
  primary_grid,
  file.path(calibration_dir, "primary_scenario_grid.csv")
)

# Remove phase-1 sample size to identify the 108 unique parameter combinations
unique_dgm_grid <- primary_grid %>%
  distinct(
    phase2_fraction,
    r2_a_marker,
    r2_y_marker,
    theta
  ) %>%
  arrange(
    phase2_fraction,
    r2_a_marker,
    r2_y_marker,
    theta
  ) %>%
  mutate(dgm_id = row_number())

# Confirm the expected number of unique DGMs
stopifnot(nrow(unique_dgm_grid) == 108L)

# Construct a unique cache filename from the four DGM factors
cache_filename <- function(
    phase2_fraction,
    r2_a_marker,
    r2_y_marker,
    theta
) {
  file.path(
    dgm_cache_dir,
    sprintf(
      "dgm_f%03d_rA%03d_rY%03d_th%03d.rds",
      round(100 * phase2_fraction),
      round(100 * r2_a_marker),
      round(100 * r2_y_marker),
      round(100 * theta)
    )
  )
}

# Store calibration checks for each unique DGM
diagnostic_rows <- vector(
  "list",
  nrow(unique_dgm_grid)
)

for (i in seq_len(nrow(unique_dgm_grid))) {
  scenario <- unique_dgm_grid[i, , drop = FALSE]
  
  cache_file <- cache_filename(
    phase2_fraction = scenario$phase2_fraction,
    r2_a_marker = scenario$r2_a_marker,
    r2_y_marker = scenario$r2_y_marker,
    theta = scenario$theta
  )
  
  cat(
    "Building DGM ",
    i,
    " / ",
    nrow(unique_dgm_grid),
    ": f2=",
    scenario$phase2_fraction,
    ", rA=",
    scenario$r2_a_marker,
    ", rY=",
    scenario$r2_y_marker,
    ", theta=",
    scenario$theta,
    "\n",
    sep = ""
  )
  
  if (file.exists(cache_file)) {
    # Reuse an existing parameter object when restarting the script
    dgm <- readRDS(cache_file)
  } else {
    # Use a fixed seed for reproducible selection-intercept calibration
    selection_seed <- as.integer(
      500000L + 10000L * i
    )
    
    dgm <- build_dgm_parameters(
      calibration = calibration,
      r2_a_marker = scenario$r2_a_marker,
      r2_y_marker = scenario$r2_y_marker,
      theta = scenario$theta,
      phase2_fraction = scenario$phase2_fraction,
      selection = "mar",
      selection_calibration_n = 50000L,
      selection_seed = selection_seed,
      marker_error = "mvn"
    )
    
    saveRDS(dgm, cache_file)
  }
  
  # Record whether the calibrated parameters reproduce their targets
  diagnostic_rows[[i]] <- data.frame(
    dgm_id = scenario$dgm_id,
    phase2_fraction = scenario$phase2_fraction,
    r2_a_marker = scenario$r2_a_marker,
    r2_y_marker = scenario$r2_y_marker,
    theta = scenario$theta,
    var_A_check = dgm$exposure$var_A_check,
    partial_r2_A_check = dgm$exposure$partial_r2_check,
    var_Y_check = dgm$outcome$var_Y_check,
    partial_r2_Y_check = dgm$outcome$partial_r2_check,
    expected_phase2_fraction_check =
      dgm$selection_parameters$expected_fraction_check,
    cache_file = basename(cache_file)
  )
}

dgm_diagnostics <- bind_rows(diagnostic_rows)

# Save the unique parameter grid and its calibration checks
write_csv(
  unique_dgm_grid,
  file.path(calibration_dir, "unique_dgm_grid.csv")
)

write_csv(
  dgm_diagnostics,
  file.path(calibration_dir, "dgm_analytic_checks.csv")
)

# Report completion and cache location
cat(
  "\nDGM cache complete.\n",
  "Unique DGM objects: ",
  nrow(unique_dgm_grid),
  "\nCache directory: ",
  dgm_cache_dir,
  "\n",
  sep = ""
)