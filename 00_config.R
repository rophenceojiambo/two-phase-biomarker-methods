# 00_config.R
#
# Defines the shared settings used throughout the MIDUS two-phase simulation.
# This includes project paths, analysis variables, simulation scenarios,
# imputation settings, Slurm chunking, random-number settings, numerical
# tolerances, and reporting thresholds.
#
# Source this file before running the remaining simulation scripts.
#
# The script uses the current working directory for local runs. On NYU Torch,
# set the project directory before running the simulation:
# export SIM_PROJECT_DIR=/scratch/$USER/midus_two_phase_sim
#
# Primary production design:
# - 324 data-generating scenarios
# - 2,000 repetitions per scenario
# - 20 imputations per incomplete dataset
# - jomo: 1,000 burn-in iterations and 1,000 iterations between imputations

# Keep character variables as character rather than converting them to factors
options(stringsAsFactors = FALSE)

# 1. Project paths

# Use SIM_PROJECT_DIR when supplied; otherwise use the working directory
project_dir <- Sys.getenv(
  "SIM_PROJECT_DIR",
  unset = normalizePath(".", winslash = "/", mustWork = FALSE)
)

# Allow the analytic dataset path to be overridden when needed
analytic_rds <- Sys.getenv(
  "SIM_ANALYTIC_RDS",
  unset = file.path(project_dir, "data", "MIDUS_discrimination_analysis.rds")
)

# Define directories for intermediate and final outputs
calibration_dir <- file.path(project_dir, "calibration")
results_dir <- file.path(project_dir, "results")
chunk_dir <- file.path(results_dir, "chunks")
combined_dir <- file.path(results_dir, "combined")
summary_dir <- file.path(results_dir, "summary")
figure_dir <- file.path(results_dir, "figures")
table_dir <- file.path(results_dir, "tables")
log_dir <- file.path(project_dir, "logs")

# Create directories that do not already exist
dirs <- c(
  calibration_dir, results_dir, chunk_dir, combined_dir,
  summary_dir, figure_dir, table_dir, log_dir
)

invisible(lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

# Define the calibration file and cached DGM directory
calibration_rds <- Sys.getenv(
  "SIM_CALIBRATION_RDS",
  unset = file.path(calibration_dir, "midus_calibration.rds")
)

dgm_cache_dir <- file.path(calibration_dir, "dgm_cache")
dir.create(dgm_cache_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Required packages

# Stop before simulation if any required package is unavailable
required_packages <- c(
  "dplyr", "tidyr", "purrr", "readr", "mvtnorm", "mice", "jomo",
  "broom", "ggplot2", "rsimsum", "generics", "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before continuing: ",
    paste(missing_packages, collapse = ", ")
  )
}

# 3. Analysis variables

# Standardized discrimination measure used as the exposure
exposure_name <- "discrimination_STD"

# Outcome used only to calibrate the direction of empirical coefficients
calibration_outcome <- "dunedinpace_STD"

# Standardized biomarkers observed only in the phase-2 sample
marker_names <- c(
  "log2_cd19_STD", "log2_cd3d_STD", "log2_cd3e_STD", "log2_cd4_STD",
  "log2_cd8a_STD", "log2_cd14_STD", "log2_fcgr3a_STD", "log2_ncam1_STD"
)

# Phase-1 covariates included in the analysis models
x_formula <- ~ age_STD + sex + race_eth

# 4. Primary factorial design

# Vary the phase-2 sampling fraction and phase-1 sample size
phase2_fractions <- c(0.25, 0.50, 0.65)
phase1_sizes <- c(500L, 800L, 1500L)

# Vary marker contributions to the outcome and exposure
r2_y_marker_levels <- c(0.00, 0.02, 0.05, 0.10)
r2_a_marker_levels <- c(0.01, 0.05, 0.10)

# Vary the true exposure effect
theta_levels <- c(0.00, 0.15, 0.30)

# Cross all design factors to create 3 × 3 × 4 × 3 × 3 = 324 scenarios
primary_grid <- expand.grid(
  phase2_fraction = phase2_fractions,
  N = phase1_sizes,
  r2_y_marker = r2_y_marker_levels,
  r2_a_marker = r2_a_marker_levels,
  theta = theta_levels,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Assign a unique identifier and retain a consistent column order
primary_grid$scenario_id <- seq_len(nrow(primary_grid))

primary_grid <- primary_grid[, c(
  "scenario_id", "N", "phase2_fraction",
  "r2_a_marker", "r2_y_marker", "theta"
)]

# Confirm that the complete factorial design contains 324 scenarios
stopifnot(nrow(primary_grid) == 324L)

# 5. Monte Carlo and imputation settings

# Number of simulation repetitions per scenario
nsim_primary <- 2000L

# Number of imputations per incomplete dataset
nimp_primary <- 20L

# MICE settings for normal linear imputation
mice_method_primary <- "norm"
mice_maxit_primary <- 10L

# Conservative JM-MI settings selected for final production
jomo_nburn_primary <- 1000L
jomo_nbetween_primary <- 1000L

# 6. Slurm chunking

# Run 100 repetitions per task, giving 20 chunks per scenario
reps_per_chunk <- 100L
chunks_per_scenario <- as.integer(ceiling(nsim_primary / reps_per_chunk))

# The full production run requires 324 × 20 = 6,480 array tasks
n_primary_array_tasks <- nrow(primary_grid) * chunks_per_scenario

# Confirm the expected chunk and task counts
stopifnot(chunks_per_scenario == 20L)
stopifnot(n_primary_array_tasks == 6480L)

# 7. Random-number settings

# Used to initialize reproducible L'Ecuyer-CMRG streams across Slurm tasks
rng_seed_master <- 20260810L

# 8. Numerical tolerances

# Safeguards used in probability calculations and matrix checks
probability_floor <- 1e-8
matrix_tolerance <- 1e-10

# 9. Reporting settings

# Nominal confidence level and corresponding significance level
nominal_level <- 0.95
alpha_level <- 0.05

# These thresholds flag concerning results but do not affect estimation
deterioration_thresholds <- list(
  abs_bias = 0.02,
  se_ratio_low = 0.90,
  se_ratio_high = 1.10,
  coverage_low = 0.925,
  coverage_high = 0.975,
  type1_high = 0.075,
  failure_rate = 0.01
)
