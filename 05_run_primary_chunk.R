# 05_run_primary_chunk.R
#
# Runs one chunk of the primary MIDUS two-phase simulation. This script is
# submitted as a Slurm array worker after the DGM cache has been created.
#
# Each array task:
# - identifies one of the 324 simulation scenarios;
# - runs one chunk of 100 repetitions;
# - applies all six analysis methods;
# - records method runtimes and random-number states;
# - saves a checkpoint every five repetitions; and
# - writes one final RDS file for the completed chunk.
#
# The 324 scenarios and 20 chunks per scenario produce 6,480 array tasks.
#
# Reproducibility uses one L'Ecuyer-CMRG stream per scenario and two
# non-overlapping substreams per repetition: one for data generation and one
# for stochastic analysis methods. The random-number mapping depends only on
# the scenario and repetition, not the order in which Slurm runs the tasks.

# Load shared settings, DGM functions, and analysis methods
source("00_config.R")
source("02_dgm.R")
source("03_methods.R")

library(dplyr)

# 1. Identify the Slurm array task

# Use the Slurm task ID on Torch or SIM_ARRAY_TASK_ID for a local test
task_id_text <- Sys.getenv(
  "SLURM_ARRAY_TASK_ID",
  unset = Sys.getenv("SIM_ARRAY_TASK_ID", unset = "")
)

if (!nzchar(task_id_text)) {
  stop(
    "No array task id found. On Torch use SLURM_ARRAY_TASK_ID. ",
    "For a local test, set SIM_ARRAY_TASK_ID."
  )
}

task_id <- as.integer(task_id_text)

# Confirm that the task ID falls within the production array
if (
  is.na(task_id) ||
  task_id < 1L ||
  task_id > n_primary_array_tasks
) {
  stop(
    "Array task id must be between 1 and ",
    n_primary_array_tasks,
    "."
  )
}

# Map the task ID to one scenario and one chunk
scenario_id <- ((task_id - 1L) %/% chunks_per_scenario) + 1L
chunk_id <- ((task_id - 1L) %% chunks_per_scenario) + 1L

# Determine the repetitions assigned to the chunk
rep_start <- ((chunk_id - 1L) * reps_per_chunk) + 1L
rep_end <- min(chunk_id * reps_per_chunk, nsim_primary)

scenario_row <- primary_grid[
  primary_grid$scenario_id == scenario_id,
  ,
  drop = FALSE
]

if (nrow(scenario_row) != 1L) {
  stop("Could not uniquely identify scenario ", scenario_id, ".")
}

# 2. Locate the calibration, cache, and output files

# Reproduce the cache naming rule used in 04_build_dgm_cache.R
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

cache_file <- cache_filename(
  phase2_fraction = scenario_row$phase2_fraction,
  r2_a_marker = scenario_row$r2_a_marker,
  r2_y_marker = scenario_row$r2_y_marker,
  theta = scenario_row$theta
)

# Confirm that the required calibration files are available
if (!file.exists(calibration_rds)) {
  stop("Calibration object not found: ", calibration_rds)
}

if (!file.exists(cache_file)) {
  stop(
    "DGM cache file not found: ",
    cache_file,
    "\nRun 04_build_dgm_cache.R first."
  )
}

calibration <- readRDS(calibration_rds)
dgm <- readRDS(cache_file)

# Create a separate output directory for each scenario
scenario_chunk_dir <- file.path(
  chunk_dir,
  sprintf("scenario_%03d", scenario_id)
)

dir.create(
  scenario_chunk_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# Define the final output and temporary checkpoint files
final_file <- file.path(
  scenario_chunk_dir,
  sprintf(
    "scenario_%03d_chunk_%02d.rds",
    scenario_id,
    chunk_id
  )
)

checkpoint_file <- file.path(
  scenario_chunk_dir,
  sprintf(
    "scenario_%03d_chunk_%02d_checkpoint.rds",
    scenario_id,
    chunk_id
  )
)

# Do not rerun a chunk that has already been completed
if (file.exists(final_file)) {
  cat(
    "Final chunk already exists; exiting:\n",
    final_file,
    "\n",
    sep = ""
  )
  
  quit(save = "no", status = 0)
}

# 3. Define reproducible random-number streams

# Create one independent L'Ecuyer-CMRG stream for each scenario
get_scenario_stream <- function(scenario_id) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(rng_seed_master)
  
  state <- .Random.seed
  
  if (scenario_id > 1L) {
    for (s in seq_len(scenario_id - 1L)) {
      state <- parallel::nextRNGStream(state)
    }
  }
  
  state
}

# Advance to the data-generation substream for a specified repetition
advance_to_repetition <- function(scenario_state, repetition) {
  state <- scenario_state
  n_advance <- 2L * (repetition - 1L)
  
  if (n_advance > 0L) {
    for (j in seq_len(n_advance)) {
      state <- parallel::nextRNGSubStream(state)
    }
  }
  
  state
}

# Convert an RNG state to text for storage in the output
state_to_string <- function(state) {
  paste(state, collapse = ",")
}

# 4. Time the six analysis methods

# Run one method, record its runtime, and retain a failure row if it errors
run_timed <- function(method_name, expr) {
  start_time <- proc.time()[["elapsed"]]
  
  result <- tryCatch(
    eval.parent(substitute(expr)),
    error = function(e) {
      failed_result(method_name, conditionMessage(e))
    }
  )
  
  end_time <- proc.time()[["elapsed"]]
  result$elapsed_seconds <- end_time - start_time
  
  result
}

# Run and time all six primary methods
run_all_methods_timed <- function(dat, marker_names) {
  results <- list(
    run_timed(
      "Naive",
      fit_naive(dat)
    ),
    run_timed(
      "CCA",
      fit_cca(dat, marker_names)
    ),
    run_timed(
      "FCS-MI",
      fit_fcs_mi(
        dat = dat,
        marker_names = marker_names,
        nimp = nimp_primary,
        maxit = mice_maxit_primary,
        method = mice_method_primary
      )
    ),
    run_timed(
      "JM-MI",
      fit_joint_mi(
        dat = dat,
        marker_names = marker_names,
        nimp = nimp_primary,
        nburn = jomo_nburn_primary,
        nbetween = jomo_nbetween_primary
      )
    ),
    run_timed(
      "IPW",
      fit_ipw(dat, marker_names)
    ),
    run_timed(
      "AIPW",
      fit_aipw(
        dat = dat,
        marker_names = marker_names,
        selection_include_y = TRUE,
        marker_include_y = TRUE,
        method_label = "AIPW"
      )
    )
  )
  
  bind_rows(results)
}

# 5. Resume an incomplete chunk

# Initialize storage for method results and RNG states
results_list <- list()
rng_list <- list()
next_rep <- rep_start

# Restore saved results when a checkpoint is available
if (file.exists(checkpoint_file)) {
  checkpoint <- readRDS(checkpoint_file)
  
  results_list <- checkpoint$results_list
  rng_list <- checkpoint$rng_list
  next_rep <- checkpoint$next_rep
  
  cat(
    "Resuming scenario ",
    scenario_id,
    ", chunk ",
    chunk_id,
    " at repetition ",
    next_rep,
    ".\n",
    sep = ""
  )
}

# 6. Run the assigned repetitions

# Initialize the stream at the next repetition to be completed
scenario_state <- get_scenario_stream(scenario_id)

state <- advance_to_repetition(
  scenario_state,
  next_rep
)

# Print the scenario and production settings
cat(
  "\n============================================================\n",
  "Task ID: ", task_id, "\n",
  "Scenario: ", scenario_id, " / 324\n",
  "Chunk: ", chunk_id, " / ", chunks_per_scenario, "\n",
  "Repetitions: ", rep_start, "-", rep_end, "\n",
  "N: ", scenario_row$N, "\n",
  "Phase-2 fraction: ", scenario_row$phase2_fraction, "\n",
  "R2 A~M|X: ", scenario_row$r2_a_marker, "\n",
  "R2 Y~M|A,X: ", scenario_row$r2_y_marker, "\n",
  "Theta: ", scenario_row$theta, "\n",
  "JM-MI nburn: ", jomo_nburn_primary, "\n",
  "JM-MI nbetween: ", jomo_nbetween_primary, "\n",
  "============================================================\n",
  sep = ""
)

if (next_rep <= rep_end) {
  for (repetition in next_rep:rep_end) {
    # Assign separate substreams to data generation and stochastic methods
    data_state <- state
    method_state <- parallel::nextRNGSubStream(data_state)
    
    # Reserve the next pair before running the current repetition
    state <- parallel::nextRNGSubStream(method_state)
    
    # Generate one observed two-phase dataset
    .Random.seed <- data_state
    
    simulated_data <- generate_two_phase_data(
      N = scenario_row$N,
      calibration = calibration,
      dgm_parameters = dgm
    )
    
    dat <- simulated_data$observed
    
    # Run the six methods using an independent random-number substream
    .Random.seed <- method_state
    
    method_results <- run_all_methods_timed(
      dat = dat,
      marker_names = calibration$marker_names
    )
    
    # Attach the scenario settings and realized sampling diagnostics
    method_results <- method_results %>%
      mutate(
        task_id = task_id,
        scenario_id = scenario_id,
        chunk_id = chunk_id,
        repetition = repetition,
        N = scenario_row$N,
        target_phase2_fraction = scenario_row$phase2_fraction,
        realized_phase2_fraction = mean(dat$phase2),
        n_phase2 = sum(dat$phase2),
        r2_a_marker = scenario_row$r2_a_marker,
        r2_y_marker = scenario_row$r2_y_marker,
        theta_true = scenario_row$theta,
        min_pi_true = min(dat$pi_true),
        p01_pi_true = unname(quantile(dat$pi_true, 0.01)),
        median_pi_true = median(dat$pi_true),
        p99_pi_true = unname(quantile(dat$pi_true, 0.99)),
        max_pi_true = max(dat$pi_true)
      )
    
    list_index <- repetition - rep_start + 1L
    results_list[[list_index]] <- method_results
    
    # Save both RNG states so any repetition can be reproduced exactly
    rng_list[[list_index]] <- data.frame(
      scenario_id = scenario_id,
      repetition = repetition,
      data_rng_state = state_to_string(data_state),
      method_rng_state = state_to_string(method_state),
      stringsAsFactors = FALSE
    )
    
    # Save progress every five repetitions and at the end of the chunk
    if (repetition %% 5L == 0L || repetition == rep_end) {
      saveRDS(
        list(
          task_id = task_id,
          scenario_id = scenario_id,
          chunk_id = chunk_id,
          rep_start = rep_start,
          rep_end = rep_end,
          next_rep = repetition + 1L,
          results_list = results_list,
          rng_list = rng_list
        ),
        checkpoint_file
      )
      
      cat(
        "Scenario ",
        scenario_id,
        ", chunk ",
        chunk_id,
        ": completed repetition ",
        repetition,
        " / ",
        rep_end,
        "\n",
        sep = ""
      )
    }
  }
}

# 7. Finalize the completed chunk

# Combine repetition-level results and stored RNG states
results <- bind_rows(results_list)
rng_states <- bind_rows(rng_list)

expected_repetitions <- seq.int(rep_start, rep_end)

# Confirm that every assigned repetition appears in the final results
if (
  !setequal(
    unique(results$repetition),
    expected_repetitions
  )
) {
  stop("Chunk did not contain the expected repetition set.")
}

# Record the scenario and method settings used for this chunk
metadata <- list(
  task_id = task_id,
  scenario_id = scenario_id,
  chunk_id = chunk_id,
  rep_start = rep_start,
  rep_end = rep_end,
  N = scenario_row$N,
  phase2_fraction = scenario_row$phase2_fraction,
  r2_a_marker = scenario_row$r2_a_marker,
  r2_y_marker = scenario_row$r2_y_marker,
  theta = scenario_row$theta,
  nimp = nimp_primary,
  mice_method = mice_method_primary,
  mice_maxit = mice_maxit_primary,
  jomo_nburn = jomo_nburn_primary,
  jomo_nbetween = jomo_nbetween_primary,
  created = as.character(Sys.time())
)

# Save the completed results, metadata, and RNG states
saveRDS(
  list(
    metadata = metadata,
    results = results,
    rng_states = rng_states
  ),
  final_file
)

# Remove the checkpoint after the final file has been written
if (file.exists(checkpoint_file)) {
  unlink(checkpoint_file)
}

cat(
  "\nChunk complete:\n",
  final_file,
  "\n",
  sep = ""
)
