################################################################################
# 25_run_empirical_residual_production_chunk.R
#
# Slurm-array worker for the empirical residual-resampling sensitivity.
#
# Defaults:
#   - four scenarios
#   - 2,000 repetitions per scenario
#   - 10 repetitions per task
#   - 800 total array tasks
#
# Seeds depend only on the matched primary scenario and repetition, making the
# results independent of array-task order and concurrency.
################################################################################

source("00_config.R")
source("02_dgm.R")
source("03_methods.R")
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

atomic_save_rds <- function(object, path) {
  temporary_path <- paste0(
    path,
    ".tmp_",
    Sys.getpid()
  )
  
  on.exit(
    unlink(temporary_path),
    add = TRUE
  )
  
  saveRDS(object, temporary_path)
  
  if (!file.rename(temporary_path, path)) {
    stop("Could not publish RDS file: ", path)
  }
  
  invisible(path)
}

nsim_empirical <- read_positive_integer_env(
  "SIM_EMPIRICAL_NSIM",
  2000L
)

reps_per_empirical_chunk <- read_positive_integer_env(
  "SIM_EMPIRICAL_REPS_PER_CHUNK",
  10L
)

chunks_per_empirical_scenario <- as.integer(
  ceiling(
    nsim_empirical /
      reps_per_empirical_chunk
  )
)

n_scenarios <- nrow(empirical_sensitivity_grid)

n_empirical_array_tasks <- (
  n_scenarios *
    chunks_per_empirical_scenario
)

# Map the array task to one scenario and chunk.
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

task_id <- suppressWarnings(
  as.integer(task_id_text)
)

if (
  is.na(task_id) ||
  task_id < 1L ||
  task_id > n_empirical_array_tasks
) {
  stop(
    "Array task ID must be between 1 and ",
    n_empirical_array_tasks,
    "; received '",
    task_id_text,
    "'."
  )
}

sensitivity_scenario_id <- (
  (task_id - 1L) %/%
    chunks_per_empirical_scenario
) + 1L

chunk_id <- (
  (task_id - 1L) %%
    chunks_per_empirical_scenario
) + 1L

rep_start <- (
  (chunk_id - 1L) *
    reps_per_empirical_chunk
) + 1L

rep_end <- min(
  chunk_id * reps_per_empirical_chunk,
  nsim_empirical
)

scenario <- empirical_sensitivity_grid[
  empirical_sensitivity_grid$sensitivity_scenario_id ==
    sensitivity_scenario_id,
  ,
  drop = FALSE
]

if (nrow(scenario) != 1L) {
  stop(
    "Could not uniquely identify empirical scenario ",
    sensitivity_scenario_id,
    "."
  )
}

# Set production paths.
empirical_seed_master <- 20260820L

empirical_results_dir <- file.path(
  results_dir,
  "empirical_residual_sensitivity"
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

empirical_chunk_dir <- file.path(
  empirical_production_dir,
  "chunks"
)

scenario_chunk_dir <- file.path(
  empirical_chunk_dir,
  sprintf(
    "scenario_%02d_%s",
    sensitivity_scenario_id,
    scenario$scenario_key
  )
)

dir.create(
  scenario_chunk_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

cache_pass_file <- file.path(
  empirical_cache_dir,
  "EMPIRICAL_RESIDUAL_PRODUCTION_CACHE_PASS.txt"
)

cache_qc_file <- file.path(
  empirical_cache_dir,
  "empirical_residual_production_cache_qc.csv"
)

empirical_calibration_file <- file.path(
  empirical_cache_dir,
  "empirical_residual_calibration.rds"
)

empirical_dgm_file <- file.path(
  empirical_dgm_dir,
  sprintf(
    "empirical_dgm_%02d_%s.rds",
    sensitivity_scenario_id,
    scenario$scenario_key
  )
)

required_cache_files <- c(
  cache_pass_file,
  cache_qc_file,
  empirical_calibration_file,
  empirical_dgm_file
)

missing_cache_files <- required_cache_files[
  !file.exists(required_cache_files)
]

if (length(missing_cache_files) > 0L) {
  stop(
    "Required production-cache files are missing:\n- ",
    paste(missing_cache_files, collapse = "\n- "),
    "\nRun 24_build_empirical_residual_production_cache.R first."
  )
}

cache_qc <- read_csv(
  cache_qc_file,
  show_col_types = FALSE
)

if (
  !all(c("passed", "overall_cache_pass") %in%
       names(cache_qc)) ||
  nrow(cache_qc) == 0L ||
  !all_true(cache_qc$passed) ||
  !all_true(cache_qc$overall_cache_pass)
) {
  stop(
    "The empirical production-cache QC file does not record a full pass."
  )
}

empirical_calibration <- readRDS(
  empirical_calibration_file
)

empirical_dgm <- readRDS(
  empirical_dgm_file
)

if (!identical(
  empirical_dgm$marker_error,
  "empirical"
)) {
  stop("Cached DGM does not use empirical marker errors.")
}

if (!identical(
  empirical_calibration$marker_names,
  marker_names
)) {
  stop(
    "Cached empirical calibration has unexpected marker names."
  )
}

final_file <- file.path(
  scenario_chunk_dir,
  sprintf(
    "empirical_scenario_%02d_chunk_%03d.rds",
    sensitivity_scenario_id,
    chunk_id
  )
)

checkpoint_file <- file.path(
  scenario_chunk_dir,
  sprintf(
    "empirical_scenario_%02d_chunk_%03d_checkpoint.rds",
    sensitivity_scenario_id,
    chunk_id
  )
)

# Safely reuse an already completed matching chunk.
if (file.exists(final_file)) {
  existing_final <- tryCatch(
    readRDS(final_file),
    error = function(e) NULL
  )
  
  existing_final_matches <- (
    is.list(existing_final) &&
      all(c(
        "metadata",
        "results",
        "rng_states",
        "dataset_qc"
      ) %in% names(existing_final)) &&
      identical(
        existing_final$metadata$
          sensitivity_scenario_id,
        sensitivity_scenario_id
      ) &&
      identical(
        existing_final$metadata$chunk_id,
        chunk_id
      ) &&
      identical(
        existing_final$metadata$rep_start,
        rep_start
      ) &&
      identical(
        existing_final$metadata$rep_end,
        rep_end
      ) &&
      identical(
        existing_final$metadata$nsim_empirical,
        nsim_empirical
      ) &&
      identical(
        existing_final$metadata$
          reps_per_empirical_chunk,
        reps_per_empirical_chunk
      ) &&
      identical(
        existing_final$metadata$marker_error,
        "empirical"
      )
  )
  
  if (!existing_final_matches) {
    stop(
      "Existing final chunk does not match this task:\n",
      final_file
    )
  }
  
  cat(
    "Final chunk already exists; exiting:\n",
    final_file,
    "\n"
  )
  
  quit(save = "no", status = 0L)
}

# Run one method while converting errors into failed-result rows.
run_timed <- function(method_name, expression) {
  start_time <- proc.time()[["elapsed"]]
  
  result <- tryCatch(
    eval.parent(substitute(expression)),
    error = function(e) {
      failed_result(
        method_name,
        conditionMessage(e)
      )
    }
  )
  
  result$elapsed_seconds <- (
    proc.time()[["elapsed"]] -
      start_time
  )
  
  result
}

# Run the six primary methods and retain per-method runtimes.
run_all_methods_timed <- function(
    dat,
    marker_names
) {
  start_time <- proc.time()[["elapsed"]]
  
  method_results <- list(
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
  
  method_results <- bind_rows(method_results)
  
  method_results$total_all_methods_elapsed_seconds <- (
    proc.time()[["elapsed"]] -
      start_time
  )
  
  method_results
}

# Initialize storage or resume from a matching checkpoint.
results_list <- list()
rng_list <- list()
dataset_qc_list <- list()
next_rep <- rep_start

if (file.exists(checkpoint_file)) {
  checkpoint <- readRDS(checkpoint_file)
  
  required_checkpoint_objects <- c(
    "task_id",
    "sensitivity_scenario_id",
    "primary_scenario_id",
    "chunk_id",
    "rep_start",
    "rep_end",
    "nsim_empirical",
    "reps_per_empirical_chunk",
    "next_rep",
    "results_list",
    "rng_list",
    "dataset_qc_list"
  )
  
  missing_checkpoint_objects <- setdiff(
    required_checkpoint_objects,
    names(checkpoint)
  )
  
  if (length(missing_checkpoint_objects) > 0L) {
    stop(
      "Checkpoint is missing: ",
      paste(
        missing_checkpoint_objects,
        collapse = ", "
      ),
      "."
    )
  }
  
  checkpoint_matches <- (
    identical(checkpoint$task_id, task_id) &&
      identical(
        checkpoint$sensitivity_scenario_id,
        sensitivity_scenario_id
      ) &&
      identical(
        checkpoint$primary_scenario_id,
        scenario$primary_scenario_id
      ) &&
      identical(checkpoint$chunk_id, chunk_id) &&
      identical(checkpoint$rep_start, rep_start) &&
      identical(checkpoint$rep_end, rep_end) &&
      identical(
        checkpoint$nsim_empirical,
        nsim_empirical
      ) &&
      identical(
        checkpoint$reps_per_empirical_chunk,
        reps_per_empirical_chunk
      )
  )
  
  if (!checkpoint_matches) {
    stop(
      "Checkpoint does not match this production task:\n",
      checkpoint_file
    )
  }
  
  next_rep <- as.integer(checkpoint$next_rep)
  
  if (
    length(next_rep) != 1L ||
    is.na(next_rep) ||
    next_rep < rep_start ||
    next_rep > rep_end + 1L
  ) {
    stop("Checkpoint contains an invalid next repetition.")
  }
  
  results_list <- checkpoint$results_list
  rng_list <- checkpoint$rng_list
  dataset_qc_list <- checkpoint$dataset_qc_list
  
  cat(
    "Resuming scenario ",
    sensitivity_scenario_id,
    ", chunk ",
    chunk_id,
    " at repetition ",
    next_rep,
    ".\n",
    sep = ""
  )
}

expected_methods <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

RNGkind("L'Ecuyer-CMRG")

cat(
  "\nEmpirical-residual production task\n",
  "Task: ", task_id, " of ",
  n_empirical_array_tasks, "\n",
  "Scenario: ", sensitivity_scenario_id,
  " of ", n_scenarios, "\n",
  "Scenario label: ", scenario$scenario_label, "\n",
  "Primary scenario: ",
  scenario$primary_scenario_id, "\n",
  "Chunk: ", chunk_id, " of ",
  chunks_per_empirical_scenario, "\n",
  "Repetitions: ", rep_start, "-", rep_end, "\n",
  "N: ", scenario$N, "\n",
  "Target Phase-2 fraction: ",
  scenario$phase2_fraction, "\n",
  "True exposure effect: ", scenario$theta, "\n",
  "Marker error: empirical\n",
  "MICE imputations: ", nimp_primary, "\n",
  "JM-MI burn-in: ", jomo_nburn_primary, "\n",
  "JM-MI spacing: ", jomo_nbetween_primary,
  "\n\n",
  sep = ""
)

# Run the repetitions assigned to this chunk.
if (next_rep <= rep_end) {
  for (repetition in seq.int(next_rep, rep_end)) {
    data_seed <- as.integer(
      empirical_seed_master +
        100000L +
        1000L * scenario$primary_scenario_id +
        2L * repetition
    )
    
    method_seed <- data_seed + 1L
    
    set.seed(data_seed)
    
    simulated_data <- generate_two_phase_data(
      N = scenario$N,
      calibration = empirical_calibration,
      dgm_parameters = empirical_dgm
    )
    
    analysis_data <- simulated_data$observed
    
    n_markers <- length(
      empirical_calibration$marker_names
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
    
    marker_blockwise <- all(
      n_markers_observed %in% c(0L, n_markers)
    )
    
    phase2_matches_markers <- identical(
      as.integer(n_markers_observed == n_markers),
      as.integer(analysis_data$phase2)
    )
    
    n_phase2 <- sum(analysis_data$phase2)
    
    sufficient_phase2_n <- (
      n_phase2 > n_markers + 10L
    )
    
    finite_probability_range <- (
      all(is.finite(analysis_data$pi_true)) &&
        all(analysis_data$pi_true > 0) &&
        all(analysis_data$pi_true < 1)
    )
    
    if (
      !marker_blockwise ||
      !phase2_matches_markers ||
      !sufficient_phase2_n ||
      !finite_probability_range
    ) {
      stop(
        "Dataset structural QC failed for scenario ",
        sensitivity_scenario_id,
        ", repetition ",
        repetition,
        "."
      )
    }
    
    list_index <- repetition - rep_start + 1L
    
    dataset_qc_list[[list_index]] <- data.frame(
      sensitivity_scenario_id =
        sensitivity_scenario_id,
      primary_scenario_id =
        scenario$primary_scenario_id,
      scenario_label = scenario$scenario_label,
      scenario_key = scenario$scenario_key,
      task_id = task_id,
      chunk_id = chunk_id,
      repetition = repetition,
      data_seed = data_seed,
      method_seed = method_seed,
      marker_blockwise = marker_blockwise,
      phase2_matches_markers =
        phase2_matches_markers,
      n_phase2 = n_phase2,
      realized_phase2_fraction =
        mean(analysis_data$phase2),
      min_pi_true = min(analysis_data$pi_true),
      max_pi_true = max(analysis_data$pi_true),
      sufficient_phase2_n =
        sufficient_phase2_n,
      finite_probability_range =
        finite_probability_range,
      stringsAsFactors = FALSE
    )
    
    set.seed(method_seed)
    
    method_results <- run_all_methods_timed(
      dat = analysis_data,
      marker_names =
        empirical_calibration$marker_names
    )
    
    method_results <- method_results %>%
      mutate(
        task_id = task_id,
        sensitivity_scenario_id =
          sensitivity_scenario_id,
        primary_scenario_id =
          scenario$primary_scenario_id,
        scenario_label = scenario$scenario_label,
        scenario_key = scenario$scenario_key,
        chunk_id = chunk_id,
        repetition = repetition,
        data_seed = data_seed,
        method_seed = method_seed,
        N = scenario$N,
        target_phase2_fraction =
          scenario$phase2_fraction,
        realized_phase2_fraction =
          mean(analysis_data$phase2),
        n_phase2 = n_phase2,
        r2_a_marker = scenario$r2_a_marker,
        r2_y_marker = scenario$r2_y_marker,
        theta_true = scenario$theta,
        min_pi_true = min(analysis_data$pi_true),
        p01_pi_true = unname(
          quantile(analysis_data$pi_true, 0.01)
        ),
        median_pi_true = median(
          analysis_data$pi_true
        ),
        p99_pi_true = unname(
          quantile(analysis_data$pi_true, 0.99)
        ),
        max_pi_true = max(analysis_data$pi_true),
        marker_error = empirical_dgm$marker_error,
        nimp = nimp_primary,
        mice_method = mice_method_primary,
        mice_maxit = mice_maxit_primary,
        jomo_nburn = jomo_nburn_primary,
        jomo_nbetween = jomo_nbetween_primary
      )
    
    results_list[[list_index]] <- method_results
    
    rng_list[[list_index]] <- data.frame(
      sensitivity_scenario_id =
        sensitivity_scenario_id,
      primary_scenario_id =
        scenario$primary_scenario_id,
      repetition = repetition,
      data_seed = data_seed,
      method_seed = method_seed,
      stringsAsFactors = FALSE
    )
    
    checkpoint <- list(
      task_id = task_id,
      sensitivity_scenario_id =
        sensitivity_scenario_id,
      primary_scenario_id =
        scenario$primary_scenario_id,
      chunk_id = chunk_id,
      rep_start = rep_start,
      rep_end = rep_end,
      nsim_empirical = nsim_empirical,
      reps_per_empirical_chunk =
        reps_per_empirical_chunk,
      next_rep = repetition + 1L,
      results_list = results_list,
      rng_list = rng_list,
      dataset_qc_list = dataset_qc_list
    )
    
    atomic_save_rds(
      checkpoint,
      checkpoint_file
    )
    
    cat(
      "Scenario ",
      sensitivity_scenario_id,
      ", chunk ",
      chunk_id,
      ": completed repetition ",
      repetition,
      " of ",
      rep_end,
      ".\n",
      sep = ""
    )
  }
}

# Validate and publish the completed chunk.
results <- bind_rows(results_list)
rng_states <- bind_rows(rng_list)
dataset_qc <- bind_rows(dataset_qc_list)

expected_repetitions <- seq.int(
  rep_start,
  rep_end
)

method_keys <- results %>%
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

rng_keys <- rng_states %>%
  count(
    sensitivity_scenario_id,
    repetition,
    name = "n_rows"
  )

chunk_checks <- c(
  result_repetitions_complete = setequal(
    results$repetition,
    expected_repetitions
  ),
  dataset_repetitions_complete = setequal(
    dataset_qc$repetition,
    expected_repetitions
  ),
  rng_repetitions_complete = setequal(
    rng_states$repetition,
    expected_repetitions
  ),
  result_row_count = (
    nrow(results) ==
      length(expected_repetitions) *
      length(expected_methods)
  ),
  unique_method_keys = (
    nrow(method_keys) == nrow(results) &&
      all_true(method_keys$n_rows == 1L)
  ),
  unique_dataset_keys = (
    nrow(dataset_keys) ==
      length(expected_repetitions) &&
      all_true(dataset_keys$n_rows == 1L)
  ),
  unique_rng_keys = (
    nrow(rng_keys) ==
      length(expected_repetitions) &&
      all_true(rng_keys$n_rows == 1L)
  ),
  expected_methods = setequal(
    results$method,
    expected_methods
  ),
  empirical_marker_error = all_true(
    !is.na(results$marker_error) &
      results$marker_error == "empirical"
  ),
  production_mi_settings = all_true(
    results$nimp == nimp_primary &
      results$mice_method ==
      mice_method_primary &
      results$mice_maxit ==
      mice_maxit_primary &
      results$jomo_nburn ==
      jomo_nburn_primary &
      results$jomo_nbetween ==
      jomo_nbetween_primary
  ),
  marker_missingness_blockwise = all_true(
    dataset_qc$marker_blockwise
  ),
  phase2_matches_markers = all_true(
    dataset_qc$phase2_matches_markers
  ),
  sufficient_phase2_n = all_true(
    dataset_qc$sufficient_phase2_n
  ),
  valid_probability_range = all_true(
    dataset_qc$finite_probability_range
  )
)

if (!all_true(chunk_checks)) {
  failed_checks <- names(chunk_checks)[
    !chunk_checks
  ]
  
  stop(
    "Completed chunk failed structural QC:\n- ",
    paste(failed_checks, collapse = "\n- ")
  )
}

metadata <- list(
  task_id = task_id,
  sensitivity_scenario_id =
    sensitivity_scenario_id,
  primary_scenario_id =
    scenario$primary_scenario_id,
  scenario_label = scenario$scenario_label,
  scenario_key = scenario$scenario_key,
  chunk_id = chunk_id,
  rep_start = rep_start,
  rep_end = rep_end,
  nsim_empirical = nsim_empirical,
  reps_per_empirical_chunk =
    reps_per_empirical_chunk,
  chunks_per_empirical_scenario =
    chunks_per_empirical_scenario,
  N = scenario$N,
  phase2_fraction = scenario$phase2_fraction,
  r2_a_marker = scenario$r2_a_marker,
  r2_y_marker = scenario$r2_y_marker,
  theta = scenario$theta,
  marker_error = empirical_dgm$marker_error,
  nimp = nimp_primary,
  mice_method = mice_method_primary,
  mice_maxit = mice_maxit_primary,
  jomo_nburn = jomo_nburn_primary,
  jomo_nbetween = jomo_nbetween_primary,
  empirical_seed_master = empirical_seed_master,
  empirical_calibration_md5 = unname(
    tools::md5sum(empirical_calibration_file)
  ),
  empirical_dgm_md5 = unname(
    tools::md5sum(empirical_dgm_file)
  ),
  created = as.character(Sys.time())
)

atomic_save_rds(
  list(
    metadata = metadata,
    results = results,
    rng_states = rng_states,
    dataset_qc = dataset_qc
  ),
  final_file
)

if (file.exists(checkpoint_file)) {
  unlink(checkpoint_file)
}

cat(
  "\nEmpirical production chunk complete:\n",
  final_file,
  "\n",
  sep = ""
)
