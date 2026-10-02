################################################################################
# 31_run_phase2_mcar_production_chunk.R
#
# Slurm-array worker for the Phase-2 MCAR sensitivity.
#
# Defaults: four scenarios, 2,000 repetitions per scenario, 10 repetitions per
# task, and 800 array tasks. Full-data and method RNG streams are paired with the
# corresponding primary scenario/repetition.
################################################################################

source("00_config.R")
source("02_dgm.R")
source("03_methods.R")
source("28_phase2_mcar_sensitivity_helpers.R")

library(dplyr)
library(readr)


# ------------------------------------------------------------------------------
# 1. Settings and array-task mapping
# ------------------------------------------------------------------------------

nsim_mcar <- read_mcar_positive_integer_env("SIM_MCAR_NSIM", 2000L)

expected_methods <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

reps_per_mcar_chunk <- read_mcar_positive_integer_env(
  "SIM_MCAR_REPS_PER_CHUNK",
  10L
)

if (nsim_mcar != 2000L || reps_per_mcar_chunk != 10L) {
  stop(
    "The prespecified MCAR production run requires SIM_MCAR_NSIM=2000 ",
    "and SIM_MCAR_REPS_PER_CHUNK=10."
  )
}

chunks_per_mcar_scenario <- as.integer(
  ceiling(nsim_mcar / reps_per_mcar_chunk)
)

n_mcar_array_tasks <-
  nrow(mcar_sensitivity_grid) * chunks_per_mcar_scenario

task_id_text <- Sys.getenv(
  "SLURM_ARRAY_TASK_ID",
  unset = Sys.getenv("SIM_ARRAY_TASK_ID", unset = "")
)

if (!nzchar(task_id_text)) {
  stop(
    "No array task id found. On Torch use SLURM_ARRAY_TASK_ID. ",
    "For a local mapping test, set SIM_ARRAY_TASK_ID."
  )
}

task_id <- suppressWarnings(as.integer(task_id_text))

if (
  is.na(task_id) ||
    task_id < 1L ||
    task_id > n_mcar_array_tasks
) {
  stop(
    "Array task id must be between 1 and ",
    n_mcar_array_tasks,
    "; received '",
    task_id_text,
    "'."
  )
}

sensitivity_scenario_id <-
  ((task_id - 1L) %/% chunks_per_mcar_scenario) + 1L

chunk_id <- ((task_id - 1L) %% chunks_per_mcar_scenario) + 1L
rep_start <- ((chunk_id - 1L) * reps_per_mcar_chunk) + 1L
rep_end <- min(chunk_id * reps_per_mcar_chunk, nsim_mcar)

scenario_row <- mcar_sensitivity_grid[
  mcar_sensitivity_grid$sensitivity_scenario_id ==
    sensitivity_scenario_id,
  ,
  drop = FALSE
]

if (nrow(scenario_row) != 1L) {
  stop(
    "Could not uniquely identify MCAR sensitivity scenario ",
    sensitivity_scenario_id,
    "."
  )
}


# ------------------------------------------------------------------------------
# 2. Paths, validated caches, and safe restart behavior
# ------------------------------------------------------------------------------

mcar_results_dir <- file.path(results_dir, "phase2_mcar_sensitivity")
mcar_production_dir <- file.path(mcar_results_dir, "production")
mcar_cache_dir <- file.path(mcar_production_dir, "cache")
mcar_dgm_dir <- file.path(mcar_cache_dir, "dgm")
mcar_chunk_dir <- file.path(mcar_production_dir, "chunks")

scenario_chunk_dir <- file.path(
  mcar_chunk_dir,
  sprintf(
    "scenario_%02d_%s",
    sensitivity_scenario_id,
    scenario_row$scenario_key
  )
)

dir.create(scenario_chunk_dir, showWarnings = FALSE, recursive = TRUE)

cache_pass_file <- file.path(
  mcar_cache_dir,
  "PHASE2_MCAR_PRODUCTION_CACHE_PASS.txt"
)

cache_qc_file <- file.path(
  mcar_cache_dir,
  "phase2_mcar_production_cache_qc.csv"
)

mcar_dgm_file <- mcar_dgm_cache_filename(scenario_row, mcar_dgm_dir)

required_files <- c(
  cache_pass_file,
  cache_qc_file,
  calibration_rds,
  mcar_dgm_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required Phase-2 MCAR production files are missing:\n- ",
    paste(missing_files, collapse = "\n- "),
    "\nRun 30_build_phase2_mcar_production_cache.R first."
  )
}

if (!any(grepl(
  "MCAR PRODUCTION CACHE: PASS",
  readLines(cache_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The Phase-2 MCAR production cache PASS marker is invalid.")
}

cache_qc <- read_csv(cache_qc_file, show_col_types = FALSE)

if (
  nrow(cache_qc) == 0L ||
    !mcar_all_true(cache_qc$passed) ||
    !mcar_all_true(cache_qc$overall_cache_pass)
) {
  stop("The MCAR production cache QC file does not record a full pass.")
}

calibration <- readRDS(calibration_rds)
dgm <- readRDS(mcar_dgm_file)

if (!identical(dgm$selection, "mcar")) {
  stop("Cached DGM does not use MCAR selection.")
}

if (!identical(dgm$marker_error, "mvn")) {
  stop("Cached MCAR DGM does not retain MVN marker errors.")
}

final_file <- file.path(
  scenario_chunk_dir,
  sprintf(
    "mcar_scenario_%02d_chunk_%03d.rds",
    sensitivity_scenario_id,
    chunk_id
  )
)

checkpoint_file <- file.path(
  scenario_chunk_dir,
  sprintf(
    "mcar_scenario_%02d_chunk_%03d_checkpoint.rds",
    sensitivity_scenario_id,
    chunk_id
  )
)

if (file.exists(final_file)) {

  existing_final <- tryCatch(readRDS(final_file), error = function(e) NULL)

  required_components <- c(
    "metadata",
    "results",
    "rng_states",
    "dataset_qc"
  )

  has_required_components <-
    !is.null(existing_final) &&
    all(required_components %in% names(existing_final))

  existing_final_matches <-
    has_required_components &&
    identical(
      existing_final$metadata$sensitivity_scenario_id,
      sensitivity_scenario_id
    ) &&
    identical(existing_final$metadata$chunk_id, chunk_id) &&
    identical(existing_final$metadata$rep_start, rep_start) &&
    identical(existing_final$metadata$rep_end, rep_end) &&
    identical(existing_final$metadata$nsim_mcar, nsim_mcar) &&
    identical(
      existing_final$metadata$reps_per_mcar_chunk,
      reps_per_mcar_chunk
    ) &&
    identical(existing_final$metadata$selection, "mcar") &&
    identical(existing_final$metadata$marker_error, "mvn") &&
    nrow(existing_final$results) ==
      length(seq.int(rep_start, rep_end)) * length(expected_methods) &&
    nrow(existing_final$rng_states) == length(seq.int(rep_start, rep_end)) &&
    nrow(existing_final$dataset_qc) == length(seq.int(rep_start, rep_end)) &&
    setequal(
      unique(existing_final$results$repetition),
      seq.int(rep_start, rep_end)
    ) &&
    setequal(unique(existing_final$results$method), expected_methods)

  if (!existing_final_matches) {
    stop(
      "Existing final chunk is unreadable or does not match this task:\n",
      final_file
    )
  }

  cat("Final chunk already exists; exiting:\n", final_file, "\n", sep = "")
  quit(save = "no", status = 0L)
}


# ------------------------------------------------------------------------------
# 3. Resume from a matching checkpoint
# ------------------------------------------------------------------------------

results_list <- list()
rng_list <- list()
dataset_qc_list <- list()
next_rep <- rep_start

if (file.exists(checkpoint_file)) {

  checkpoint <- tryCatch(
    readRDS(checkpoint_file),
    error = function(e) NULL
  )

  required_checkpoint_components <- c(
    "sensitivity_scenario_id",
    "chunk_id",
    "rep_start",
    "rep_end",
    "nsim_mcar",
    "reps_per_mcar_chunk",
    "next_rep",
    "results_list",
    "rng_list",
    "dataset_qc_list"
  )

  checkpoint_complete <-
    !is.null(checkpoint) &&
    all(required_checkpoint_components %in% names(checkpoint))

  checkpoint_matches <-
    checkpoint_complete &&
    identical(checkpoint$sensitivity_scenario_id, sensitivity_scenario_id) &&
    identical(checkpoint$chunk_id, chunk_id) &&
    identical(checkpoint$rep_start, rep_start) &&
    identical(checkpoint$rep_end, rep_end) &&
    identical(checkpoint$nsim_mcar, nsim_mcar) &&
    identical(
      checkpoint$reps_per_mcar_chunk,
      reps_per_mcar_chunk
    ) &&
    checkpoint$next_rep >= rep_start &&
    checkpoint$next_rep <= rep_end + 1L &&
    is.list(checkpoint$results_list) &&
    is.list(checkpoint$rng_list) &&
    is.list(checkpoint$dataset_qc_list)

  if (!checkpoint_matches) {
    stop(
      "Existing checkpoint does not match the requested task:\n",
      checkpoint_file
    )
  }

  results_list <- checkpoint$results_list
  rng_list <- checkpoint$rng_list
  dataset_qc_list <- checkpoint$dataset_qc_list
  next_rep <- checkpoint$next_rep

  cat(
    "Resuming MCAR scenario ",
    sensitivity_scenario_id,
    ", chunk ",
    chunk_id,
    " at repetition ",
    next_rep,
    ".\n",
    sep = ""
  )
}


# ------------------------------------------------------------------------------
# 4. Run assigned repetitions
# ------------------------------------------------------------------------------

cat(
  "\n============================================================\n",
  "Task ID: ", task_id, " / ", n_mcar_array_tasks, "\n",
  "Sensitivity scenario: ", sensitivity_scenario_id, " / 4\n",
  "Primary scenario ID: ", scenario_row$primary_scenario_id, "\n",
  "Scenario: ", scenario_row$scenario_label, "\n",
  "Chunk: ", chunk_id, " / ", chunks_per_mcar_scenario, "\n",
  "Repetitions: ", rep_start, "-", rep_end, "\n",
  "N: ", scenario_row$N, "\n",
  "Phase-2 fraction: ", scenario_row$phase2_fraction, "\n",
  "R2 A~M|X: ", scenario_row$r2_a_marker, "\n",
  "R2 Y~M|A,X: ", scenario_row$r2_y_marker, "\n",
  "Theta: ", scenario_row$theta, "\n",
  "Selection: MCAR\n",
  "Marker error: MVN\n",
  "MICE imputations: ", nimp_primary, "\n",
  "MICE maxit: ", mice_maxit_primary, "\n",
  "JM-MI nburn: ", jomo_nburn_primary, "\n",
  "JM-MI nbetween: ", jomo_nbetween_primary, "\n",
  "============================================================\n",
  sep = ""
)

if (next_rep <= rep_end) {

  for (repetition in seq.int(next_rep, rep_end)) {

    states <- mcar_get_repetition_states(
      primary_scenario_id = scenario_row$primary_scenario_id,
      repetition = repetition
    )

    .Random.seed <- states$data_state

    simulated <- generate_two_phase_data(
      N = scenario_row$N,
      calibration = calibration,
      dgm_parameters = dgm
    )

    dat <- simulated$observed

    n_marker_observed <- rowSums(
      !is.na(dat[, calibration$marker_names, drop = FALSE])
    )

    marker_blockwise <- all(
      n_marker_observed %in% c(0L, length(calibration$marker_names))
    )

    phase2_matches_markers <- identical(
      as.integer(n_marker_observed == length(calibration$marker_names)),
      as.integer(dat$phase2)
    )

    constant_true_probability <- max(
      abs(dat$pi_true - scenario_row$phase2_fraction)
    ) <= 1e-15

    sufficient_phase2_n <-
      sum(dat$phase2) > length(calibration$marker_names) + 10L

    if (
      !marker_blockwise ||
        !phase2_matches_markers ||
        !constant_true_probability ||
        !sufficient_phase2_n
    ) {
      stop(
        "Structural DGM QC failed for MCAR scenario ",
        sensitivity_scenario_id,
        ", repetition ",
        repetition,
        "."
      )
    }

    list_index <- repetition - rep_start + 1L

    dataset_qc_list[[list_index]] <- data.frame(
      sensitivity_scenario_id = sensitivity_scenario_id,
      primary_scenario_id = scenario_row$primary_scenario_id,
      scenario_label = scenario_row$scenario_label,
      scenario_key = scenario_row$scenario_key,
      task_id = task_id,
      chunk_id = chunk_id,
      repetition = repetition,
      marker_blockwise = marker_blockwise,
      phase2_matches_markers = phase2_matches_markers,
      n_phase2 = sum(dat$phase2),
      realized_phase2_fraction = mean(dat$phase2),
      min_pi_true = min(dat$pi_true),
      max_pi_true = max(dat$pi_true),
      max_abs_pi_deviation = max(
        abs(dat$pi_true - scenario_row$phase2_fraction)
      ),
      constant_true_probability = constant_true_probability,
      sufficient_phase2_n = sufficient_phase2_n,
      stringsAsFactors = FALSE
    )

    .Random.seed <- states$method_state

    ans <- run_all_mcar_methods_timed(
      dat = dat,
      marker_names = calibration$marker_names
    )

    ans <- ans %>%
      mutate(
        task_id = task_id,
        sensitivity_scenario_id = sensitivity_scenario_id,
        primary_scenario_id = scenario_row$primary_scenario_id,
        scenario_label = scenario_row$scenario_label,
        scenario_key = scenario_row$scenario_key,
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
        max_pi_true = max(dat$pi_true),
        selection = dgm$selection,
        marker_error = dgm$marker_error,
        nimp = nimp_primary,
        mice_method = mice_method_primary,
        mice_maxit = mice_maxit_primary,
        jomo_nburn = jomo_nburn_primary,
        jomo_nbetween = jomo_nbetween_primary
      )

    results_list[[list_index]] <- ans

    rng_list[[list_index]] <- data.frame(
      sensitivity_scenario_id = sensitivity_scenario_id,
      primary_scenario_id = scenario_row$primary_scenario_id,
      repetition = repetition,
      data_rng_state = mcar_state_to_string(states$data_state),
      method_rng_state = mcar_state_to_string(states$method_state),
      stringsAsFactors = FALSE
    )

    checkpoint <- list(
      task_id = task_id,
      sensitivity_scenario_id = sensitivity_scenario_id,
      primary_scenario_id = scenario_row$primary_scenario_id,
      chunk_id = chunk_id,
      rep_start = rep_start,
      rep_end = rep_end,
      nsim_mcar = nsim_mcar,
      reps_per_mcar_chunk = reps_per_mcar_chunk,
      next_rep = repetition + 1L,
      results_list = results_list,
      rng_list = rng_list,
      dataset_qc_list = dataset_qc_list
    )

    mcar_atomic_save_rds(checkpoint, checkpoint_file, compress = FALSE)

    cat(
      "MCAR scenario ",
      sensitivity_scenario_id,
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


# ------------------------------------------------------------------------------
# 5. Validate and publish completed chunk
# ------------------------------------------------------------------------------

results <- bind_rows(results_list)
rng_states <- bind_rows(rng_list)
dataset_qc <- bind_rows(dataset_qc_list)
expected_repetitions <- seq.int(rep_start, rep_end)

method_keys <- results %>%
  count(
    sensitivity_scenario_id,
    repetition,
    method,
    name = "n_rows"
  )

chunk_checks <- c(
  setequal(unique(results$repetition), expected_repetitions),
  setequal(unique(dataset_qc$repetition), expected_repetitions),
  setequal(unique(rng_states$repetition), expected_repetitions),
  nrow(results) == length(expected_repetitions) * length(expected_methods),
  nrow(dataset_qc) == length(expected_repetitions),
  nrow(rng_states) == length(expected_repetitions),
  nrow(method_keys) == nrow(results) && all(method_keys$n_rows == 1L),
  setequal(unique(results$method), expected_methods),
  mcar_all_true(results$selection == "mcar"),
  mcar_all_true(results$marker_error == "mvn"),
  mcar_all_true(results$nimp == nimp_primary),
  mcar_all_true(results$mice_method == mice_method_primary),
  mcar_all_true(results$mice_maxit == mice_maxit_primary),
  mcar_all_true(results$jomo_nburn == jomo_nburn_primary),
  mcar_all_true(results$jomo_nbetween == jomo_nbetween_primary),
  mcar_all_true(dataset_qc$marker_blockwise),
  mcar_all_true(dataset_qc$phase2_matches_markers),
  mcar_all_true(dataset_qc$constant_true_probability),
  mcar_all_true(dataset_qc$sufficient_phase2_n)
)

if (!mcar_all_true(chunk_checks)) {
  stop("Completed Phase-2 MCAR chunk failed final structural QC.")
}

metadata <- list(
  task_id = task_id,
  sensitivity_scenario_id = sensitivity_scenario_id,
  primary_scenario_id = scenario_row$primary_scenario_id,
  scenario_label = scenario_row$scenario_label,
  scenario_key = scenario_row$scenario_key,
  chunk_id = chunk_id,
  rep_start = rep_start,
  rep_end = rep_end,
  nsim_mcar = nsim_mcar,
  reps_per_mcar_chunk = reps_per_mcar_chunk,
  chunks_per_mcar_scenario = chunks_per_mcar_scenario,
  N = scenario_row$N,
  phase2_fraction = scenario_row$phase2_fraction,
  r2_a_marker = scenario_row$r2_a_marker,
  r2_y_marker = scenario_row$r2_y_marker,
  theta = scenario_row$theta,
  selection = dgm$selection,
  marker_error = dgm$marker_error,
  nimp = nimp_primary,
  mice_method = mice_method_primary,
  mice_maxit = mice_maxit_primary,
  jomo_nburn = jomo_nburn_primary,
  jomo_nbetween = jomo_nbetween_primary,
  primary_rng_seed_master = rng_seed_master,
  calibration_md5 = unname(tools::md5sum(calibration_rds)),
  mcar_dgm_md5 = unname(tools::md5sum(mcar_dgm_file)),
  created = as.character(Sys.time())
)

mcar_atomic_save_rds(
  list(
    metadata = metadata,
    results = results,
    rng_states = rng_states,
    dataset_qc = dataset_qc
  ),
  final_file,
  compress = FALSE
)

if (file.exists(checkpoint_file)) {
  unlink(checkpoint_file)
}

cat(
  "\nPhase-2 MCAR production chunk complete:\n",
  final_file,
  "\n",
  sep = ""
)
