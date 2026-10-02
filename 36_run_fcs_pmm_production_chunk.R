################################################################################
# 36_run_fcs_pmm_production_chunk.R
#
# Slurm-array worker for the FCS-PMM sensitivity.
#
# Defaults: four scenarios, 2,000 repetitions per scenario, 10 repetitions per
# task, and 800 array tasks. Only FCS-PMM is fitted in production.
################################################################################

source("00_config.R")
source("02_dgm.R")
source("03_methods.R")
source("34_fcs_pmm_sensitivity_helpers.R")

library(dplyr)


# ------------------------------------------------------------------------------
# 1. Settings and array-task mapping
# ------------------------------------------------------------------------------

nsim_fcs_pmm <- read_fcs_pmm_positive_integer_env(
  "SIM_FCS_PMM_NSIM",
  2000L
)

reps_per_fcs_pmm_chunk <- read_fcs_pmm_positive_integer_env(
  "SIM_FCS_PMM_REPS_PER_CHUNK",
  10L
)

if (nsim_fcs_pmm != 2000L || reps_per_fcs_pmm_chunk != 10L) {
  stop(
    "The prespecified FCS-PMM production run requires ",
    "SIM_FCS_PMM_NSIM=2000 and SIM_FCS_PMM_REPS_PER_CHUNK=10."
  )
}

chunks_per_fcs_pmm_scenario <- as.integer(
  ceiling(nsim_fcs_pmm / reps_per_fcs_pmm_chunk)
)

n_fcs_pmm_array_tasks <-
  nrow(fcs_pmm_sensitivity_grid) * chunks_per_fcs_pmm_scenario

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
    task_id > n_fcs_pmm_array_tasks
) {
  stop(
    "Array task id must be between 1 and ",
    n_fcs_pmm_array_tasks,
    "; received '",
    task_id_text,
    "'."
  )
}

sensitivity_scenario_id <-
  ((task_id - 1L) %/% chunks_per_fcs_pmm_scenario) + 1L

chunk_id <- ((task_id - 1L) %% chunks_per_fcs_pmm_scenario) + 1L
rep_start <- ((chunk_id - 1L) * reps_per_fcs_pmm_chunk) + 1L
rep_end <- min(chunk_id * reps_per_fcs_pmm_chunk, nsim_fcs_pmm)

scenario_row <- fcs_pmm_sensitivity_grid[
  fcs_pmm_sensitivity_grid$sensitivity_scenario_id ==
    sensitivity_scenario_id,
  ,
  drop = FALSE
]

if (nrow(scenario_row) != 1L) {
  stop(
    "Could not uniquely identify FCS-PMM sensitivity scenario ",
    sensitivity_scenario_id,
    "."
  )
}


# ------------------------------------------------------------------------------
# 2. Gates, primary DGM, paths, and safe restart behavior
# ------------------------------------------------------------------------------

fcs_pmm_results_dir <- file.path(results_dir, "fcs_pmm_sensitivity")
fcs_pmm_canary_dir <- file.path(fcs_pmm_results_dir, "canary")
fcs_pmm_production_dir <- file.path(fcs_pmm_results_dir, "production")
fcs_pmm_chunk_dir <- file.path(fcs_pmm_production_dir, "chunks")

scenario_chunk_dir <- file.path(
  fcs_pmm_chunk_dir,
  sprintf(
    "scenario_%02d_%s",
    sensitivity_scenario_id,
    scenario_row$scenario_key
  )
)

dir.create(scenario_chunk_dir, showWarnings = FALSE, recursive = TRUE)

mcar_summary_pass_file <- file.path(
  results_dir,
  "phase2_mcar_sensitivity",
  "production",
  "summary",
  "PHASE2_MCAR_PRODUCTION_SUMMARY_PASS.txt"
)

canary_pass_file <- file.path(
  fcs_pmm_canary_dir,
  "FCS_PMM_CANARY_PASS.txt"
)

canary_qc_file <- file.path(
  fcs_pmm_canary_dir,
  "fcs_pmm_canary_qc.csv"
)

canary_defaults_file <- file.path(
  fcs_pmm_canary_dir,
  "fcs_pmm_default_settings.csv"
)

dgm_file <- fcs_pmm_primary_cache_filename(scenario_row)

required_files <- c(
  mcar_summary_pass_file,
  canary_pass_file,
  canary_qc_file,
  canary_defaults_file,
  calibration_rds,
  dgm_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required FCS-PMM production files are missing:\n- ",
    paste(missing_files, collapse = "\n- "),
    "\nRun and review 35_run_fcs_pmm_canary.R first."
  )
}

if (!any(grepl(
  "MCAR PRODUCTION SUMMARY: PASS",
  readLines(mcar_summary_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The Phase-2 MCAR production summary PASS marker is invalid.")
}

if (!any(grepl(
  "FCS-PMM CANARY: PASS",
  readLines(canary_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The FCS-PMM canary PASS marker is invalid.")
}

canary_qc <- readr::read_csv(canary_qc_file, show_col_types = FALSE)
canary_defaults <- readr::read_csv(
  canary_defaults_file,
  show_col_types = FALSE
)

if (
  nrow(canary_qc) == 0L ||
    !fcs_pmm_all_true(canary_qc$passed) ||
    !fcs_pmm_all_true(canary_qc$overall_canary_pass)
) {
  stop("The FCS-PMM canary QC file does not record a full pass.")
}

calibration <- readRDS(calibration_rds)
dgm <- readRDS(dgm_file)
dgm_qc <- validate_fcs_pmm_primary_dgm(dgm, scenario_row)
pmm_defaults <- fcs_pmm_default_settings()

if (!fcs_pmm_all_true(dgm_qc$dgm_matches_exactly)) {
  stop("The matched primary DGM cache failed FCS-PMM structural validation.")
}

if (!fcs_pmm_all_true(pmm_defaults$defaults_match_prespecified)) {
  stop("The installed mice PMM defaults do not match the validated canary.")
}

if (!isTRUE(all.equal(
  pmm_defaults,
  canary_defaults,
  tolerance = 0,
  check.attributes = FALSE
))) {
  stop("The current mice version or PMM defaults differ from the canary.")
}

final_file <- file.path(
  scenario_chunk_dir,
  sprintf(
    "fcs_pmm_scenario_%02d_chunk_%03d.rds",
    sensitivity_scenario_id,
    chunk_id
  )
)

checkpoint_file <- file.path(
  scenario_chunk_dir,
  sprintf(
    "fcs_pmm_scenario_%02d_chunk_%03d_checkpoint.rds",
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

  expected_repetitions <- seq.int(rep_start, rep_end)

  existing_final_matches <-
    has_required_components &&
    identical(
      existing_final$metadata$sensitivity_scenario_id,
      sensitivity_scenario_id
    ) &&
    identical(existing_final$metadata$chunk_id, chunk_id) &&
    identical(existing_final$metadata$rep_start, rep_start) &&
    identical(existing_final$metadata$rep_end, rep_end) &&
    identical(existing_final$metadata$nsim_fcs_pmm, nsim_fcs_pmm) &&
    identical(
      existing_final$metadata$reps_per_fcs_pmm_chunk,
      reps_per_fcs_pmm_chunk
    ) &&
    identical(existing_final$metadata$mice_method, "pmm") &&
    nrow(existing_final$results) == length(expected_repetitions) &&
    nrow(existing_final$rng_states) == length(expected_repetitions) &&
    nrow(existing_final$dataset_qc) == length(expected_repetitions) &&
    setequal(
      unique(existing_final$results$repetition),
      expected_repetitions
    ) &&
    setequal(
      unique(existing_final$rng_states$repetition),
      expected_repetitions
    ) &&
    setequal(
      unique(existing_final$dataset_qc$repetition),
      expected_repetitions
    ) &&
    identical(unique(existing_final$results$method), fcs_pmm_label)

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
    "nsim_fcs_pmm",
    "reps_per_fcs_pmm_chunk",
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
    identical(checkpoint$nsim_fcs_pmm, nsim_fcs_pmm) &&
    identical(
      checkpoint$reps_per_fcs_pmm_chunk,
      reps_per_fcs_pmm_chunk
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
    "Resuming FCS-PMM scenario ",
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
  "Task ID: ", task_id, " / ", n_fcs_pmm_array_tasks, "\n",
  "Sensitivity scenario: ", sensitivity_scenario_id, " / 4\n",
  "Primary scenario ID: ", scenario_row$primary_scenario_id, "\n",
  "Scenario: ", scenario_row$scenario_label, "\n",
  "Chunk: ", chunk_id, " / ", chunks_per_fcs_pmm_scenario, "\n",
  "Repetitions: ", rep_start, "-", rep_end, "\n",
  "N: ", scenario_row$N, "\n",
  "Phase-2 fraction: ", scenario_row$phase2_fraction, "\n",
  "R2 A~M|X: ", scenario_row$r2_a_marker, "\n",
  "R2 Y~M|A,X: ", scenario_row$r2_y_marker, "\n",
  "Theta: ", scenario_row$theta, "\n",
  "Selection: MAR\n",
  "Marker error: MVN\n",
  "FCS method: PMM\n",
  "MICE imputations: ", nimp_primary, "\n",
  "MICE maxit: ", mice_maxit_primary, "\n",
  "PMM donors: ", pmm_defaults$donors, "\n",
  "============================================================\n",
  sep = ""
)

if (next_rep <= rep_end) {

  for (repetition in seq.int(next_rep, rep_end)) {

    states <- fcs_pmm_get_repetition_states(
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

    structural_qc <- fcs_pmm_dataset_qc(
      dat = dat,
      marker_names = calibration$marker_names,
      row = scenario_row
    )

    if (
      !structural_qc$marker_blockwise ||
        !structural_qc$phase2_matches_markers ||
        !structural_qc$finite_valid_probabilities ||
        !structural_qc$sufficient_phase2_n
    ) {
      stop(
        "Structural DGM QC failed for FCS-PMM scenario ",
        sensitivity_scenario_id,
        ", repetition ",
        repetition,
        "."
      )
    }

    list_index <- repetition - rep_start + 1L

    dataset_qc_list[[list_index]] <- bind_cols(
      data.frame(
        sensitivity_scenario_id = sensitivity_scenario_id,
        primary_scenario_id = scenario_row$primary_scenario_id,
        scenario_label = scenario_row$scenario_label,
        scenario_key = scenario_row$scenario_key,
        task_id = task_id,
        chunk_id = chunk_id,
        repetition = repetition,
        stringsAsFactors = FALSE
      ),
      structural_qc
    )

    .Random.seed <- states$method_state

    ans <- fit_fcs_pmm_sensitivity(
      dat = dat,
      marker_names = calibration$marker_names
    ) %>%
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
        mice_method = "pmm",
        nimp = nimp_primary,
        mice_maxit = mice_maxit_primary,
        pmm_donors = pmm_defaults$donors,
        pmm_matchtype = pmm_defaults$matchtype,
        pmm_ridge = pmm_defaults$ridge,
        mice_version = pmm_defaults$mice_version
      )

    results_list[[list_index]] <- ans

    rng_list[[list_index]] <- data.frame(
      sensitivity_scenario_id = sensitivity_scenario_id,
      primary_scenario_id = scenario_row$primary_scenario_id,
      repetition = repetition,
      data_rng_state = fcs_pmm_state_to_string(states$data_state),
      method_rng_state = fcs_pmm_state_to_string(states$method_state),
      stringsAsFactors = FALSE
    )

    checkpoint <- list(
      task_id = task_id,
      sensitivity_scenario_id = sensitivity_scenario_id,
      primary_scenario_id = scenario_row$primary_scenario_id,
      chunk_id = chunk_id,
      rep_start = rep_start,
      rep_end = rep_end,
      nsim_fcs_pmm = nsim_fcs_pmm,
      reps_per_fcs_pmm_chunk = reps_per_fcs_pmm_chunk,
      next_rep = repetition + 1L,
      results_list = results_list,
      rng_list = rng_list,
      dataset_qc_list = dataset_qc_list
    )

    fcs_pmm_atomic_save_rds(
      checkpoint,
      checkpoint_file,
      compress = FALSE
    )

    cat(
      "FCS-PMM scenario ",
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

result_keys <- results %>%
  count(sensitivity_scenario_id, repetition, method, name = "n_rows")

chunk_checks <- c(
  setequal(unique(results$repetition), expected_repetitions),
  setequal(unique(dataset_qc$repetition), expected_repetitions),
  setequal(unique(rng_states$repetition), expected_repetitions),
  nrow(results) == length(expected_repetitions),
  nrow(dataset_qc) == length(expected_repetitions),
  nrow(rng_states) == length(expected_repetitions),
  nrow(result_keys) == nrow(results) &&
    fcs_pmm_all_true(result_keys$n_rows == 1L),
  identical(unique(results$method), fcs_pmm_label),
  fcs_pmm_all_true(results$mice_method == "pmm"),
  fcs_pmm_all_true(results$nimp == nimp_primary),
  fcs_pmm_all_true(results$mice_maxit == mice_maxit_primary),
  fcs_pmm_all_true(results$mice_version == pmm_defaults$mice_version),
  fcs_pmm_all_true(results$selection == "mar"),
  fcs_pmm_all_true(results$marker_error == "mvn"),
  fcs_pmm_all_true(dataset_qc$marker_blockwise),
  fcs_pmm_all_true(dataset_qc$phase2_matches_markers),
  fcs_pmm_all_true(dataset_qc$finite_valid_probabilities),
  fcs_pmm_all_true(dataset_qc$sufficient_phase2_n)
)

if (!fcs_pmm_all_true(chunk_checks)) {
  stop("Completed FCS-PMM chunk failed final structural QC.")
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
  nsim_fcs_pmm = nsim_fcs_pmm,
  reps_per_fcs_pmm_chunk = reps_per_fcs_pmm_chunk,
  chunks_per_fcs_pmm_scenario = chunks_per_fcs_pmm_scenario,
  N = scenario_row$N,
  phase2_fraction = scenario_row$phase2_fraction,
  r2_a_marker = scenario_row$r2_a_marker,
  r2_y_marker = scenario_row$r2_y_marker,
  theta = scenario_row$theta,
  selection = dgm$selection,
  marker_error = dgm$marker_error,
  mice_method = "pmm",
  nimp = nimp_primary,
  mice_maxit = mice_maxit_primary,
  pmm_donors = pmm_defaults$donors,
  pmm_matchtype = pmm_defaults$matchtype,
  pmm_ridge = pmm_defaults$ridge,
  mice_version = pmm_defaults$mice_version,
  primary_rng_seed_master = rng_seed_master,
  calibration_md5 = unname(tools::md5sum(calibration_rds)),
  primary_dgm_md5 = unname(tools::md5sum(dgm_file)),
  created = as.character(Sys.time())
)

fcs_pmm_atomic_save_rds(
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
  "\nFCS-PMM production chunk complete:\n",
  final_file,
  "\n",
  sep = ""
)
