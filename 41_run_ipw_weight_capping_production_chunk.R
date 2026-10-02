################################################################################
# 41_run_ipw_weight_capping_production_chunk.R
#
# Slurm-array worker for the IPW 1st/99th-percentile capping sensitivity.
#
# Defaults: four scenarios, 2,000 repetitions per scenario, 100 repetitions per
# task, and 80 array tasks. Only capped IPW is fitted in production.
################################################################################

source("00_config.R")
source("02_dgm.R")
source("03_methods.R")
source("39_ipw_weight_capping_helpers.R")

library(dplyr)


# ------------------------------------------------------------------------------
# 1. Settings and array-task mapping
# ------------------------------------------------------------------------------

nsim_ipw_cap <- read_ipw_cap_positive_integer_env(
  "SIM_IPW_CAP_NSIM",
  2000L
)

reps_per_ipw_cap_chunk <- read_ipw_cap_positive_integer_env(
  "SIM_IPW_CAP_REPS_PER_CHUNK",
  100L
)

if (nsim_ipw_cap != 2000L || reps_per_ipw_cap_chunk != 100L) {
  stop(
    "The prespecified IPW-capping production run requires ",
    "SIM_IPW_CAP_NSIM=2000 and SIM_IPW_CAP_REPS_PER_CHUNK=100."
  )
}

chunks_per_ipw_cap_scenario <- as.integer(
  ceiling(nsim_ipw_cap / reps_per_ipw_cap_chunk)
)

n_ipw_cap_array_tasks <-
  nrow(ipw_cap_sensitivity_grid) * chunks_per_ipw_cap_scenario

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

task_id_numeric <- suppressWarnings(as.numeric(task_id_text))

if (
  !is.finite(task_id_numeric) ||
    task_id_numeric < 1 ||
    task_id_numeric != floor(task_id_numeric) ||
    task_id_numeric > n_ipw_cap_array_tasks
) {
  stop(
    "Array task id must be between 1 and ",
    n_ipw_cap_array_tasks,
    "; received '",
    task_id_text,
    "'."
  )
}

task_id <- as.integer(task_id_numeric)

sensitivity_scenario_id <-
  ((task_id - 1L) %/% chunks_per_ipw_cap_scenario) + 1L

chunk_id <- ((task_id - 1L) %% chunks_per_ipw_cap_scenario) + 1L
rep_start <- ((chunk_id - 1L) * reps_per_ipw_cap_chunk) + 1L
rep_end <- min(chunk_id * reps_per_ipw_cap_chunk, nsim_ipw_cap)

scenario_row <- ipw_cap_sensitivity_grid[
  ipw_cap_sensitivity_grid$sensitivity_scenario_id ==
    sensitivity_scenario_id,
  ,
  drop = FALSE
]

if (nrow(scenario_row) != 1L) {
  stop(
    "Could not uniquely identify IPW-capping sensitivity scenario ",
    sensitivity_scenario_id,
    "."
  )
}


# ------------------------------------------------------------------------------
# 2. Gates, primary DGM, paths, and safe restart behavior
# ------------------------------------------------------------------------------

ipw_cap_results_dir <- file.path(results_dir, "ipw_weight_capping_sensitivity")
ipw_cap_canary_dir <- file.path(ipw_cap_results_dir, "canary")
ipw_cap_production_dir <- file.path(ipw_cap_results_dir, "production")
ipw_cap_chunk_dir <- file.path(ipw_cap_production_dir, "chunks")

scenario_chunk_dir <- file.path(
  ipw_cap_chunk_dir,
  sprintf(
    "scenario_%02d_%s",
    sensitivity_scenario_id,
    scenario_row$scenario_key
  )
)

dir.create(scenario_chunk_dir, showWarnings = FALSE, recursive = TRUE)

fcs_pmm_summary_pass_file <- file.path(
  results_dir,
  "fcs_pmm_sensitivity",
  "production",
  "summary",
  "FCS_PMM_PRODUCTION_SUMMARY_PASS.txt"
)

fcs_pmm_summary_qc_file <- file.path(
  results_dir,
  "fcs_pmm_sensitivity",
  "production",
  "summary",
  "fcs_pmm_summary_qc.csv"
)

canary_pass_file <- file.path(
  ipw_cap_canary_dir,
  "IPW_WEIGHT_CAPPING_CANARY_PASS.txt"
)

canary_qc_file <- file.path(
  ipw_cap_canary_dir,
  "ipw_weight_capping_canary_qc.csv"
)

canary_specification_file <- file.path(
  ipw_cap_canary_dir,
  "ipw_weight_capping_specification.csv"
)

dgm_file <- ipw_cap_primary_cache_filename(scenario_row)

required_files <- c(
  fcs_pmm_summary_pass_file,
  fcs_pmm_summary_qc_file,
  canary_pass_file,
  canary_qc_file,
  canary_specification_file,
  calibration_rds,
  dgm_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required IPW-capping production files are missing:\n- ",
    paste(missing_files, collapse = "\n- "),
    "\nRun and review 40_run_ipw_weight_capping_canary.R first."
  )
}

if (!any(grepl(
  "FCS-PMM PRODUCTION SUMMARY: PASS",
  readLines(fcs_pmm_summary_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The FCS-PMM production summary PASS marker is invalid.")
}

fcs_pmm_summary_qc <- readr::read_csv(
  fcs_pmm_summary_qc_file,
  show_col_types = FALSE
)

if (
  nrow(fcs_pmm_summary_qc) == 0L ||
    !ipw_cap_all_true(fcs_pmm_summary_qc$passed) ||
    !ipw_cap_all_true(fcs_pmm_summary_qc$overall_summary_pass)
) {
  stop("The FCS-PMM summary QC file does not record a full pass.")
}

if (!any(grepl(
  "IPW WEIGHT CAPPING CANARY: PASS",
  readLines(canary_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The IPW weight-capping canary PASS marker is invalid.")
}

canary_qc <- readr::read_csv(canary_qc_file, show_col_types = FALSE)
canary_specification <- readr::read_csv(
  canary_specification_file,
  show_col_types = FALSE
)

if (
  nrow(canary_qc) == 0L ||
    !ipw_cap_all_true(canary_qc$passed) ||
    !ipw_cap_all_true(canary_qc$overall_canary_pass)
) {
  stop("The IPW-capping canary QC file does not record a full pass.")
}

current_specification <- data.frame(
  lower_probability = ipw_cap_lower_probability,
  upper_probability = ipw_cap_upper_probability,
  quantile_type = ipw_cap_quantile_type,
  threshold_population = "estimated raw weights among Phase-2 observations",
  operation = "two-sided winsorization",
  stringsAsFactors = FALSE
)

if (!isTRUE(all.equal(
  current_specification,
  canary_specification,
  tolerance = 0,
  check.attributes = FALSE
))) {
  stop("The current IPW cap specification differs from the validated canary.")
}

calibration <- readRDS(calibration_rds)
dgm <- readRDS(dgm_file)
dgm_qc <- validate_ipw_cap_primary_dgm(dgm, scenario_row)

if (!ipw_cap_all_true(dgm_qc$dgm_matches_exactly)) {
  stop("The matched primary DGM cache failed IPW-capping validation.")
}

final_file <- file.path(
  scenario_chunk_dir,
  sprintf(
    "ipw_cap_scenario_%02d_chunk_%02d.rds",
    sensitivity_scenario_id,
    chunk_id
  )
)

checkpoint_file <- file.path(
  scenario_chunk_dir,
  sprintf(
    "ipw_cap_scenario_%02d_chunk_%02d_checkpoint.rds",
    sensitivity_scenario_id,
    chunk_id
  )
)

if (file.exists(final_file)) {

  existing_final <- tryCatch(readRDS(final_file), error = function(e) NULL)

  required_final_components <- c(
    "metadata", "results", "rng_states", "dataset_qc"
  )
  expected_repetitions <- seq.int(rep_start, rep_end)
  existing_has_components <-
    !is.null(existing_final) &&
    all(required_final_components %in% names(existing_final)) &&
    is.list(existing_final$metadata) &&
    is.data.frame(existing_final$results) &&
    is.data.frame(existing_final$rng_states) &&
    is.data.frame(existing_final$dataset_qc)

  existing_final_matches <-
    existing_has_components &&
    identical(
      existing_final$metadata$sensitivity_scenario_id,
      sensitivity_scenario_id
    ) &&
    identical(existing_final$metadata$chunk_id, chunk_id) &&
    identical(existing_final$metadata$rep_start, rep_start) &&
    identical(existing_final$metadata$rep_end, rep_end) &&
    identical(existing_final$metadata$nsim_ipw_cap, nsim_ipw_cap) &&
    identical(
      existing_final$metadata$reps_per_ipw_cap_chunk,
      reps_per_ipw_cap_chunk
    ) &&
    identical(
      existing_final$metadata$cap_probabilities,
      c(ipw_cap_lower_probability, ipw_cap_upper_probability)
    ) &&
    identical(
      existing_final$metadata$cap_quantile_type,
      ipw_cap_quantile_type
    ) &&
    nrow(existing_final$results) == length(expected_repetitions) &&
    nrow(existing_final$rng_states) == length(expected_repetitions) &&
    nrow(existing_final$dataset_qc) == length(expected_repetitions) &&
    setequal(existing_final$results$repetition, expected_repetitions) &&
    setequal(existing_final$rng_states$repetition, expected_repetitions) &&
    setequal(existing_final$dataset_qc$repetition, expected_repetitions) &&
    identical(unique(existing_final$results$method), ipw_capped_label)

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

  required_checkpoint_fields <- c(
    "sensitivity_scenario_id", "chunk_id", "rep_start", "rep_end",
    "nsim_ipw_cap", "reps_per_ipw_cap_chunk", "next_rep",
    "results_list", "rng_list", "dataset_qc_list"
  )

  checkpoint_has_fields <-
    !is.null(checkpoint) &&
    all(required_checkpoint_fields %in% names(checkpoint)) &&
    is.list(checkpoint$results_list) &&
    is.list(checkpoint$rng_list) &&
    is.list(checkpoint$dataset_qc_list)

  checkpoint_matches <-
    checkpoint_has_fields &&
    identical(checkpoint$sensitivity_scenario_id, sensitivity_scenario_id) &&
    identical(checkpoint$chunk_id, chunk_id) &&
    identical(checkpoint$rep_start, rep_start) &&
    identical(checkpoint$rep_end, rep_end) &&
    identical(checkpoint$nsim_ipw_cap, nsim_ipw_cap) &&
    identical(
      checkpoint$reps_per_ipw_cap_chunk,
      reps_per_ipw_cap_chunk
    ) &&
    length(checkpoint$next_rep) == 1L &&
    is.numeric(checkpoint$next_rep) &&
    is.finite(checkpoint$next_rep) &&
    checkpoint$next_rep >= rep_start &&
    checkpoint$next_rep <= rep_end + 1L &&
    checkpoint$next_rep == floor(checkpoint$next_rep)

  if (!checkpoint_matches) {
    stop(
      "Existing checkpoint does not match the requested task:\n",
      checkpoint_file
    )
  }

  results_list <- checkpoint$results_list
  rng_list <- checkpoint$rng_list
  dataset_qc_list <- checkpoint$dataset_qc_list
  next_rep <- as.integer(checkpoint$next_rep)

  completed_repetitions <- next_rep - rep_start
  populated_results <- sum(!vapply(results_list, is.null, logical(1)))
  populated_rng <- sum(!vapply(rng_list, is.null, logical(1)))
  populated_qc <- sum(!vapply(dataset_qc_list, is.null, logical(1)))

  if (
    populated_results != completed_repetitions ||
      populated_rng != completed_repetitions ||
      populated_qc != completed_repetitions
  ) {
    stop("The checkpoint contents do not agree with its next repetition.")
  }

  cat(
    "Resuming IPW-capping scenario ",
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
  "Task ID: ", task_id, " / ", n_ipw_cap_array_tasks, "\n",
  "Sensitivity scenario: ", sensitivity_scenario_id, " / 4\n",
  "Primary scenario ID: ", scenario_row$primary_scenario_id, "\n",
  "Scenario: ", scenario_row$scenario_label, "\n",
  "Chunk: ", chunk_id, " / ", chunks_per_ipw_cap_scenario, "\n",
  "Repetitions: ", rep_start, "-", rep_end, "\n",
  "N: ", scenario_row$N, "\n",
  "Phase-2 fraction: ", scenario_row$phase2_fraction, "\n",
  "R2 A~M|X: ", scenario_row$r2_a_marker, "\n",
  "R2 Y~M|A,X: ", scenario_row$r2_y_marker, "\n",
  "Theta: ", scenario_row$theta, "\n",
  "Selection: MAR\n",
  "Marker error: MVN\n",
  "Method: IPW with 1st/99th-percentile weight winsorization\n",
  "Quantile type: 7\n",
  "============================================================\n",
  sep = ""
)

if (next_rep <= rep_end) {

  for (repetition in seq.int(next_rep, rep_end)) {

    states <- ipw_cap_get_repetition_states(
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

    structural_qc <- ipw_cap_dataset_qc(
      dat = dat,
      marker_names = calibration$marker_names,
      row = scenario_row
    )

    structural_checks <- unlist(
      structural_qc[
        c(
          "marker_blockwise",
          "phase2_matches_markers",
          "finite_valid_probabilities",
          "sufficient_phase2_n"
        )
      ],
      use.names = FALSE
    )

    if (!ipw_cap_all_true(structural_checks)) {
      stop(
        "Structural DGM QC failed for IPW-capping scenario ",
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

    ans <- fit_ipw_capped_sensitivity(
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
        marker_error = dgm$marker_error
      )

    results_list[[list_index]] <- ans

    rng_list[[list_index]] <- data.frame(
      sensitivity_scenario_id = sensitivity_scenario_id,
      primary_scenario_id = scenario_row$primary_scenario_id,
      repetition = repetition,
      data_rng_state = ipw_cap_state_to_string(states$data_state),
      method_rng_state = ipw_cap_state_to_string(states$method_state),
      stringsAsFactors = FALSE
    )

    checkpoint <- list(
      task_id = task_id,
      sensitivity_scenario_id = sensitivity_scenario_id,
      primary_scenario_id = scenario_row$primary_scenario_id,
      chunk_id = chunk_id,
      rep_start = rep_start,
      rep_end = rep_end,
      nsim_ipw_cap = nsim_ipw_cap,
      reps_per_ipw_cap_chunk = reps_per_ipw_cap_chunk,
      next_rep = repetition + 1L,
      results_list = results_list,
      rng_list = rng_list,
      dataset_qc_list = dataset_qc_list
    )

    ipw_cap_atomic_save_rds(
      checkpoint,
      checkpoint_file,
      compress = FALSE
    )

    cat(
      "IPW-capping scenario ",
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

dataset_keys <- dataset_qc %>%
  count(sensitivity_scenario_id, repetition, name = "n_rows")

rng_keys <- rng_states %>%
  count(sensitivity_scenario_id, repetition, name = "n_rows")

successful_results <- results %>%
  filter(status == "ok")

chunk_checks <- c(
  setequal(unique(results$repetition), expected_repetitions),
  setequal(unique(dataset_qc$repetition), expected_repetitions),
  setequal(unique(rng_states$repetition), expected_repetitions),
  nrow(results) == length(expected_repetitions),
  nrow(dataset_qc) == length(expected_repetitions),
  nrow(rng_states) == length(expected_repetitions),
  nrow(result_keys) == nrow(results) &&
    ipw_cap_all_true(result_keys$n_rows == 1L),
  nrow(dataset_keys) == nrow(dataset_qc) &&
    ipw_cap_all_true(dataset_keys$n_rows == 1L),
  nrow(rng_keys) == nrow(rng_states) &&
    ipw_cap_all_true(rng_keys$n_rows == 1L),
  identical(unique(results$method), ipw_capped_label),
  ipw_cap_all_true(results$selection == "mar"),
  ipw_cap_all_true(results$marker_error == "mvn"),
  ipw_cap_all_true(dataset_qc$marker_blockwise),
  ipw_cap_all_true(dataset_qc$phase2_matches_markers),
  ipw_cap_all_true(dataset_qc$finite_valid_probabilities),
  ipw_cap_all_true(dataset_qc$sufficient_phase2_n),
  nrow(successful_results) == 0L || ipw_cap_all_true(
    successful_results$cap_lower_probability ==
      ipw_cap_lower_probability &
      successful_results$cap_upper_probability ==
        ipw_cap_upper_probability &
      successful_results$cap_quantile_type == ipw_cap_quantile_type
  ),
  nrow(successful_results) == 0L || ipw_cap_all_true(
    successful_results$capped_weight_min >=
      successful_results$cap_lower_threshold - 1e-12 &
      successful_results$capped_weight_max <=
        successful_results$cap_upper_threshold + 1e-12
  )
)

if (!ipw_cap_all_true(chunk_checks)) {
  stop("Completed IPW-capping chunk failed final structural QC.")
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
  nsim_ipw_cap = nsim_ipw_cap,
  reps_per_ipw_cap_chunk = reps_per_ipw_cap_chunk,
  chunks_per_ipw_cap_scenario = chunks_per_ipw_cap_scenario,
  N = scenario_row$N,
  phase2_fraction = scenario_row$phase2_fraction,
  r2_a_marker = scenario_row$r2_a_marker,
  r2_y_marker = scenario_row$r2_y_marker,
  theta = scenario_row$theta,
  selection = dgm$selection,
  marker_error = dgm$marker_error,
  cap_probabilities = c(
    ipw_cap_lower_probability,
    ipw_cap_upper_probability
  ),
  cap_quantile_type = ipw_cap_quantile_type,
  cap_threshold_population =
    "estimated raw weights among Phase-2 observations",
  cap_operation = "two-sided winsorization",
  primary_rng_seed_master = rng_seed_master,
  calibration_md5 = unname(tools::md5sum(calibration_rds)),
  primary_dgm_md5 = unname(tools::md5sum(dgm_file)),
  created = as.character(Sys.time())
)

ipw_cap_atomic_save_rds(
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
  "\nIPW weight-capping production chunk complete:\n",
  final_file,
  "\n",
  sep = ""
)
