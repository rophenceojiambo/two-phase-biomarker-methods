################################################################################
# 34_fcs_pmm_sensitivity_helpers.R
#
# Shared definitions for the FCS predictive-mean-matching sensitivity.
#
# The validated primary MAR/MVN data-generating mechanisms are reused without
# modification. The sole analysis change is method = "pmm" instead of "norm"
# for the eight marker equations in fit_fcs_mi().
################################################################################


# ------------------------------------------------------------------------------
# 1. Targeted scenarios and method labels
# ------------------------------------------------------------------------------

fcs_pmm_sensitivity_grid <- data.frame(
  sensitivity_scenario_id = 1:4,
  primary_scenario_id = c(60L, 168L, 100L, 208L),
  scenario_label = c(
    "MIDUS-like null",
    "MIDUS-like alternative",
    "Stress null",
    "Stress alternative"
  ),
  scenario_key = c(
    "midus_null",
    "midus_alternative",
    "stress_null",
    "stress_alternative"
  ),
  N = c(800L, 800L, 500L, 500L),
  phase2_fraction = c(0.65, 0.65, 0.25, 0.25),
  r2_a_marker = c(0.05, 0.05, 0.10, 0.10),
  r2_y_marker = c(0.05, 0.05, 0.10, 0.10),
  theta = c(0, 0.15, 0, 0.15),
  selection = "mar",
  marker_error = "mvn",
  stringsAsFactors = FALSE
)


match_fcs_pmm_primary_scenario_id <- function(
  row,
  primary_grid,
  tolerance = 1e-12
) {

  required_columns <- c(
    "N",
    "phase2_fraction",
    "r2_a_marker",
    "r2_y_marker",
    "theta"
  )

  if (nrow(row) != 1L) {
    stop("row must contain exactly one sensitivity scenario.")
  }

  if (!all(required_columns %in% names(row))) {
    stop("row does not contain all required scenario variables.")
  }

  hit <- primary_grid$N == row$N &
    abs(primary_grid$phase2_fraction - row$phase2_fraction) < tolerance &
    abs(primary_grid$r2_a_marker - row$r2_a_marker) < tolerance &
    abs(primary_grid$r2_y_marker - row$r2_y_marker) < tolerance &
    abs(primary_grid$theta - row$theta) < tolerance

  ids <- primary_grid$scenario_id[hit]

  if (length(ids) != 1L) {
    stop(
      "Could not identify exactly one matching primary scenario for ",
      row$scenario_label,
      "."
    )
  }

  ids[[1]]
}


matched_primary_ids <- vapply(
  seq_len(nrow(fcs_pmm_sensitivity_grid)),
  function(i) {
    match_fcs_pmm_primary_scenario_id(
      fcs_pmm_sensitivity_grid[i, , drop = FALSE],
      primary_grid
    )
  },
  integer(1)
)

stopifnot(
  identical(
    matched_primary_ids,
    fcs_pmm_sensitivity_grid$primary_scenario_id
  )
)

fcs_pmm_canary_grid <- fcs_pmm_sensitivity_grid[
  fcs_pmm_sensitivity_grid$theta == 0.15,
  ,
  drop = FALSE
]

stopifnot(nrow(fcs_pmm_canary_grid) == 2L)

fcs_norm_label <- "FCS-Norm"
fcs_pmm_label <- "FCS-PMM"


# Return TRUE only for a nonempty logical vector containing no FALSE or NA.
fcs_pmm_all_true <- function(x) {
  length(x) > 0L && isTRUE(all(!is.na(x) & x))
}


# ------------------------------------------------------------------------------
# 2. Environment readers, file paths, and atomic writers
# ------------------------------------------------------------------------------

read_fcs_pmm_positive_integer_env <- function(name, default) {

  value_text <- Sys.getenv(name, unset = as.character(default))
  value_numeric <- suppressWarnings(as.numeric(value_text))

  valid <- length(value_numeric) == 1L &&
    is.finite(value_numeric) &&
    value_numeric >= 1 &&
    value_numeric == floor(value_numeric) &&
    value_numeric <= .Machine$integer.max

  if (!valid) {
    stop(name, " must be a positive integer; received '", value_text, "'.")
  }

  as.integer(value_numeric)
}


fcs_pmm_primary_cache_filename <- function(row) {

  file.path(
    dgm_cache_dir,
    sprintf(
      "dgm_f%03d_rA%03d_rY%03d_th%03d.rds",
      round(100 * row$phase2_fraction),
      round(100 * row$r2_a_marker),
      round(100 * row$r2_y_marker),
      round(100 * row$theta)
    )
  )
}


fcs_pmm_atomic_save_rds <- function(object, path, compress = FALSE) {

  temporary_path <- paste0(path, ".tmp_", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)

  saveRDS(object, temporary_path, compress = compress)

  if (!file.rename(temporary_path, path)) {
    stop("Could not atomically publish RDS file: ", path)
  }

  invisible(path)
}


fcs_pmm_atomic_write_csv <- function(object, path) {

  temporary_path <- paste0(path, ".tmp_", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)

  readr::write_csv(object, temporary_path)

  if (!file.rename(temporary_path, path)) {
    stop("Could not atomically publish CSV file: ", path)
  }

  invisible(path)
}


fcs_pmm_atomic_write_lines <- function(text, path) {

  temporary_path <- paste0(path, ".tmp_", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)

  writeLines(text, temporary_path)

  if (!file.rename(temporary_path, path)) {
    stop("Could not atomically publish text file: ", path)
  }

  invisible(path)
}


# ------------------------------------------------------------------------------
# 3. Primary-paired L'Ecuyer-CMRG streams
# ------------------------------------------------------------------------------

fcs_pmm_get_primary_scenario_stream <- function(primary_scenario_id) {

  if (
    length(primary_scenario_id) != 1L ||
      is.na(primary_scenario_id) ||
      primary_scenario_id < 1L
  ) {
    stop("primary_scenario_id must be a positive integer.")
  }

  RNGkind("L'Ecuyer-CMRG")
  set.seed(rng_seed_master)
  state <- .Random.seed

  if (primary_scenario_id > 1L) {
    for (s in seq_len(primary_scenario_id - 1L)) {
      state <- parallel::nextRNGStream(state)
    }
  }

  state
}


fcs_pmm_get_repetition_states <- function(primary_scenario_id, repetition) {

  if (
    length(repetition) != 1L ||
      is.na(repetition) ||
      repetition < 1L
  ) {
    stop("repetition must be a positive integer.")
  }

  state <- fcs_pmm_get_primary_scenario_stream(primary_scenario_id)
  n_advance <- 2L * (repetition - 1L)

  if (n_advance > 0L) {
    for (j in seq_len(n_advance)) {
      state <- parallel::nextRNGSubStream(state)
    }
  }

  data_state <- state
  method_state <- parallel::nextRNGSubStream(data_state)

  list(
    data_state = data_state,
    method_state = method_state,
    next_state = parallel::nextRNGSubStream(method_state)
  )
}


fcs_pmm_state_to_string <- function(state) {
  paste(state, collapse = ",")
}


# ------------------------------------------------------------------------------
# 4. Exact primary-DGM and PMM-default validation
# ------------------------------------------------------------------------------

validate_fcs_pmm_primary_dgm <- function(dgm, row) {

  data.frame(
    sensitivity_scenario_id = row$sensitivity_scenario_id,
    primary_scenario_id = row$primary_scenario_id,
    scenario_label = row$scenario_label,
    cache_file = basename(fcs_pmm_primary_cache_filename(row)),
    N = row$N,
    expected_selection = row$selection,
    observed_selection = dgm$selection,
    expected_marker_error = row$marker_error,
    observed_marker_error = dgm$marker_error,
    expected_phase2_fraction = row$phase2_fraction,
    observed_phase2_fraction = dgm$phase2_fraction,
    expected_r2_a_marker = row$r2_a_marker,
    observed_r2_a_marker = dgm$r2_a_marker,
    expected_r2_y_marker = row$r2_y_marker,
    observed_r2_y_marker = dgm$r2_y_marker,
    expected_theta = row$theta,
    observed_theta = dgm$theta,
    dgm_matches_exactly =
      identical(dgm$selection, "mar") &&
      identical(dgm$marker_error, "mvn") &&
      isTRUE(all.equal(dgm$phase2_fraction, row$phase2_fraction, tolerance = 0)) &&
      isTRUE(all.equal(dgm$r2_a_marker, row$r2_a_marker, tolerance = 0)) &&
      isTRUE(all.equal(dgm$r2_y_marker, row$r2_y_marker, tolerance = 0)) &&
      isTRUE(all.equal(dgm$theta, row$theta, tolerance = 0)),
    stringsAsFactors = FALSE
  )
}


fcs_pmm_default_settings <- function() {

  pmm_function <- mice::mice.impute.pmm
  defaults <- formals(pmm_function)
  function_environment <- environment(pmm_function)

  evaluate_default <- function(argument_name) {
    eval(defaults[[argument_name]], envir = function_environment)
  }

  settings <- data.frame(
    mice_version = as.character(utils::packageVersion("mice")),
    donors = as.integer(evaluate_default("donors")),
    matchtype = as.integer(evaluate_default("matchtype")),
    ridge = as.numeric(evaluate_default("ridge")),
    expected_donors = 5L,
    expected_matchtype = 1L,
    expected_ridge = 1e-5,
    defaults_match_prespecified =
      as.integer(evaluate_default("donors")) == 5L &&
      as.integer(evaluate_default("matchtype")) == 1L &&
      isTRUE(all.equal(
        as.numeric(evaluate_default("ridge")),
        1e-5,
        tolerance = 0
      )),
    stringsAsFactors = FALSE
  )

  if (!fcs_pmm_all_true(settings$defaults_match_prespecified)) {
    warning(
      "The installed mice PMM defaults differ from the prespecified values."
    )
  }

  settings
}


# ------------------------------------------------------------------------------
# 5. Thin wrappers around the validated FCS implementation
# ------------------------------------------------------------------------------

fcs_pmm_run_timed <- function(method_label, expr) {

  start_time <- proc.time()[["elapsed"]]

  ans <- tryCatch(
    eval.parent(substitute(expr)),
    error = function(e) failed_result(method_label, conditionMessage(e))
  )

  ans$method <- method_label
  ans$elapsed_seconds <- proc.time()[["elapsed"]] - start_time
  ans
}


fit_fcs_norm_reference <- function(dat, marker_names) {

  fcs_pmm_run_timed(
    fcs_norm_label,
    fit_fcs_mi(
      dat = dat,
      marker_names = marker_names,
      nimp = nimp_primary,
      maxit = mice_maxit_primary,
      method = "norm"
    )
  )
}


fit_fcs_pmm_sensitivity <- function(dat, marker_names) {

  fcs_pmm_run_timed(
    fcs_pmm_label,
    fit_fcs_mi(
      dat = dat,
      marker_names = marker_names,
      nimp = nimp_primary,
      maxit = mice_maxit_primary,
      method = "pmm"
    )
  )
}


fcs_pmm_dataset_qc <- function(dat, marker_names, row) {

  n_marker_observed <- rowSums(
    !is.na(dat[, marker_names, drop = FALSE])
  )

  data.frame(
    marker_blockwise = all(
      n_marker_observed %in% c(0L, length(marker_names))
    ),
    phase2_matches_markers = identical(
      as.integer(n_marker_observed == length(marker_names)),
      as.integer(dat$phase2)
    ),
    n_phase2 = sum(dat$phase2),
    realized_phase2_fraction = mean(dat$phase2),
    min_pi_true = min(dat$pi_true),
    max_pi_true = max(dat$pi_true),
    finite_valid_probabilities =
      fcs_pmm_all_true(is.finite(dat$pi_true)) &&
      fcs_pmm_all_true(dat$pi_true > 0) &&
      fcs_pmm_all_true(dat$pi_true < 1),
    sufficient_phase2_n =
      sum(dat$phase2) > length(marker_names) + 10L,
    target_phase2_fraction = row$phase2_fraction,
    stringsAsFactors = FALSE
  )
}
