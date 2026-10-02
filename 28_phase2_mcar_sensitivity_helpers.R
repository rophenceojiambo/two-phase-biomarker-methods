################################################################################
# 28_phase2_mcar_sensitivity_helpers.R
#
# Shared definitions for the Phase-2 MCAR sensitivity.
#
# This is a thin extension of the validated primary simulation framework. It
# does not modify 00_config.R, 02_dgm.R, or 03_methods.R. The only DGM change is
# selection = "mcar" instead of selection = "mar"; marker_error remains "mvn".
################################################################################


# ------------------------------------------------------------------------------
# 1. Targeted sensitivity scenarios
# ------------------------------------------------------------------------------

mcar_sensitivity_grid <- data.frame(
  sensitivity_scenario_id = 1:4,
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
  stringsAsFactors = FALSE
)


match_mcar_primary_scenario_id <- function(
  row,
  primary_grid,
  tolerance = 1e-12
) {

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


mcar_sensitivity_grid$primary_scenario_id <- vapply(
  seq_len(nrow(mcar_sensitivity_grid)),
  function(i) {
    match_mcar_primary_scenario_id(
      mcar_sensitivity_grid[i, , drop = FALSE],
      primary_grid
    )
  },
  integer(1)
)

mcar_sensitivity_grid <- mcar_sensitivity_grid[
  ,
  c(
    "sensitivity_scenario_id",
    "primary_scenario_id",
    "scenario_label",
    "scenario_key",
    "N",
    "phase2_fraction",
    "r2_a_marker",
    "r2_y_marker",
    "theta"
  )
]

stopifnot(
  identical(
    mcar_sensitivity_grid$primary_scenario_id,
    c(60L, 168L, 100L, 208L)
  )
)

# The computational canary covers both sample-size/Phase-2-fraction regimes.
# Fast large-sample DGM validation is performed for all four scenarios.
mcar_canary_grid <- mcar_sensitivity_grid[
  mcar_sensitivity_grid$theta == 0.15,
  ,
  drop = FALSE
]

stopifnot(nrow(mcar_canary_grid) == 2L)


# Return TRUE only for a nonempty logical vector containing no FALSE or NA.
mcar_all_true <- function(x) {
  length(x) > 0L && isTRUE(all(!is.na(x) & x))
}


# ------------------------------------------------------------------------------
# 2. Common readers, paths, and atomic writers
# ------------------------------------------------------------------------------

read_mcar_positive_integer_env <- function(name, default) {

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


mcar_primary_cache_filename <- function(row) {

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


mcar_dgm_cache_filename <- function(row, dgm_directory) {

  file.path(
    dgm_directory,
    sprintf(
      "mcar_dgm_%02d_%s.rds",
      row$sensitivity_scenario_id,
      row$scenario_key
    )
  )
}


mcar_atomic_save_rds <- function(object, path, compress = FALSE) {

  temporary_path <- paste0(path, ".tmp_", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)

  saveRDS(object, temporary_path, compress = compress)

  if (!file.rename(temporary_path, path)) {
    stop("Could not atomically publish RDS file: ", path)
  }

  invisible(path)
}


mcar_atomic_write_csv <- function(object, path) {

  temporary_path <- paste0(path, ".tmp_", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)

  readr::write_csv(object, temporary_path)

  if (!file.rename(temporary_path, path)) {
    stop("Could not atomically publish CSV file: ", path)
  }

  invisible(path)
}


mcar_atomic_write_lines <- function(text, path) {

  temporary_path <- paste0(path, ".tmp_", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)

  writeLines(text, temporary_path)

  if (!file.rename(temporary_path, path)) {
    stop("Could not atomically publish text file: ", path)
  }

  invisible(path)
}


# Return the largest finite absolute value, or NA if none are finite.
mcar_safe_max_abs <- function(x) {

  values <- abs(x[is.finite(x)])

  if (length(values) == 0L) {
    return(NA_real_)
  }

  max(values)
}


# ------------------------------------------------------------------------------
# 3. Primary-paired L'Ecuyer-CMRG streams
# ------------------------------------------------------------------------------

# Reusing the primary scenario/repetition stream makes the full Phase-1 data
# identical to the matched primary simulation. Only the subsequent Phase-2
# selection draw changes because the DGM selection mechanism is MCAR.

mcar_get_primary_scenario_stream <- function(primary_scenario_id) {

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


mcar_get_repetition_states <- function(primary_scenario_id, repetition) {

  if (
    length(repetition) != 1L ||
      is.na(repetition) ||
      repetition < 1L
  ) {
    stop("repetition must be a positive integer.")
  }

  state <- mcar_get_primary_scenario_stream(primary_scenario_id)
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


mcar_state_to_string <- function(state) {
  paste(state, collapse = ",")
}


# ------------------------------------------------------------------------------
# 4. DGM construction and exact comparison with the primary DGM
# ------------------------------------------------------------------------------

build_mcar_dgm <- function(calibration, row) {

  dgm <- build_dgm_parameters(
    calibration = calibration,
    r2_a_marker = row$r2_a_marker,
    r2_y_marker = row$r2_y_marker,
    theta = row$theta,
    phase2_fraction = row$phase2_fraction,
    selection = "mcar",
    marker_error = "mvn"
  )

  if (!identical(dgm$selection, "mcar")) {
    stop("DGM was not built with selection = 'mcar'.")
  }

  if (!identical(dgm$marker_error, "mvn")) {
    stop("MCAR sensitivity must retain marker_error = 'mvn'.")
  }

  dgm
}


compare_mcar_with_primary_dgm <- function(
  row,
  mcar_dgm,
  primary_dgm,
  primary_cache_file,
  mcar_cache_file = NA_character_
) {

  common_design_components <- c(
    "r2_a_marker",
    "r2_y_marker",
    "theta",
    "phase2_fraction",
    "marker_error"
  )

  data.frame(
    sensitivity_scenario_id = row$sensitivity_scenario_id,
    primary_scenario_id = row$primary_scenario_id,
    scenario_label = row$scenario_label,
    scenario_key = row$scenario_key,
    N = row$N,
    phase2_fraction = row$phase2_fraction,
    r2_a_marker = row$r2_a_marker,
    r2_y_marker = row$r2_y_marker,
    theta = row$theta,
    primary_cache_file = basename(primary_cache_file),
    mcar_cache_file = if (is.na(mcar_cache_file)) {
      ""
    } else {
      basename(mcar_cache_file)
    },
    exposure_parameters_identical = isTRUE(
      all.equal(
        mcar_dgm$exposure,
        primary_dgm$exposure,
        tolerance = 0,
        check.attributes = TRUE
      )
    ),
    outcome_parameters_identical = isTRUE(
      all.equal(
        mcar_dgm$outcome,
        primary_dgm$outcome,
        tolerance = 0,
        check.attributes = TRUE
      )
    ),
    common_design_values_identical = isTRUE(
      all.equal(
        unlist(mcar_dgm[common_design_components]),
        unlist(primary_dgm[common_design_components]),
        tolerance = 0,
        check.attributes = TRUE
      )
    ),
    primary_selection = primary_dgm$selection,
    sensitivity_selection = mcar_dgm$selection,
    selection_changed_as_intended =
      identical(primary_dgm$selection, "mar") &&
      identical(mcar_dgm$selection, "mcar"),
    primary_marker_error = primary_dgm$marker_error,
    sensitivity_marker_error = mcar_dgm$marker_error,
    marker_error_unchanged =
      identical(primary_dgm$marker_error, "mvn") &&
      identical(mcar_dgm$marker_error, "mvn"),
    primary_selection_intercept =
      primary_dgm$selection_parameters$eta0,
    mcar_selection_intercept =
      mcar_dgm$selection_parameters$eta0,
    mcar_intercept_not_required = is.na(
      mcar_dgm$selection_parameters$eta0
    ),
    mcar_expected_fraction_exact = abs(
      mcar_dgm$selection_parameters$expected_fraction_check -
        row$phase2_fraction
    ) <= 1e-15,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------------------
# 5. Large-sample DGM and MCAR-independence validation
# ------------------------------------------------------------------------------

validate_mcar_dgm <- function(
  calibration,
  dgm,
  row,
  validation_n,
  validation_seed
) {

  validation <- check_dgm(
    calibration = calibration,
    dgm_parameters = dgm,
    N_check = validation_n,
    seed = validation_seed
  )

  validation$selection_target_check <-
    dgm$selection_parameters$expected_fraction_check

  validation$dgm_selection <- dgm$selection
  validation$dgm_marker_error <- dgm$marker_error
  validation$validation_seed <- validation_seed

  validation <- dplyr::bind_cols(
    row[rep(1, nrow(validation)), , drop = FALSE],
    validation
  ) %>%
    dplyr::mutate(
      abs_error_r2_A = abs(
        observed_r2_A_M_given_X - target_r2_A_M_given_X
      ),
      abs_error_r2_Y = abs(
        observed_r2_Y_M_given_AX - target_r2_Y_M_given_AX
      ),
      abs_error_theta = abs(observed_theta_full_model - target_theta),
      abs_error_phase2_fraction = abs(
        realized_phase2_fraction - target_phase2_fraction
      ),
      abs_error_selection_target = abs(
        selection_target_check - target_phase2_fraction
      ),
      pass_r2_A = abs_error_r2_A <= 0.01,
      pass_r2_Y = abs_error_r2_Y <= 0.01,
      pass_theta = abs_error_theta <= 0.01,
      pass_mean_A = abs(mean_A) <= 0.02,
      pass_sd_A = abs(sd_A - 1) <= 0.02,
      pass_mean_Y = abs(mean_Y) <= 0.02,
      pass_sd_Y = abs(sd_Y - 1) <= 0.02,
      pass_phase2_fraction = abs_error_phase2_fraction <= 0.01,
      pass_selection_target = abs_error_selection_target <= 1e-15,
      pass_constant_probability =
        abs(min_pi - target_phase2_fraction) <= 1e-15 &
        abs(p01_pi - target_phase2_fraction) <= 1e-15 &
        abs(median_pi - target_phase2_fraction) <= 1e-15 &
        abs(p99_pi - target_phase2_fraction) <= 1e-15 &
        abs(max_pi - target_phase2_fraction) <= 1e-15,
      pass_selection = dgm_selection == "mcar",
      pass_marker_error = dgm_marker_error == "mvn",
      pass_all =
        pass_r2_A &
        pass_r2_Y &
        pass_theta &
        pass_mean_A &
        pass_sd_A &
        pass_mean_Y &
        pass_sd_Y &
        pass_phase2_fraction &
        pass_selection_target &
        pass_constant_probability &
        pass_selection &
        pass_marker_error
    )

  # Independent verification using a separately generated large dataset.
  RNGkind("L'Ecuyer-CMRG")
  set.seed(validation_seed + 1L)

  simulated <- generate_two_phase_data(
    N = validation_n,
    calibration = calibration,
    dgm_parameters = dgm
  )

  dat <- simulated$observed
  Z <- selection_matrix(dat, include_y = TRUE)

  selection_fit <- glm.fit(
    x = Z,
    y = dat$phase2,
    family = binomial()
  )

  nonintercept_columns <- seq.int(2L, ncol(Z))

  selection_correlations <- vapply(
    nonintercept_columns,
    function(j) {
      suppressWarnings(cor(dat$phase2, Z[, j]))
    },
    numeric(1)
  )

  # Compute these summaries before data.frame(). Arguments to data.frame() are
  # evaluated in the calling environment, so a later column cannot refer to a
  # column defined earlier in the same data.frame() call.
  max_abs_pi_deviation_value <- max(
    abs(dat$pi_true - row$phase2_fraction)
  )

  max_abs_selection_correlation_value <- mcar_safe_max_abs(
    selection_correlations
  )

  max_abs_fitted_selection_slope_value <- mcar_safe_max_abs(
    selection_fit$coefficients[nonintercept_columns]
  )

  independence <- data.frame(
    sensitivity_scenario_id = row$sensitivity_scenario_id,
    primary_scenario_id = row$primary_scenario_id,
    scenario_label = row$scenario_label,
    scenario_key = row$scenario_key,
    validation_n = validation_n,
    independence_seed = validation_seed + 1L,
    target_phase2_fraction = row$phase2_fraction,
    realized_phase2_fraction = mean(dat$phase2),
    min_pi_true = min(dat$pi_true),
    max_pi_true = max(dat$pi_true),
    max_abs_pi_deviation = max_abs_pi_deviation_value,
    max_abs_selection_correlation =
      max_abs_selection_correlation_value,
    max_abs_fitted_selection_slope =
      max_abs_fitted_selection_slope_value,
    pass_constant_probability =
      max_abs_pi_deviation_value <= 1e-15,
    pass_empirical_independence =
      isTRUE(selection_fit$converged) &
      is.finite(max_abs_selection_correlation_value) &
      is.finite(max_abs_fitted_selection_slope_value) &
      max_abs_selection_correlation_value <= 0.03 &
      max_abs_fitted_selection_slope_value <= 0.08,
    stringsAsFactors = FALSE
  )

  independence$pass_all <-
    independence$pass_constant_probability &
    independence$pass_empirical_independence

  list(
    dgm_validation = validation,
    independence_validation = independence
  )
}


# ------------------------------------------------------------------------------
# 6. Timed execution of the six unchanged primary methods
# ------------------------------------------------------------------------------

mcar_run_timed <- function(method_name, expr) {

  start_time <- proc.time()[["elapsed"]]

  ans <- tryCatch(
    eval.parent(substitute(expr)),
    error = function(e) failed_result(method_name, conditionMessage(e))
  )

  ans$elapsed_seconds <- proc.time()[["elapsed"]] - start_time
  ans
}


run_all_mcar_methods_timed <- function(dat, marker_names) {

  dplyr::bind_rows(
    mcar_run_timed("Naive", fit_naive(dat)),
    mcar_run_timed("CCA", fit_cca(dat, marker_names)),
    mcar_run_timed(
      "FCS-MI",
      fit_fcs_mi(
        dat = dat,
        marker_names = marker_names,
        nimp = nimp_primary,
        maxit = mice_maxit_primary,
        method = mice_method_primary
      )
    ),
    mcar_run_timed(
      "JM-MI",
      fit_joint_mi(
        dat = dat,
        marker_names = marker_names,
        nimp = nimp_primary,
        nburn = jomo_nburn_primary,
        nbetween = jomo_nbetween_primary
      )
    ),
    mcar_run_timed("IPW", fit_ipw(dat, marker_names)),
    mcar_run_timed(
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
}
