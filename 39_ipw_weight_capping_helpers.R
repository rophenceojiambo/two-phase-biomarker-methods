################################################################################
# 39_ipw_weight_capping_helpers.R
#
# Shared definitions for the IPW 1st/99th-percentile weight-capping sensitivity.
#
# The validated primary MAR/MVN data-generating mechanisms and selection model
# are unchanged. Raw inverse-probability weights among Phase-2 observations are
# winsorized within each simulated dataset at their empirical 1st and 99th
# percentiles (R quantile type 7).
################################################################################


# ------------------------------------------------------------------------------
# 1. Targeted scenarios and prespecified cap
# ------------------------------------------------------------------------------

ipw_cap_sensitivity_grid <- data.frame(
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

ipw_cap_lower_probability <- 0.01
ipw_cap_upper_probability <- 0.99
ipw_cap_quantile_type <- 7L

ipw_uncapped_label <- "IPW-Uncapped"
ipw_capped_label <- "IPW-Capped-1-99"


match_ipw_cap_primary_scenario_id <- function(
  row,
  primary_grid,
  tolerance = 1e-12
) {

  row_columns <- c(
    "N", "phase2_fraction", "r2_a_marker", "r2_y_marker", "theta",
    "scenario_label"
  )
  grid_columns <- c(
    "scenario_id", "N", "phase2_fraction", "r2_a_marker",
    "r2_y_marker", "theta"
  )

  if (nrow(row) != 1L || !all(row_columns %in% names(row))) {
    stop("row must contain exactly one complete sensitivity scenario.")
  }

  if (!all(grid_columns %in% names(primary_grid))) {
    stop("primary_grid is missing columns required for scenario matching.")
  }

  if (
    length(tolerance) != 1L ||
      !is.finite(tolerance) ||
      tolerance < 0
  ) {
    stop("tolerance must be one finite, non-negative number.")
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
  seq_len(nrow(ipw_cap_sensitivity_grid)),
  function(i) {
    match_ipw_cap_primary_scenario_id(
      ipw_cap_sensitivity_grid[i, , drop = FALSE],
      primary_grid
    )
  },
  integer(1)
)

stopifnot(
  identical(
    matched_primary_ids,
    ipw_cap_sensitivity_grid$primary_scenario_id
  )
)

ipw_cap_canary_grid <- ipw_cap_sensitivity_grid[
  ipw_cap_sensitivity_grid$theta == 0.15,
  ,
  drop = FALSE
]

stopifnot(nrow(ipw_cap_canary_grid) == 2L)


# ------------------------------------------------------------------------------
# 2. Environment readers, paths, and atomic writers
# ------------------------------------------------------------------------------

ipw_cap_all_true <- function(x) {
  length(x) > 0L && !anyNA(x) && all(x)
}


read_ipw_cap_positive_integer_env <- function(name, default) {

  value_text <- Sys.getenv(name, unset = as.character(default))
  value_numeric <- suppressWarnings(as.numeric(value_text))

  if (
    length(value_numeric) != 1L ||
      !is.finite(value_numeric) ||
      value_numeric < 1 ||
      value_numeric != floor(value_numeric) ||
      value_numeric > .Machine$integer.max
  ) {
    stop(name, " must be a positive integer; received '", value_text, "'.")
  }

  as.integer(value_numeric)
}


ipw_cap_primary_cache_filename <- function(row) {

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


ipw_cap_atomic_save_rds <- function(object, path, compress = FALSE) {

  temporary_path <- paste0(path, ".tmp_", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)

  saveRDS(object, temporary_path, compress = compress)

  if (!file.rename(temporary_path, path)) {
    stop("Could not atomically publish RDS file: ", path)
  }

  invisible(path)
}


ipw_cap_atomic_write_csv <- function(object, path) {

  temporary_path <- paste0(path, ".tmp_", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)

  readr::write_csv(object, temporary_path)

  if (!file.rename(temporary_path, path)) {
    stop("Could not atomically publish CSV file: ", path)
  }

  invisible(path)
}


ipw_cap_atomic_write_lines <- function(text, path) {

  temporary_path <- paste0(path, ".tmp_", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)

  writeLines(text, temporary_path)

  if (!file.rename(temporary_path, path)) {
    stop("Could not atomically publish text file: ", path)
  }

  invisible(path)
}


# ------------------------------------------------------------------------------
# 3. Exact primary L'Ecuyer-CMRG streams
# ------------------------------------------------------------------------------

ipw_cap_get_primary_scenario_stream <- function(primary_scenario_id) {

  if (
    length(primary_scenario_id) != 1L ||
      !is.finite(primary_scenario_id) ||
      primary_scenario_id < 1 ||
      primary_scenario_id != floor(primary_scenario_id)
  ) {
    stop("primary_scenario_id must be a positive integer.")
  }

  primary_scenario_id <- as.integer(primary_scenario_id)

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


ipw_cap_get_repetition_states <- function(primary_scenario_id, repetition) {

  if (
    length(repetition) != 1L ||
      !is.finite(repetition) ||
      repetition < 1 ||
      repetition != floor(repetition)
  ) {
    stop("repetition must be a positive integer.")
  }

  repetition <- as.integer(repetition)

  state <- ipw_cap_get_primary_scenario_stream(primary_scenario_id)
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


ipw_cap_state_to_string <- function(state) {
  paste(state, collapse = ",")
}


# ------------------------------------------------------------------------------
# 4. Exact primary-DGM validation
# ------------------------------------------------------------------------------

validate_ipw_cap_primary_dgm <- function(dgm, row) {

  required_dgm_components <- c(
    "selection", "marker_error", "phase2_fraction", "r2_a_marker",
    "r2_y_marker", "theta"
  )

  if (!all(required_dgm_components %in% names(dgm))) {
    stop("The DGM object is missing components required for validation.")
  }

  if (nrow(row) != 1L) {
    stop("row must contain exactly one sensitivity scenario.")
  }

  data.frame(
    sensitivity_scenario_id = row$sensitivity_scenario_id,
    primary_scenario_id = row$primary_scenario_id,
    scenario_label = row$scenario_label,
    cache_file = basename(ipw_cap_primary_cache_filename(row)),
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


# ------------------------------------------------------------------------------
# 5. Type-7 quantile and derivative with respect to selection coefficients
# ------------------------------------------------------------------------------

ipw_cap_type7_quantile_derivative <- function(values, derivatives, probability) {

  derivatives <- as.matrix(derivatives)
  n <- length(values)

  if (
    n < 2L ||
      !is.numeric(values) ||
      any(!is.finite(values)) ||
      !is.numeric(derivatives) ||
      nrow(derivatives) != n ||
      ncol(derivatives) < 1L ||
      any(!is.finite(derivatives)) ||
      length(probability) != 1L ||
      !is.finite(probability) ||
      probability < 0 ||
      probability > 1
  ) {
    stop("Invalid inputs to the type-7 quantile derivative helper.")
  }

  ordering <- order(values)
  ordered_values <- values[ordering]
  ordered_derivatives <- derivatives[ordering, , drop = FALSE]

  h <- 1 + (n - 1) * probability
  lower_index <- floor(h)
  interpolation_fraction <- h - lower_index
  upper_index <- min(lower_index + 1L, n)

  value <-
    (1 - interpolation_fraction) * ordered_values[lower_index] +
    interpolation_fraction * ordered_values[upper_index]

  derivative <-
    (1 - interpolation_fraction) *
      ordered_derivatives[lower_index, , drop = FALSE] +
    interpolation_fraction *
      ordered_derivatives[upper_index, , drop = FALSE]

  list(
    value = as.numeric(value),
    derivative = as.numeric(derivative),
    lower_index = lower_index,
    upper_index = upper_index,
    interpolation_fraction = interpolation_fraction
  )
}


ipw_cap_weight_information <- function(
  alpha,
  Z_obs,
  cap_probabilities = c(
    ipw_cap_lower_probability,
    ipw_cap_upper_probability
  )
) {

  Z_obs <- as.matrix(Z_obs)

  if (
    !is.numeric(alpha) ||
      any(!is.finite(alpha)) ||
      !is.numeric(Z_obs) ||
      any(!is.finite(Z_obs)) ||
      nrow(Z_obs) < 2L ||
      ncol(Z_obs) != length(alpha)
  ) {
    stop("alpha and Z_obs are not valid for weight calculation.")
  }

  pi_obs <- as.vector(plogis(Z_obs %*% alpha))

  if (any(!is.finite(pi_obs)) || any(pi_obs <= 0)) {
    stop("Observed selection probabilities are non-finite or zero.")
  }

  raw_weights <- 1 / pi_obs

  raw_weight_derivatives <- Z_obs * (
    -(1 - pi_obs) / pi_obs
  )

  if (is.null(cap_probabilities)) {
    return(list(
      pi_obs = pi_obs,
      raw_weights = raw_weights,
      analysis_weights = raw_weights,
      analysis_weight_derivatives = raw_weight_derivatives,
      lower_threshold = -Inf,
      upper_threshold = Inf,
      lower_probability = NA_real_,
      upper_probability = NA_real_,
      lower_capped = rep(FALSE, length(raw_weights)),
      upper_capped = rep(FALSE, length(raw_weights))
    ))
  }

  if (
    length(cap_probabilities) != 2L ||
      !is.numeric(cap_probabilities) ||
      any(!is.finite(cap_probabilities)) ||
      cap_probabilities[[1]] <= 0 ||
      cap_probabilities[[2]] >= 1 ||
      cap_probabilities[[1]] >= cap_probabilities[[2]]
  ) {
    stop("cap_probabilities must contain valid lower and upper probabilities.")
  }

  lower <- ipw_cap_type7_quantile_derivative(
    values = raw_weights,
    derivatives = raw_weight_derivatives,
    probability = cap_probabilities[[1]]
  )

  upper <- ipw_cap_type7_quantile_derivative(
    values = raw_weights,
    derivatives = raw_weight_derivatives,
    probability = cap_probabilities[[2]]
  )

  lower_capped <- raw_weights < lower$value
  upper_capped <- raw_weights > upper$value

  analysis_weights <- pmin(
    pmax(raw_weights, lower$value),
    upper$value
  )

  analysis_weight_derivatives <- raw_weight_derivatives

  if (any(lower_capped)) {
    analysis_weight_derivatives[lower_capped, ] <- matrix(
      lower$derivative,
      nrow = sum(lower_capped),
      ncol = ncol(Z_obs),
      byrow = TRUE
    )
  }

  if (any(upper_capped)) {
    analysis_weight_derivatives[upper_capped, ] <- matrix(
      upper$derivative,
      nrow = sum(upper_capped),
      ncol = ncol(Z_obs),
      byrow = TRUE
    )
  }

  list(
    pi_obs = pi_obs,
    raw_weights = raw_weights,
    analysis_weights = analysis_weights,
    analysis_weight_derivatives = analysis_weight_derivatives,
    lower_threshold = lower$value,
    upper_threshold = upper$value,
    lower_probability = cap_probabilities[[1]],
    upper_probability = cap_probabilities[[2]],
    lower_capped = lower_capped,
    upper_capped = upper_capped
  )
}


# ------------------------------------------------------------------------------
# 6. IPW estimator with optional percentile capping
# ------------------------------------------------------------------------------

fit_ipw_percentile_cap <- function(
  dat,
  marker_names,
  cap_probabilities = c(
    ipw_cap_lower_probability,
    ipw_cap_upper_probability
  ),
  method_label = ipw_capped_label,
  prob_floor = probability_floor
) {

  required_columns <- c("phase2", "Y", marker_names)

  if (!all(required_columns %in% names(dat))) {
    stop("The analysis dataset is missing variables required for IPW.")
  }

  if (
    length(prob_floor) != 1L ||
      !is.finite(prob_floor) ||
      prob_floor <= 0 ||
      prob_floor >= 0.5
  ) {
    stop("prob_floor must be one finite number between 0 and 0.5.")
  }

  S <- dat$phase2

  if (anyNA(S) || !all(S %in% c(0, 1))) {
    stop("phase2 must be a complete binary indicator.")
  }

  N <- nrow(dat)
  Z <- selection_matrix(dat, include_y = TRUE)

  sel_fit <- glm.fit(
    x = Z,
    y = S,
    family = binomial()
  )

  alpha <- sel_fit$coefficients

  if (any(!is.finite(alpha))) {
    stop("Non-finite coefficient in capped-IPW selection model.")
  }

  pi_hat <- as.vector(plogis(Z %*% alpha))

  if (
    any(!is.finite(pi_hat)) ||
      any(pi_hat < prob_floor | pi_hat > 1 - prob_floor)
  ) {
    stop(
      "Estimated Phase-2 probabilities too close to 0/1 for stable capped IPW."
    )
  }

  obs <- which(S == 1)

  if (length(obs) <= length(marker_names) + 10L) {
    stop("Too few Phase-2 observations for capped IPW.")
  }

  C_obs <- model.matrix(
    full_formula(marker_names),
    data = dat[obs, , drop = FALSE]
  )

  Y_obs <- dat$Y[obs]
  Z_obs <- Z[obs, , drop = FALSE]

  weight_info <- ipw_cap_weight_information(
    alpha = alpha,
    Z_obs = Z_obs,
    cap_probabilities = cap_probabilities
  )

  w <- weight_info$analysis_weights
  raw_w <- weight_info$raw_weights

  XtWX <- crossprod(C_obs, C_obs * w)
  XtWY <- crossprod(C_obs, Y_obs * w)

  beta <- as.vector(solve(XtWX, XtWY))
  names(beta) <- colnames(C_obs)

  q <- ncol(Z)
  d <- ncol(C_obs)

  psi_alpha <- Z * (S - pi_hat)
  psi_beta <- matrix(0, nrow = N, ncol = d)
  colnames(psi_beta) <- colnames(C_obs)

  resid_obs <- Y_obs - as.vector(C_obs %*% beta)
  G_obs <- C_obs * resid_obs
  psi_beta[obs, ] <- G_obs * w

  J_aa <- -crossprod(
    Z,
    Z * (pi_hat * (1 - pi_hat))
  )

  J_bb <- -crossprod(C_obs, C_obs * w)

  J_ba <- crossprod(
    G_obs,
    weight_info$analysis_weight_derivatives
  )

  J <- rbind(
    cbind(J_aa, matrix(0, q, d)),
    cbind(J_ba, J_bb)
  )

  Psi <- cbind(psi_alpha, psi_beta)
  meat <- crossprod(Psi)

  J_inv <- solve(J)
  V_stack <- J_inv %*% meat %*% t(J_inv)

  beta_indices <- q + seq_len(d)
  V_beta <- V_stack[beta_indices, beta_indices, drop = FALSE]

  A_index <- match("A", names(beta))

  if (is.na(A_index)) {
    stop("A coefficient not found in capped-IPW outcome design.")
  }

  est <- beta[A_index]
  se <- sqrt(V_beta[A_index, A_index])

  if (!is.finite(se) || se <= 0) {
    stop("Capped-IPW standard error is non-finite or non-positive.")
  }

  stat <- est / se
  p <- 2 * pnorm(-abs(stat))
  crit <- qnorm(0.975)

  data.frame(
    method = method_label,
    estimate = est,
    se = se,
    df = Inf,
    p_value = p,
    conf_low = est - crit * se,
    conf_high = est + crit * se,
    status = "ok",
    message = NA_character_,
    cap_lower_probability = weight_info$lower_probability,
    cap_upper_probability = weight_info$upper_probability,
    cap_quantile_type = ipw_cap_quantile_type,
    cap_lower_threshold = weight_info$lower_threshold,
    cap_upper_threshold = weight_info$upper_threshold,
    n_capped_lower = sum(weight_info$lower_capped),
    n_capped_upper = sum(weight_info$upper_capped),
    proportion_capped_lower = mean(weight_info$lower_capped),
    proportion_capped_upper = mean(weight_info$upper_capped),
    raw_weight_min = min(raw_w),
    raw_weight_p01 = unname(quantile(
      raw_w,
      ipw_cap_lower_probability,
      type = ipw_cap_quantile_type
    )),
    raw_weight_p99 = unname(quantile(
      raw_w,
      ipw_cap_upper_probability,
      type = ipw_cap_quantile_type
    )),
    raw_weight_max = max(raw_w),
    raw_weight_cv = sd(raw_w) / mean(raw_w),
    raw_weight_ess = sum(raw_w)^2 / sum(raw_w^2),
    capped_weight_min = min(w),
    capped_weight_p01 = unname(quantile(
      w,
      ipw_cap_lower_probability,
      type = ipw_cap_quantile_type
    )),
    capped_weight_p99 = unname(quantile(
      w,
      ipw_cap_upper_probability,
      type = ipw_cap_quantile_type
    )),
    capped_weight_max = max(w),
    capped_weight_cv = sd(w) / mean(w),
    capped_weight_ess = sum(w)^2 / sum(w^2),
    weight_min = min(w),
    weight_p99 = unname(quantile(
      w,
      ipw_cap_upper_probability,
      type = ipw_cap_quantile_type
    )),
    weight_max = max(w),
    weight_cv = sd(w) / mean(w),
    weight_ess = sum(w)^2 / sum(w^2),
    stringsAsFactors = FALSE
  )
}


ipw_cap_run_timed <- function(method_label, expr) {

  start_time <- proc.time()[["elapsed"]]

  ans <- tryCatch(
    eval.parent(substitute(expr)),
    error = function(e) failed_result(method_label, conditionMessage(e))
  )

  ans$method <- method_label
  ans$elapsed_seconds <- proc.time()[["elapsed"]] - start_time
  ans
}


fit_ipw_uncapped_reference <- function(dat, marker_names) {

  ipw_cap_run_timed(
    ipw_uncapped_label,
    fit_ipw_percentile_cap(
      dat = dat,
      marker_names = marker_names,
      cap_probabilities = NULL,
      method_label = ipw_uncapped_label
    )
  )
}


fit_ipw_capped_sensitivity <- function(dat, marker_names) {

  ipw_cap_run_timed(
    ipw_capped_label,
    fit_ipw_percentile_cap(
      dat = dat,
      marker_names = marker_names,
      cap_probabilities = c(
        ipw_cap_lower_probability,
        ipw_cap_upper_probability
      ),
      method_label = ipw_capped_label
    )
  )
}


# ------------------------------------------------------------------------------
# 7. Numerical check of the capped-weight alpha-to-beta Jacobian block
# ------------------------------------------------------------------------------

ipw_cap_jacobian_check <- function(
  dat,
  marker_names,
  relative_step = 1e-6
) {

  S <- dat$phase2
  Z <- selection_matrix(dat, include_y = TRUE)
  obs <- which(S == 1)
  Z_obs <- Z[obs, , drop = FALSE]

  C_obs <- model.matrix(
    full_formula(marker_names),
    data = dat[obs, , drop = FALSE]
  )

  Y_obs <- dat$Y[obs]

  selection_fit <- glm.fit(
    x = Z,
    y = S,
    family = binomial()
  )

  alpha <- selection_fit$coefficients
  weight_info <- ipw_cap_weight_information(alpha, Z_obs)
  w <- weight_info$analysis_weights

  beta <- as.vector(
    solve(
      crossprod(C_obs, C_obs * w),
      crossprod(C_obs, Y_obs * w)
    )
  )

  residuals <- Y_obs - as.vector(C_obs %*% beta)
  G_obs <- C_obs * residuals

  analytic <- crossprod(
    G_obs,
    weight_info$analysis_weight_derivatives
  )

  score_at_alpha <- function(alpha_value) {
    weights <- ipw_cap_weight_information(alpha_value, Z_obs)$analysis_weights
    colSums(G_obs * weights)
  }

  numerical <- matrix(NA_real_, nrow = ncol(C_obs), ncol = ncol(Z))

  for (j in seq_along(alpha)) {
    step <- relative_step * (1 + abs(alpha[[j]]))
    alpha_plus <- alpha
    alpha_minus <- alpha
    alpha_plus[[j]] <- alpha_plus[[j]] + step
    alpha_minus[[j]] <- alpha_minus[[j]] - step

    numerical[, j] <- (
      score_at_alpha(alpha_plus) - score_at_alpha(alpha_minus)
    ) / (2 * step)
  }

  difference <- analytic - numerical
  scale <- pmax(1, abs(analytic), abs(numerical))

  data.frame(
    max_absolute_difference = max(abs(difference)),
    max_scaled_difference = max(abs(difference) / scale),
    analytic_max_absolute = max(abs(analytic)),
    numerical_max_absolute = max(abs(numerical)),
    relative_step = relative_step,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------------------------
# 8. Dataset QC
# ------------------------------------------------------------------------------

ipw_cap_dataset_qc <- function(dat, marker_names, row) {

  required_columns <- c("phase2", "pi_true", marker_names)

  if (!all(required_columns %in% names(dat))) {
    stop("The generated dataset is missing variables required for QC.")
  }

  if (nrow(row) != 1L || !("phase2_fraction" %in% names(row))) {
    stop("row must contain one scenario with phase2_fraction.")
  }

  n_marker_observed <- rowSums(
    !is.na(dat[, marker_names, drop = FALSE])
  )

  data.frame(
    marker_blockwise = ipw_cap_all_true(
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
      ipw_cap_all_true(is.finite(dat$pi_true)) &&
      ipw_cap_all_true(dat$pi_true > 0) &&
      ipw_cap_all_true(dat$pi_true < 1),
    sufficient_phase2_n =
      sum(dat$phase2) > length(marker_names) + 10L,
    target_phase2_fraction = row$phase2_fraction,
    stringsAsFactors = FALSE
  )
}
