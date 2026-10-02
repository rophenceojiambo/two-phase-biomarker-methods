################################################################################
# 16_aipw_strong_misspec_helpers.R
#
# Defines the strong-misspecification AIPW sensitivity analysis.
#
# The primary robustness analysis omitted Y from the misspecified nuisance
# models. This sensitivity omits both A and Y:
#
#   Correct selection model:       S | A, Y, X
#   Misspecified selection model:  S | X
#
#   Correct marker model:          M | A, Y, X
#   Misspecified marker model:     M | X
#
# Here, X contains age, sex, and race/ethnicity. The DGM, generating settings,
# primary datasets, and both-correct AIPW estimator remain unchanged.
################################################################################

source("00_config.R")
source("02_dgm.R")
source("03_methods.R")

library(dplyr)

# Prespecified generating settings from the primary simulation.
aipw_strong_settings <- tibble::tribble(
  ~setting_id, ~setting, ~primary_scenario_id, ~N,
  ~phase2_fraction, ~r2_a_marker, ~r2_y_marker, ~theta,
  1L, "MIDUS-like", 168L, 800L, 0.65, 0.05, 0.05, 0.15,
  2L, "Stress", 208L, 500L, 0.25, 0.10, 0.10, 0.15
)

# Confirm that each sensitivity setting matches its primary scenario.
for (i in seq_len(nrow(aipw_strong_settings))) {
  setting <- aipw_strong_settings[i, , drop = FALSE]
  
  primary_scenario <- primary_grid[
    primary_grid$scenario_id == setting$primary_scenario_id,
    ,
    drop = FALSE
  ]
  
  if (nrow(primary_scenario) != 1L) {
    stop(
      "Primary scenario ID not uniquely found: ",
      setting$primary_scenario_id,
      "."
    )
  }
  
  setting_matches <- c(
    primary_scenario$N == setting$N,
    abs(
      primary_scenario$phase2_fraction -
        setting$phase2_fraction
    ) < 1e-12,
    abs(
      primary_scenario$r2_a_marker -
        setting$r2_a_marker
    ) < 1e-12,
    abs(
      primary_scenario$r2_y_marker -
        setting$r2_y_marker
    ) < 1e-12,
    abs(
      primary_scenario$theta -
        setting$theta
    ) < 1e-12
  )
  
  if (!all(setting_matches)) {
    stop(
      "Strong-misspecification setting does not match primary scenario ",
      setting$primary_scenario_id,
      "."
    )
  }
}

# Four combinations of correct and misspecified nuisance models.
aipw_strong_specs <- tibble::tribble(
  ~spec_id, ~specification, ~selection_spec, ~augmentation_spec,
  ~selection_correct, ~augmentation_correct,
  1L, "Both correct", "full", "full", TRUE, TRUE,
  2L, "Selection correct only", "full", "x_only", TRUE, FALSE,
  3L, "Augmentation correct only", "x_only", "full", FALSE, TRUE,
  4L, "Both misspecified", "x_only", "x_only", FALSE, FALSE
)

# Construct a nuisance-model design matrix.
#
# full:   A + Y + age + sex + race/ethnicity
# x_only: age + sex + race/ethnicity
strong_nuisance_matrix <- function(
    dat,
    spec = c("full", "x_only")
) {
  spec <- match.arg(spec)
  
  if (spec == "full") {
    model.matrix(
      ~ A + Y + age_STD + sex + race_eth,
      data = dat
    )
  } else {
    model.matrix(
      ~ age_STD + sex + race_eth,
      data = dat
    )
  }
}

# Retain separate function names for selection and marker models.
strong_selection_matrix <- function(
    dat,
    spec = c("full", "x_only")
) {
  strong_nuisance_matrix(
    dat = dat,
    spec = match.arg(spec)
  )
}

strong_marker_predictor_matrix <- function(
    dat,
    spec = c("full", "x_only")
) {
  strong_nuisance_matrix(
    dat = dat,
    spec = match.arg(spec)
  )
}

# Fit the sensitivity-specific AIPW estimator.
#
# The estimating equations and stacked sandwich variance are unchanged from
# fit_aipw(). Only the available nuisance-model specifications differ.
fit_aipw_strong_misspec <- function(
    dat,
    marker_names,
    selection_spec = c("full", "x_only"),
    augmentation_spec = c("full", "x_only"),
    method_label = "AIPW",
    prob_floor = probability_floor
) {
  selection_spec <- match.arg(selection_spec)
  augmentation_spec <- match.arg(augmentation_spec)
  
  phase2 <- dat$phase2
  n <- nrow(dat)
  observed <- which(phase2 == 1)
  n_observed <- length(observed)
  
  if (n_observed <= length(marker_names) + 10L) {
    stop("Too few Phase-2 observations for AIPW.")
  }
  
  # Fit the Phase-2 selection model.
  selection_design <- strong_selection_matrix(
    dat,
    spec = selection_spec
  )
  
  selection_fit <- glm.fit(
    x = selection_design,
    y = phase2,
    family = binomial()
  )
  
  if (!isTRUE(selection_fit$converged)) {
    stop("AIPW selection model did not converge.")
  }
  
  alpha <- selection_fit$coefficients
  
  if (any(!is.finite(alpha))) {
    stop("Non-finite coefficient in AIPW selection model.")
  }
  
  pi_hat <- as.vector(
    plogis(selection_design %*% alpha)
  )
  
  if (
    any(!is.finite(pi_hat)) ||
    any(pi_hat < prob_floor | pi_hat > 1 - prob_floor)
  ) {
    stop(
      "Estimated Phase-2 probabilities are too close to zero or one."
    )
  }
  
  inverse_probability_factor <- phase2 / pi_hat
  
  # Fit the multivariate marker model among Phase-2 participants.
  marker_design <- strong_marker_predictor_matrix(
    dat,
    spec = augmentation_spec
  )
  
  marker_design_observed <- marker_design[
    observed,
    ,
    drop = FALSE
  ]
  
  marker_observed <- as.matrix(
    dat[observed, marker_names, drop = FALSE]
  )
  
  marker_coefficients <- qr.solve(
    marker_design_observed,
    marker_observed
  )
  
  if (any(!is.finite(marker_coefficients))) {
    stop("Non-finite coefficient in the AIPW marker model.")
  }
  
  marker_mean <- marker_design %*% marker_coefficients
  colnames(marker_mean) <- marker_names
  
  marker_residuals <- (
    marker_observed -
      marker_design_observed %*% marker_coefficients
  )
  
  # The MLE covariance makes the covariance estimating equations sum to zero.
  marker_covariance <- (
    crossprod(marker_residuals) /
      n_observed
  )
  
  covariance_eigenvalues <- eigen(
    marker_covariance,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  if (
    any(!is.finite(covariance_eigenvalues)) ||
    min(covariance_eigenvalues) <= matrix_tolerance
  ) {
    stop("Estimated AIPW marker covariance is singular or nearly singular.")
  }
  
  # Construct the target outcome model.
  outcome_design <- target_w_matrix(dat)
  
  n_outcome_terms <- ncol(outcome_design)
  n_markers <- length(marker_names)
  n_target_terms <- n_outcome_terms + n_markers
  
  expected_complete_design <- cbind(
    outcome_design,
    marker_mean
  )
  
  complete_design <- expected_complete_design
  complete_design[
    observed,
    (n_outcome_terms + 1L):n_target_terms
  ] <- marker_observed
  
  marker_variance_block <- matrix(
    0,
    nrow = n_target_terms,
    ncol = n_target_terms
  )
  
  marker_variance_block[
    (n_outcome_terms + 1L):n_target_terms,
    (n_outcome_terms + 1L):n_target_terms
  ] <- marker_covariance
  
  sqrt_ip_factor <- sqrt(inverse_probability_factor)
  
  estimating_matrix <- (
    crossprod(expected_complete_design) +
      n * marker_variance_block +
      crossprod(complete_design * sqrt_ip_factor) -
      crossprod(expected_complete_design * sqrt_ip_factor) -
      sum(inverse_probability_factor) * marker_variance_block
  )
  
  estimating_vector <- (
    colSums(expected_complete_design * dat$Y) +
      colSums(
        (complete_design - expected_complete_design) *
          (inverse_probability_factor * dat$Y)
      )
  )
  
  target_coefficients <- solve(
    estimating_matrix,
    estimating_vector
  )
  
  target_coefficients <- as.vector(target_coefficients)
  names(target_coefficients) <- colnames(expected_complete_design)
  
  # Construct the unit-level target estimating equations.
  outcome_coefficients <- target_coefficients[
    seq_len(n_outcome_terms)
  ]
  
  marker_effects <- target_coefficients[
    (n_outcome_terms + 1L):n_target_terms
  ]
  
  expected_residual <- (
    dat$Y -
      as.vector(outcome_design %*% outcome_coefficients) -
      as.vector(marker_mean %*% marker_effects)
  )
  
  expected_top <- outcome_design * expected_residual
  
  covariance_marker_effect <- as.vector(
    marker_covariance %*% marker_effects
  )
  
  expected_bottom <- sweep(
    marker_mean * expected_residual,
    MARGIN = 2,
    STATS = covariance_marker_effect,
    FUN = "-"
  )
  
  expected_score <- cbind(
    expected_top,
    expected_bottom
  )
  
  complete_residual <- (
    dat$Y -
      as.vector(complete_design %*% target_coefficients)
  )
  
  complete_score <- complete_design * complete_residual
  
  target_score <- (
    expected_score +
      (complete_score - expected_score) *
      inverse_probability_factor
  )
  
  # Construct the nuisance-model estimating equations.
  selection_score <- (
    selection_design *
      (phase2 - pi_hat)
  )
  
  n_marker_predictors <- ncol(marker_design)
  
  residuals_all <- matrix(
    0,
    nrow = n,
    ncol = n_markers
  )
  
  residuals_all[observed, ] <- marker_residuals
  
  marker_coefficient_scores <- lapply(
    seq_len(n_markers),
    function(j) {
      marker_design * (
        phase2 * residuals_all[, j]
      )
    }
  )
  
  marker_coefficient_score <- do.call(
    cbind,
    marker_coefficient_scores
  )
  
  covariance_pairs <- vech_pairs(n_markers)
  n_covariance_terms <- length(covariance_pairs)
  
  covariance_score <- matrix(
    0,
    nrow = n,
    ncol = n_covariance_terms
  )
  
  for (j in seq_along(covariance_pairs)) {
    row_index <- covariance_pairs[[j]][1]
    column_index <- covariance_pairs[[j]][2]
    
    covariance_score[, j] <- phase2 * (
      residuals_all[, row_index] *
        residuals_all[, column_index] -
        marker_covariance[row_index, column_index]
    )
  }
  
  # Construct the stacked Jacobian.
  n_selection_terms <- ncol(selection_design)
  
  n_marker_coefficient_terms <- (
    n_marker_predictors * n_markers
  )
  
  jacobian_selection <- -crossprod(
    selection_design,
    selection_design * (
      pi_hat * (1 - pi_hat)
    )
  )
  
  jacobian_marker_coefficients <- -kronecker(
    diag(n_markers),
    crossprod(marker_design_observed)
  )
  
  jacobian_covariance <- (
    -n_observed * diag(n_covariance_terms)
  )
  
  derivative_ip_alpha <- selection_design * (
    -phase2 * (1 - pi_hat) / pi_hat
  )
  
  jacobian_target_selection <- crossprod(
    complete_score - expected_score,
    derivative_ip_alpha
  )
  
  jacobian_target_marker <- matrix(
    0,
    nrow = n_target_terms,
    ncol = n_marker_coefficient_terms
  )
  
  one_minus_ip_factor <- (
    1 - inverse_probability_factor
  )
  
  for (j in seq_len(n_markers)) {
    derivative_mean <- matrix(
      0,
      nrow = n,
      ncol = n_target_terms
    )
    
    derivative_mean[, seq_len(n_outcome_terms)] <- (
      -outcome_design * marker_effects[j]
    )
    
    derivative_bottom <- (
      -marker_mean * marker_effects[j]
    )
    
    derivative_bottom[, j] <- (
      derivative_bottom[, j] +
        expected_residual
    )
    
    derivative_mean[
      ,
      (n_outcome_terms + 1L):n_target_terms
    ] <- derivative_bottom
    
    marker_columns <- (
      ((j - 1L) * n_marker_predictors + 1L):
        (j * n_marker_predictors)
    )
    
    jacobian_target_marker[, marker_columns] <- crossprod(
      derivative_mean * one_minus_ip_factor,
      marker_design
    )
  }
  
  covariance_derivative <- matrix(
    0,
    nrow = n_target_terms,
    ncol = n_covariance_terms
  )
  
  for (j in seq_along(covariance_pairs)) {
    row_index <- covariance_pairs[[j]][1]
    column_index <- covariance_pairs[[j]][2]
    
    derivative_bottom <- rep(
      0,
      n_markers
    )
    
    if (row_index == column_index) {
      derivative_bottom[row_index] <- (
        -marker_effects[row_index]
      )
    } else {
      derivative_bottom[row_index] <- (
        -marker_effects[column_index]
      )
      
      derivative_bottom[column_index] <- (
        -marker_effects[row_index]
      )
    }
    
    covariance_derivative[
      (n_outcome_terms + 1L):n_target_terms,
      j
    ] <- derivative_bottom
  }
  
  jacobian_target_covariance <- (
    sum(one_minus_ip_factor) *
      covariance_derivative
  )
  
  jacobian_target <- -estimating_matrix
  
  n_parameters <- (
    n_selection_terms +
      n_marker_coefficient_terms +
      n_covariance_terms +
      n_target_terms
  )
  
  jacobian <- matrix(
    0,
    nrow = n_parameters,
    ncol = n_parameters
  )
  
  selection_indices <- seq_len(n_selection_terms)
  
  marker_indices <- (
    n_selection_terms +
      seq_len(n_marker_coefficient_terms)
  )
  
  covariance_indices <- (
    n_selection_terms +
      n_marker_coefficient_terms +
      seq_len(n_covariance_terms)
  )
  
  target_indices <- (
    n_selection_terms +
      n_marker_coefficient_terms +
      n_covariance_terms +
      seq_len(n_target_terms)
  )
  
  jacobian[
    selection_indices,
    selection_indices
  ] <- jacobian_selection
  
  jacobian[
    marker_indices,
    marker_indices
  ] <- jacobian_marker_coefficients
  
  jacobian[
    covariance_indices,
    covariance_indices
  ] <- jacobian_covariance
  
  # At the OLS solution, the covariance-score derivative with respect to the
  # marker coefficients sums to zero because Z'e = 0 for each marker.
  jacobian[
    target_indices,
    selection_indices
  ] <- jacobian_target_selection
  
  jacobian[
    target_indices,
    marker_indices
  ] <- jacobian_target_marker
  
  jacobian[
    target_indices,
    covariance_indices
  ] <- jacobian_target_covariance
  
  jacobian[
    target_indices,
    target_indices
  ] <- jacobian_target
  
  stacked_scores <- cbind(
    selection_score,
    marker_coefficient_score,
    covariance_score,
    target_score
  )
  
  sandwich_meat <- crossprod(stacked_scores)
  
  inverse_jacobian <- solve(jacobian)
  
  stacked_variance <- (
    inverse_jacobian %*%
      sandwich_meat %*%
      t(inverse_jacobian)
  )
  
  target_variance <- stacked_variance[
    target_indices,
    target_indices,
    drop = FALSE
  ]
  
  exposure_index <- match(
    "A",
    names(target_coefficients)
  )
  
  if (is.na(exposure_index)) {
    stop("A coefficient not found in the AIPW target model.")
  }
  
  estimate <- target_coefficients[exposure_index]
  
  estimate_variance <- target_variance[
    exposure_index,
    exposure_index
  ]
  
  if (
    !is.finite(estimate) ||
    !is.finite(estimate_variance) ||
    estimate_variance <= 0
  ) {
    stop("Invalid AIPW estimate or sandwich variance.")
  }
  
  standard_error <- sqrt(estimate_variance)
  statistic <- estimate / standard_error
  p_value <- 2 * pnorm(-abs(statistic))
  critical_value <- qnorm(1 - alpha_level / 2)
  
  observed_weights <- 1 / pi_hat[observed]
  
  data.frame(
    method = method_label,
    estimate = estimate,
    se = standard_error,
    df = Inf,
    p_value = p_value,
    conf_low = estimate - critical_value * standard_error,
    conf_high = estimate + critical_value * standard_error,
    status = "ok",
    message = NA_character_,
    weight_min = min(observed_weights),
    weight_p99 = unname(
      quantile(observed_weights, 0.99)
    ),
    weight_max = max(observed_weights),
    weight_cv = (
      sd(observed_weights) /
        mean(observed_weights)
    ),
    weight_ess = (
      sum(observed_weights)^2 /
        sum(observed_weights^2)
    )
  )
}

# Return the cached DGM filename for a generating setting.
strong_cache_filename <- function(
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

# Reconstruct the primary RNG stream for a scenario.
get_strong_primary_scenario_stream <- function(scenario_id) {
  scenario_id <- as.integer(scenario_id)
  
  if (
    length(scenario_id) != 1L ||
    is.na(scenario_id) ||
    scenario_id < 1L
  ) {
    stop("scenario_id must be a positive integer.")
  }
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(rng_seed_master)
  
  state <- .Random.seed
  
  if (scenario_id > 1L) {
    for (i in seq_len(scenario_id - 1L)) {
      state <- parallel::nextRNGStream(state)
    }
  }
  
  state
}

# Return the data and method RNG states for one primary repetition.
get_strong_primary_repetition_states <- function(
    scenario_id,
    repetition
) {
  repetition <- as.integer(repetition)
  
  if (
    length(repetition) != 1L ||
    is.na(repetition) ||
    repetition < 1L ||
    repetition > nsim_primary
  ) {
    stop(
      "repetition must be between 1 and ",
      nsim_primary,
      "."
    )
  }
  
  state <- get_strong_primary_scenario_stream(
    scenario_id
  )
  
  n_substreams <- 2L * (repetition - 1L)
  
  if (n_substreams > 0L) {
    for (i in seq_len(n_substreams)) {
      state <- parallel::nextRNGSubStream(state)
    }
  }
  
  data_state <- state
  method_state <- parallel::nextRNGSubStream(
    data_state
  )
  
  list(
    data_state = data_state,
    method_state = method_state
  )
}

# Store an RNG state in a CSV-compatible form.
strong_state_to_string <- function(state) {
  paste(state, collapse = ",")
}

# Fit a selection model used only for misspecification diagnostics.
fit_strong_selection_diagnostic <- function(
    dat,
    spec,
    prob_floor = probability_floor
) {
  design <- strong_selection_matrix(
    dat,
    spec = spec
  )
  
  fit <- glm.fit(
    x = design,
    y = dat$phase2,
    family = binomial()
  )
  
  if (!isTRUE(fit$converged)) {
    stop("Strong selection diagnostic did not converge.")
  }
  
  coefficients <- fit$coefficients
  
  if (any(!is.finite(coefficients))) {
    stop("Non-finite coefficient in strong selection diagnostic.")
  }
  
  pi_hat <- as.vector(
    plogis(design %*% coefficients)
  )
  
  if (
    any(!is.finite(pi_hat)) ||
    any(pi_hat < prob_floor | pi_hat > 1 - prob_floor)
  ) {
    stop(
      "Strong diagnostic selection probabilities are too close to zero or one."
    )
  }
  
  list(
    alpha = coefficients,
    pi_hat = pi_hat
  )
}

# Fit a marker model used only for misspecification diagnostics.
fit_strong_marker_diagnostic <- function(
    dat,
    marker_names,
    spec
) {
  observed <- which(dat$phase2 == 1)
  
  marker_design <- strong_marker_predictor_matrix(
    dat,
    spec = spec
  )
  
  marker_design_observed <- marker_design[
    observed,
    ,
    drop = FALSE
  ]
  
  if (length(observed) <= ncol(marker_design_observed)) {
    stop(
      "Too few Phase-2 observations for the marker diagnostic."
    )
  }
  
  marker_observed <- as.matrix(
    dat[observed, marker_names, drop = FALSE]
  )
  
  marker_coefficients <- qr.solve(
    marker_design_observed,
    marker_observed
  )
  
  if (any(!is.finite(marker_coefficients))) {
    stop("Non-finite coefficient in strong marker diagnostic.")
  }
  
  marker_mean <- marker_design %*% marker_coefficients
  colnames(marker_mean) <- marker_names
  
  marker_residuals <- (
    marker_observed -
      marker_design_observed %*% marker_coefficients
  )
  
  marker_covariance <- (
    crossprod(marker_residuals) /
      length(observed)
  )
  
  list(
    Bhat = marker_coefficients,
    mu = marker_mean,
    Sigma_hat = marker_covariance
  )
}

# Compare the correct and X-only nuisance-model fits.
compute_strong_misspec_diagnostics <- function(
    dat,
    marker_names
) {
  selection_full <- fit_strong_selection_diagnostic(
    dat,
    spec = "full"
  )
  
  selection_x_only <- fit_strong_selection_diagnostic(
    dat,
    spec = "x_only"
  )
  
  marker_full <- fit_strong_marker_diagnostic(
    dat,
    marker_names,
    spec = "full"
  )
  
  marker_x_only <- fit_strong_marker_diagnostic(
    dat,
    marker_names,
    spec = "x_only"
  )
  
  probability_difference <- (
    selection_full$pi_hat -
      selection_x_only$pi_hat
  )
  
  marker_mean_difference <- (
    marker_full$mu -
      marker_x_only$mu
  )
  
  covariance_difference <- (
    marker_full$Sigma_hat -
      marker_x_only$Sigma_hat
  )
  
  covariance_denominator <- sqrt(
    sum(marker_full$Sigma_hat^2)
  )
  
  relative_covariance_difference <- if (
    is.finite(covariance_denominator) &&
    covariance_denominator > 0
  ) {
    sqrt(
      sum(covariance_difference^2)
    ) / covariance_denominator
  } else {
    NA_real_
  }
  
  data.frame(
    mean_abs_pi_difference = mean(
      abs(probability_difference)
    ),
    max_abs_pi_difference = max(
      abs(probability_difference)
    ),
    rmse_pi_difference = sqrt(
      mean(probability_difference^2)
    ),
    cor_pi_correct_misspecified = suppressWarnings(
      cor(
        selection_full$pi_hat,
        selection_x_only$pi_hat
      )
    ),
    mean_abs_marker_mean_difference = mean(
      abs(marker_mean_difference)
    ),
    max_abs_marker_mean_difference = max(
      abs(marker_mean_difference)
    ),
    rmse_marker_mean_difference = sqrt(
      mean(marker_mean_difference^2)
    ),
    relative_frobenius_sigma_difference =
      relative_covariance_difference,
    pi_correct_min = min(selection_full$pi_hat),
    pi_correct_max = max(selection_full$pi_hat),
    pi_misspecified_min = min(
      selection_x_only$pi_hat
    ),
    pi_misspecified_max = max(
      selection_x_only$pi_hat
    )
  )
}

# Run one nuisance-model specification and convert errors into failed rows.
run_one_strong_aipw_spec <- function(
    dat,
    marker_names,
    spec_row
) {
  method_label <- spec_row$specification[[1]]
  start_time <- proc.time()[["elapsed"]]
  
  result <- tryCatch(
    fit_aipw_strong_misspec(
      dat = dat,
      marker_names = marker_names,
      selection_spec =
        spec_row$selection_spec[[1]],
      augmentation_spec =
        spec_row$augmentation_spec[[1]],
      method_label = method_label
    ),
    error = function(e) {
      failed_result(
        method_label,
        conditionMessage(e)
      )
    }
  )
  
  result$elapsed_seconds <- (
    proc.time()[["elapsed"]] -
      start_time
  )
  
  result$spec_id <- spec_row$spec_id[[1]]
  result$specification <- spec_row$specification[[1]]
  result$selection_spec <- spec_row$selection_spec[[1]]
  result$augmentation_spec <- spec_row$augmentation_spec[[1]]
  result$selection_correct <- spec_row$selection_correct[[1]]
  result$augmentation_correct <- spec_row$augmentation_correct[[1]]
  
  result
}

# Run all four specifications on the same generated dataset.
run_strong_aipw_specs <- function(
    dat,
    marker_names
) {
  results <- vector(
    "list",
    nrow(aipw_strong_specs)
  )
  
  for (i in seq_len(nrow(aipw_strong_specs))) {
    results[[i]] <- run_one_strong_aipw_spec(
      dat = dat,
      marker_names = marker_names,
      spec_row = aipw_strong_specs[
        i,
        ,
        drop = FALSE
      ]
    )
  }
  
  results <- bind_rows(results)
  
  if (nrow(results) != nrow(aipw_strong_specs)) {
    stop(
      "Strong AIPW sensitivity did not return exactly ",
      nrow(aipw_strong_specs),
      " rows."
    )
  }
  
  results
}

