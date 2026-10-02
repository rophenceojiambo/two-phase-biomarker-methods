# 10_aipw_robustness_helpers.R
#
# Defines the settings and helper functions for the AIPW double-robustness
# experiment.
#
# The primary data-generating mechanism and validated AIPW estimator are not
# modified. The experiment uses the two existing model switches in fit_aipw():
#
# - selection_include_y controls whether the selection model includes Y;
# - marker_include_y controls whether the marker model includes Y.
#
# Each generated dataset is analyzed under four working-model specifications:
#
# 1. Both correct
#    Selection: S | A, Y, X
#    Augmentation: M | A, Y, X
#
# 2. Selection correct only
#    Selection: S | A, Y, X
#    Augmentation: M | A, X
#
# 3. Augmentation correct only
#    Selection: S | A, X
#    Augmentation: M | A, Y, X
#
# 4. Both misspecified
#    Selection: S | A, X
#    Augmentation: M | A, X
#
# The same generated dataset is used for all four specifications.

# Load the primary simulation settings, DGM, and analysis methods
source("00_config.R")
source("02_dgm.R")
source("03_methods.R")

library(dplyr)

# 1. Prespecified generating settings

# These are primary scenarios whose RNG streams and cached DGM parameters
# already exist. Reusing them permits exact comparison with the primary AIPW
# estimates.
aipw_robustness_settings <- tibble::tribble(
  ~setting_id, ~setting, ~primary_scenario_id, ~N,
  ~phase2_fraction, ~r2_a_marker, ~r2_y_marker, ~theta,
  
  1L, "MIDUS-like", 168L, 800L,
  0.65, 0.05, 0.05, 0.15,
  
  2L, "Stress", 208L, 500L,
  0.25, 0.10, 0.10, 0.15
)

# Confirm that the stored scenario IDs still agree with the primary design
for (i in seq_len(nrow(aipw_robustness_settings))) {
  setting_row <- aipw_robustness_settings[
    i,
    ,
    drop = FALSE
  ]
  
  primary_row <- primary_grid[
    primary_grid$scenario_id ==
      setting_row$primary_scenario_id,
    ,
    drop = FALSE
  ]
  
  if (nrow(primary_row) != 1L) {
    stop(
      "Primary scenario ID not found: ",
      setting_row$primary_scenario_id
    )
  }
  
  setting_matches <- c(
    primary_row$N == setting_row$N,
    
    abs(
      primary_row$phase2_fraction -
        setting_row$phase2_fraction
    ) < 1e-12,
    
    abs(
      primary_row$r2_a_marker -
        setting_row$r2_a_marker
    ) < 1e-12,
    
    abs(
      primary_row$r2_y_marker -
        setting_row$r2_y_marker
    ) < 1e-12,
    
    abs(
      primary_row$theta -
        setting_row$theta
    ) < 1e-12
  )
  
  if (!all(setting_matches)) {
    stop(
      "Robustness setting does not agree with primary scenario ",
      setting_row$primary_scenario_id,
      "."
    )
  }
}

# 2. AIPW working-model specifications

aipw_robustness_specs <- tibble::tribble(
  ~spec_id, ~specification,
  ~selection_correct, ~augmentation_correct,
  
  1L, "Both correct",
  TRUE, TRUE,
  
  2L, "Selection correct only",
  TRUE, FALSE,
  
  3L, "Augmentation correct only",
  FALSE, TRUE,
  
  4L, "Both misspecified",
  FALSE, FALSE
)

# 3. DGM-cache filename

# Reproduce the cache naming rule used in 04_build_dgm_cache.R
robustness_cache_filename <- function(
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

# 4. Primary random-number streams

# Reproduce the scenario-level stream from 05_run_primary_chunk.R
get_primary_scenario_stream <- function(scenario_id) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(rng_seed_master)
  
  state <- .Random.seed
  
  if (scenario_id > 1L) {
    for (stream_index in seq_len(scenario_id - 1L)) {
      state <- parallel::nextRNGStream(state)
    }
  }
  
  state
}

# Return the data and method substreams for one primary repetition
get_primary_repetition_states <- function(
    scenario_id,
    repetition
) {
  state <- get_primary_scenario_stream(scenario_id)
  n_advance <- 2L * (repetition - 1L)
  
  if (n_advance > 0L) {
    for (substream_index in seq_len(n_advance)) {
      state <- parallel::nextRNGSubStream(state)
    }
  }
  
  data_state <- state
  method_state <- parallel::nextRNGSubStream(data_state)
  
  list(
    data_state = data_state,
    method_state = method_state
  )
}

# Convert an RNG state to text for storage in the canary output
state_to_string <- function(state) {
  paste(state, collapse = ",")
}

# 5. Nuisance-model diagnostics

# These models quantify how much omitting Y changes the nuisance estimates.
# They do not replace the models fitted internally by fit_aipw().

fit_selection_diagnostic <- function(
    dat,
    include_y,
    prob_floor = probability_floor
) {
  design_matrix <- selection_matrix(
    dat,
    include_y = include_y
  )
  
  selection_fit <- glm.fit(
    x = design_matrix,
    y = dat$phase2,
    family = binomial()
  )
  
  coefficients <- selection_fit$coefficients
  
  if (any(!is.finite(coefficients))) {
    stop("Non-finite coefficient in diagnostic selection model.")
  }
  
  probabilities <- as.vector(
    plogis(design_matrix %*% coefficients)
  )
  
  if (
    any(
      probabilities < prob_floor |
      probabilities > 1 - prob_floor
    )
  ) {
    stop("Diagnostic selection probabilities too close to 0/1.")
  }
  
  list(
    alpha = coefficients,
    pi_hat = probabilities
  )
}

fit_marker_diagnostic <- function(
    dat,
    marker_names,
    include_y
) {
  observed <- which(dat$phase2 == 1)
  
  marker_design <- marker_predictor_matrix(
    dat,
    include_y = include_y
  )
  
  marker_design_observed <- marker_design[
    observed,
    ,
    drop = FALSE
  ]
  
  markers_observed <- as.matrix(
    dat[
      observed,
      marker_names,
      drop = FALSE
    ]
  )
  
  if (length(observed) <= ncol(marker_design_observed)) {
    stop("Too few phase-2 observations for diagnostic marker model.")
  }
  
  coefficients <- qr.solve(
    marker_design_observed,
    markers_observed
  )
  
  marker_means <- marker_design %*% coefficients
  colnames(marker_means) <- marker_names
  
  marker_residuals <- markers_observed -
    marker_design_observed %*% coefficients
  
  residual_covariance <- crossprod(
    marker_residuals
  ) / length(observed)
  
  list(
    Bhat = coefficients,
    mu = marker_means,
    Sigma_hat = residual_covariance
  )
}

# Compare the correct and misspecified nuisance-model estimates
compute_misspecification_diagnostics <- function(
    dat,
    marker_names
) {
  selection_correct <- fit_selection_diagnostic(
    dat,
    include_y = TRUE
  )
  
  selection_misspecified <- fit_selection_diagnostic(
    dat,
    include_y = FALSE
  )
  
  marker_correct <- fit_marker_diagnostic(
    dat,
    marker_names,
    include_y = TRUE
  )
  
  marker_misspecified <- fit_marker_diagnostic(
    dat,
    marker_names,
    include_y = FALSE
  )
  
  probability_difference <-
    selection_correct$pi_hat -
    selection_misspecified$pi_hat
  
  marker_mean_difference <-
    marker_correct$mu -
    marker_misspecified$mu
  
  covariance_difference <-
    marker_correct$Sigma_hat -
    marker_misspecified$Sigma_hat
  
  covariance_denominator <- sqrt(
    sum(marker_correct$Sigma_hat^2)
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
  
  probability_correlation <- suppressWarnings(
    cor(
      selection_correct$pi_hat,
      selection_misspecified$pi_hat
    )
  )
  
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
    
    cor_pi_correct_misspecified =
      probability_correlation,
    
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
    
    pi_correct_min = min(
      selection_correct$pi_hat
    ),
    
    pi_correct_max = max(
      selection_correct$pi_hat
    ),
    
    pi_misspecified_min = min(
      selection_misspecified$pi_hat
    ),
    
    pi_misspecified_max = max(
      selection_misspecified$pi_hat
    )
  )
}

# 6. Fit the four AIPW specifications

# Fit one working-model specification and retain a standardized failure row
run_one_aipw_spec <- function(
    dat,
    marker_names,
    spec_row
) {
  method_label <- spec_row$specification[[1]]
  start_time <- proc.time()[["elapsed"]]
  
  result <- tryCatch(
    fit_aipw(
      dat = dat,
      marker_names = marker_names,
      selection_include_y =
        spec_row$selection_correct[[1]],
      marker_include_y =
        spec_row$augmentation_correct[[1]],
      method_label = method_label
    ),
    error = function(e) {
      failed_result(
        method_label,
        conditionMessage(e)
      )
    }
  )
  
  end_time <- proc.time()[["elapsed"]]
  
  result$elapsed_seconds <- end_time - start_time
  result$spec_id <- spec_row$spec_id[[1]]
  result$specification <- spec_row$specification[[1]]
  result$selection_correct <- spec_row$selection_correct[[1]]
  result$augmentation_correct <-
    spec_row$augmentation_correct[[1]]
  
  result
}

# Apply all four working-model specifications to the same dataset
run_aipw_robustness_specs <- function(
    dat,
    marker_names
) {
  results <- vector(
    "list",
    nrow(aipw_robustness_specs)
  )
  
  for (spec_index in seq_len(nrow(aipw_robustness_specs))) {
    results[[spec_index]] <- run_one_aipw_spec(
      dat = dat,
      marker_names = marker_names,
      spec_row = aipw_robustness_specs[
        spec_index,
        ,
        drop = FALSE
      ]
    )
  }
  
  combined_results <- bind_rows(results)
  
  if (
    nrow(combined_results) !=
    nrow(aipw_robustness_specs)
  ) {
    stop(
      "AIPW robustness analysis did not return exactly four rows."
    )
  }
  
  combined_results
}
