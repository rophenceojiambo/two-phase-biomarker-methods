################################################################################
# 22_empirical_residual_sensitivity_helpers.R
#
# Defines scenarios and calibration tools for the empirical marker-residual
# sensitivity analysis.
#
# Instead of generating multivariate-normal marker errors, this sensitivity
# resamples empirical MIDUS residual vectors. The residual pool is transformed
# to retain the marker covariance used in the primary DGM, so the sensitivity
# changes residual distributional shape without changing the target covariance.
#
# This file does not modify the primary DGM or analysis methods.
################################################################################

source("00_config.R")

# Define null and alternative scenarios for two generating settings.
empirical_sensitivity_grid <- data.frame(
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

# Find the primary scenario with the same generating parameters.
match_primary_scenario_id <- function(
    scenario,
    primary_grid,
    tolerance = 1e-12
) {
  if (!is.data.frame(scenario) || nrow(scenario) != 1L) {
    stop("scenario must contain exactly one row.")
  }
  
  required_scenario_columns <- c(
    "scenario_label",
    "N",
    "phase2_fraction",
    "r2_a_marker",
    "r2_y_marker",
    "theta"
  )
  
  required_primary_columns <- c(
    "scenario_id",
    "N",
    "phase2_fraction",
    "r2_a_marker",
    "r2_y_marker",
    "theta"
  )
  
  missing_scenario_columns <- setdiff(
    required_scenario_columns,
    names(scenario)
  )
  
  missing_primary_columns <- setdiff(
    required_primary_columns,
    names(primary_grid)
  )
  
  if (length(missing_scenario_columns) > 0L) {
    stop(
      "Sensitivity scenario is missing: ",
      paste(missing_scenario_columns, collapse = ", "),
      "."
    )
  }
  
  if (length(missing_primary_columns) > 0L) {
    stop(
      "Primary grid is missing: ",
      paste(missing_primary_columns, collapse = ", "),
      "."
    )
  }
  
  matching_row <- (
    primary_grid$N == scenario$N &
      abs(
        primary_grid$phase2_fraction -
          scenario$phase2_fraction
      ) < tolerance &
      abs(
        primary_grid$r2_a_marker -
          scenario$r2_a_marker
      ) < tolerance &
      abs(
        primary_grid$r2_y_marker -
          scenario$r2_y_marker
      ) < tolerance &
      abs(
        primary_grid$theta -
          scenario$theta
      ) < tolerance
  )
  
  matching_ids <- primary_grid$scenario_id[
    matching_row
  ]
  
  if (length(matching_ids) != 1L) {
    stop(
      "Could not identify exactly one primary scenario for ",
      scenario$scenario_label,
      "."
    )
  }
  
  as.integer(matching_ids[[1]])
}

# Attach the corresponding primary scenario ID to each sensitivity scenario.
empirical_sensitivity_grid$primary_scenario_id <- vapply(
  seq_len(nrow(empirical_sensitivity_grid)),
  function(i) {
    match_primary_scenario_id(
      scenario = empirical_sensitivity_grid[
        i,
        ,
        drop = FALSE
      ],
      primary_grid = primary_grid
    )
  },
  integer(1)
)

empirical_sensitivity_grid <- empirical_sensitivity_grid[
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

# Protect against accidental changes to the primary scenario mapping.
expected_primary_scenario_ids <- c(
  60L,
  168L,
  100L,
  208L
)

if (
  !identical(
    empirical_sensitivity_grid$primary_scenario_id,
    expected_primary_scenario_ids
  )
) {
  stop(
    "The empirical sensitivity scenarios no longer match the expected ",
    "primary scenario IDs."
  )
}

# The canary uses the two alternative scenarios to test all six methods.
empirical_canary_grid <- empirical_sensitivity_grid[
  empirical_sensitivity_grid$theta == 0.15,
  ,
  drop = FALSE
]

if (nrow(empirical_canary_grid) != 2L) {
  stop("The empirical-residual canary must contain exactly two scenarios.")
}

# Center and covariance-match the empirical marker-residual vectors.
#
# Uniform resampling from a finite residual pool uses crossprod(U) / n as its
# population covariance. The primary calibration uses crossprod(U) / residual
# degrees of freedom. The transformation below maps the resampling covariance
# exactly to Sigma_marker while retaining the empirical residual shape.
prepare_empirical_residual_calibration <- function(
    calibration,
    tolerance = 1e-10
) {
  required_components <- c(
    "marker_names",
    "marker_residuals",
    "Sigma_marker"
  )
  
  missing_components <- setdiff(
    required_components,
    names(calibration)
  )
  
  if (length(missing_components) > 0L) {
    stop(
      "Calibration object is missing: ",
      paste(missing_components, collapse = ", "),
      "."
    )
  }
  
  marker_names <- calibration$marker_names
  residuals_raw <- as.matrix(
    calibration$marker_residuals
  )
  
  target_covariance <- as.matrix(
    calibration$Sigma_marker
  )
  
  n_markers <- length(marker_names)
  
  if (
    length(marker_names) == 0L ||
    anyNA(marker_names) ||
    anyDuplicated(marker_names)
  ) {
    stop("marker_names must contain unique, nonmissing names.")
  }
  
  if (ncol(residuals_raw) != n_markers) {
    stop(
      "marker_residuals must contain one column per marker."
    )
  }
  
  if (
    !identical(
      dim(target_covariance),
      c(n_markers, n_markers)
    )
  ) {
    stop("Sigma_marker has the wrong dimensions.")
  }
  
  if (nrow(residuals_raw) <= n_markers) {
    stop(
      "Too few empirical residual vectors for resampling."
    )
  }
  
  if (
    anyNA(residuals_raw) ||
    any(!is.finite(residuals_raw))
  ) {
    stop(
      "marker_residuals contains missing or non-finite values."
    )
  }
  
  if (
    anyNA(target_covariance) ||
    any(!is.finite(target_covariance))
  ) {
    stop(
      "Sigma_marker contains missing or non-finite values."
    )
  }
  
  # Align residual columns with marker_names when names are available.
  if (!is.null(colnames(residuals_raw))) {
    if (anyDuplicated(colnames(residuals_raw))) {
      stop(
        "marker_residuals contains duplicated column names."
      )
    }
    
    missing_residual_markers <- setdiff(
      marker_names,
      colnames(residuals_raw)
    )
    
    if (length(missing_residual_markers) > 0L) {
      stop(
        "marker_residuals is missing: ",
        paste(missing_residual_markers, collapse = ", "),
        "."
      )
    }
    
    residuals_raw <- residuals_raw[
      ,
      marker_names,
      drop = FALSE
    ]
  }
  
  # Align covariance rows and columns when complete names are available.
  if (
    !is.null(rownames(target_covariance)) &&
    !is.null(colnames(target_covariance))
  ) {
    if (
      anyDuplicated(rownames(target_covariance)) ||
      anyDuplicated(colnames(target_covariance))
    ) {
      stop(
        "Sigma_marker contains duplicated row or column names."
      )
    }
    
    missing_covariance_rows <- setdiff(
      marker_names,
      rownames(target_covariance)
    )
    
    missing_covariance_columns <- setdiff(
      marker_names,
      colnames(target_covariance)
    )
    
    if (
      length(missing_covariance_rows) > 0L ||
      length(missing_covariance_columns) > 0L
    ) {
      stop(
        "Sigma_marker names do not contain all marker_names."
      )
    }
    
    target_covariance <- target_covariance[
      marker_names,
      marker_names,
      drop = FALSE
    ]
  }
  
  colnames(residuals_raw) <- marker_names
  rownames(residuals_raw) <- NULL
  
  if (
    max(
      abs(
        target_covariance -
        t(target_covariance)
      )
    ) > tolerance
  ) {
    stop("Sigma_marker is not symmetric within tolerance.")
  }
  
  # Remove small numerical asymmetry before eigendecomposition and Cholesky.
  target_covariance <- (
    target_covariance +
      t(target_covariance)
  ) / 2
  
  residuals_centered <- sweep(
    residuals_raw,
    MARGIN = 2,
    STATS = colMeans(residuals_raw),
    FUN = "-"
  )
  
  n_residual_vectors <- nrow(
    residuals_centered
  )
  
  pool_covariance_before <- (
    crossprod(residuals_centered) /
      n_residual_vectors
  )
  
  pool_covariance_before <- (
    pool_covariance_before +
      t(pool_covariance_before)
  ) / 2
  
  pool_eigenvalues <- eigen(
    pool_covariance_before,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  target_eigenvalues <- eigen(
    target_covariance,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  if (
    any(!is.finite(pool_eigenvalues)) ||
    min(pool_eigenvalues) <= tolerance
  ) {
    stop(
      "The residual-pool covariance is not positive definite. ",
      "Minimum eigenvalue: ",
      signif(min(pool_eigenvalues), 6),
      "."
    )
  }
  
  if (
    any(!is.finite(target_eigenvalues)) ||
    min(target_eigenvalues) <= tolerance
  ) {
    stop(
      "Sigma_marker is not positive definite. ",
      "Minimum eigenvalue: ",
      signif(min(target_eigenvalues), 6),
      "."
    )
  }
  
  # If R'R = Sigma, then R_pool^{-1}R_target maps the pool covariance
  # to the target covariance.
  pool_cholesky <- chol(pool_covariance_before)
  target_cholesky <- chol(target_covariance)
  
  covariance_map <- solve(
    pool_cholesky,
    target_cholesky
  )
  
  residuals_matched <- (
    residuals_centered %*%
      covariance_map
  )
  
  colnames(residuals_matched) <- marker_names
  
  pool_covariance_after <- (
    crossprod(residuals_matched) /
      n_residual_vectors
  )
  
  matched_means <- colMeans(residuals_matched)
  
  covariance_error <- (
    pool_covariance_after -
      target_covariance
  )
  
  max_abs_mean_after <- max(
    abs(matched_means)
  )
  
  max_abs_covariance_error <- max(
    abs(covariance_error)
  )
  
  covariance_scale <- max(
    1,
    max(abs(target_covariance))
  )
  
  max_relative_covariance_error <- (
    max_abs_covariance_error /
      covariance_scale
  )
  
  if (max_abs_mean_after > 1e-10) {
    stop(
      "Matched residuals are not sufficiently centered. ",
      "Maximum absolute mean: ",
      signif(max_abs_mean_after, 6),
      "."
    )
  }
  
  if (max_relative_covariance_error > 1e-8) {
    stop(
      "Matched residuals do not reproduce Sigma_marker. ",
      "Maximum relative covariance error: ",
      signif(max_relative_covariance_error, 6),
      "."
    )
  }
  
  # Replace only the residual pool; all other calibration values are retained.
  empirical_calibration <- calibration
  empirical_calibration$marker_residuals <-
    residuals_matched
  
  pool_summary <- data.frame(
    n_residual_vectors = n_residual_vectors,
    n_markers = n_markers,
    max_abs_raw_mean = max(
      abs(colMeans(residuals_raw))
    ),
    max_abs_matched_mean = max_abs_mean_after,
    max_abs_covariance_error_before = max(
      abs(
        pool_covariance_before -
          target_covariance
      )
    ),
    max_abs_covariance_error_after =
      max_abs_covariance_error,
    max_relative_covariance_error_after =
      max_relative_covariance_error,
    min_eigenvalue_pool_before =
      min(pool_eigenvalues),
    min_eigenvalue_target =
      min(target_eigenvalues),
    stringsAsFactors = FALSE
  )
  
  # Describe the marginal shape of each matched residual distribution.
  marginal_diagnostics <- do.call(
    rbind,
    lapply(
      seq_len(n_markers),
      function(j) {
        residual_values <- residuals_matched[, j]
        
        centered_values <- (
          residual_values -
            mean(residual_values)
        )
        
        population_sd <- sqrt(
          mean(centered_values^2)
        )
        
        if (
          !is.finite(population_sd) ||
          population_sd <= 0
        ) {
          stop(
            "Zero or invalid empirical residual variance for marker ",
            marker_names[[j]],
            "."
          )
        }
        
        quantiles <- unname(
          quantile(
            residual_values,
            probs = c(
              0.01,
              0.05,
              0.50,
              0.95,
              0.99
            ),
            names = FALSE
          )
        )
        
        data.frame(
          marker = marker_names[[j]],
          mean = mean(residual_values),
          sd = sd(residual_values),
          skewness = mean(
            (centered_values / population_sd)^3
          ),
          excess_kurtosis = (
            mean(
              (centered_values / population_sd)^4
            ) - 3
          ),
          q01 = quantiles[[1]],
          q05 = quantiles[[2]],
          median = quantiles[[3]],
          q95 = quantiles[[4]],
          q99 = quantiles[[5]],
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  list(
    calibration = empirical_calibration,
    pool_summary = pool_summary,
    marginal_diagnostics = marginal_diagnostics,
    covariance_map = covariance_map
  )
}
