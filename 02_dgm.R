# 02_dgm.R
#
# Defines the data-generating mechanism for the MIDUS two-phase simulation.
# The functions in this script:
# - resample phase-1 covariates from the empirical MIDUS distribution;
# - generate correlated RNA markers conditional on those covariates;
# - generate the exposure and outcome with target marker partial R-squared values;
# - select participants into phase 2 under MAR or MCAR; and
# - set the RNA markers to missing outside the phase-2 sample.
#
# Exposure and outcome parameters are calibrated to produce approximately
# mean-zero, unit-variance variables while preserving the specified exposure
# effect and marker partial R-squared values.
#
# Sourcing this script defines the functions but does not generate data.

# Load shared settings
source("00_config.R")

library(mvtnorm)
library(dplyr)

# 1. Helper functions

# Construct an X matrix with the factor levels and columns used in calibration
make_x_matrix <- function(x_data, calibration) {
  x_data <- x_data %>%
    mutate(
      sex = factor(sex, levels = calibration$sex_levels),
      race_eth = factor(race_eth, levels = calibration$race_levels)
    )
  
  X <- model.matrix(calibration$x_formula, data = x_data)
  
  # Add zero-filled columns when a factor level is absent from the current sample
  missing_cols <- setdiff(calibration$x_columns, colnames(X))
  
  if (length(missing_cols) > 0L) {
    X <- cbind(
      X,
      matrix(
        0,
        nrow = nrow(X),
        ncol = length(missing_cols),
        dimnames = list(NULL, missing_cols)
      )
    )
  }
  
  # Return columns in the same order used during calibration
  X[, calibration$x_columns, drop = FALSE]
}

# Return the smallest positive solution to a quadratic equation
solve_positive_quadratic <- function(a, b, c, label = "quadratic") {
  if (!is.finite(a) || !is.finite(b) || !is.finite(c) || a <= 0) {
    stop("Invalid coefficients in ", label, ".")
  }
  
  discriminant <- b^2 - 4 * a * c
  
  if (discriminant < -1e-10) {
    stop(
      label,
      " has no real solution. Discriminant = ",
      signif(discriminant, 6)
    )
  }
  
  # Treat a very small negative value as numerical rounding
  discriminant <- max(discriminant, 0)
  
  roots <- c(
    (-b + sqrt(discriminant)) / (2 * a),
    (-b - sqrt(discriminant)) / (2 * a)
  )
  
  positive_roots <- roots[is.finite(roots) & roots > 0]
  
  if (length(positive_roots) == 0L) {
    stop(label, " has no positive solution.")
  }
  
  min(positive_roots)
}

# 2. Calibrate the exposure model

# Choose exposure parameters that produce the target partial R-squared for
# M in the model A | X, M while setting the marginal variance of A to one
calibrate_exposure_parameters <- function(calibration, r2_a_marker) {
  if (r2_a_marker <= 0 || r2_a_marker >= 1) {
    stop("r2_a_marker must lie strictly between 0 and 1.")
  }
  
  X <- make_x_matrix(calibration$x_empirical, calibration)
  X_no <- X[, calibration$x_no_intercept, drop = FALSE]
  mu_M <- X %*% calibration$B_marker
  
  gamma_X <- calibration$gamma_X_empirical
  dA <- calibration$gamma_M_direction
  Sigma <- calibration$Sigma_marker
  
  # Linear predictors contributed by X and the conditional marker means
  x_lp <- as.vector(X_no %*% gamma_X)
  z_mu <- as.vector(mu_M %*% dA)
  
  # Residual marker variance in the empirical marker direction
  v_u <- as.numeric(t(dA) %*% Sigma %*% dA)
  
  var_z <- var(z_mu) + v_u
  cov_x_z <- cov(x_lp, z_mu)
  var_x <- var(x_lp)
  
  # Solve for the marker scale that gives Var(A) = 1
  quad_a <- var_z + v_u * (1 - r2_a_marker) / r2_a_marker
  quad_b <- 2 * cov_x_z
  quad_c <- var_x - 1
  
  marker_scale <- solve_positive_quadratic(
    quad_a,
    quad_b,
    quad_c,
    label = "Exposure calibration"
  )
  
  gamma_M <- marker_scale * dA
  
  # Set residual variance to achieve the target marker partial R-squared
  sigma_A2 <- marker_scale^2 * v_u *
    (1 - r2_a_marker) / r2_a_marker
  
  # Center the exposure at zero
  gamma0 <- -mean(x_lp + as.vector(mu_M %*% gamma_M))
  
  a_mean_x <- gamma0 + x_lp + as.vector(mu_M %*% gamma_M)
  
  # Check the resulting variance and partial R-squared
  marker_signal_A <- as.numeric(
    t(gamma_M) %*% Sigma %*% gamma_M
  )
  
  var_A_check <- var(a_mean_x) + marker_signal_A + sigma_A2
  
  partial_r2_check <- marker_signal_A /
    (marker_signal_A + sigma_A2)
  
  list(
    gamma0 = gamma0,
    gamma_X = gamma_X,
    gamma_M = gamma_M,
    sigma_A2 = sigma_A2,
    var_A_check = var_A_check,
    partial_r2_check = partial_r2_check
  )
}

# 3. Calibrate the outcome model

# Choose outcome parameters that produce the target partial R-squared for
# M in the model Y | A, X, M while setting the marginal variance of Y to one
calibrate_outcome_parameters <- function(
    calibration,
    exposure_parameters,
    r2_y_marker,
    theta
) {
  if (r2_y_marker < 0 || r2_y_marker >= 1) {
    stop("r2_y_marker must lie in [0, 1).")
  }
  
  X <- make_x_matrix(calibration$x_empirical, calibration)
  X_no <- X[, calibration$x_no_intercept, drop = FALSE]
  mu_M <- X %*% calibration$B_marker
  Sigma <- calibration$Sigma_marker
  
  gamma_X <- exposure_parameters$gamma_X
  gamma_M <- exposure_parameters$gamma_M
  sigma_A2 <- exposure_parameters$sigma_A2
  gamma0 <- exposure_parameters$gamma0
  
  beta_X <- calibration$beta_X_empirical
  dY <- calibration$beta_M_direction
  
  # Conditional means contributed by the exposure and phase-1 covariates
  x_gamma <- as.vector(X_no %*% gamma_X)
  a_mean_x <- gamma0 + x_gamma + as.vector(mu_M %*% gamma_M)
  x_beta <- as.vector(X_no %*% beta_X)
  zy_mu <- as.vector(mu_M %*% dY)
  
  # Derive the residual marker covariance conditional on A and X
  var_A_given_X <- as.numeric(
    t(gamma_M) %*% Sigma %*% gamma_M
  ) + sigma_A2
  
  cov_MA_given_X <- Sigma %*% gamma_M
  
  Sigma_M_given_AX <- Sigma -
    cov_MA_given_X %*% t(cov_MA_given_X) / var_A_given_X
  
  v_res_marker_direction <- as.numeric(
    t(dY) %*% Sigma_M_given_AX %*% dY
  )
  
  # Variance contributed by theta*A and X*beta
  var_A <- exposure_parameters$var_A_check
  var_x_beta <- var(x_beta)
  cov_A_xbeta <- cov(a_mean_x, x_beta)
  
  var_base <- theta^2 * var_A +
    var_x_beta +
    2 * theta * cov_A_xbeta
  
  # Total marker variance and its covariance with the base outcome predictor
  v_zy_u <- as.numeric(t(dY) %*% Sigma %*% dY)
  var_zy <- var(zy_mu) + v_zy_u
  
  cov_A_zy <- cov(a_mean_x, zy_mu) +
    as.numeric(t(gamma_M) %*% Sigma %*% dY)
  
  cov_xbeta_zy <- cov(x_beta, zy_mu)
  cov_base_zy <- theta * cov_A_zy + cov_xbeta_zy
  
  if (r2_y_marker == 0) {
    # Remove the marker contribution when its target partial R-squared is zero
    beta_M <- setNames(rep(0, length(dY)), names(dY))
    sigma_Y2 <- 1 - var_base
    
    if (sigma_Y2 <= 0) {
      stop(
        "Outcome calibration with r2_y_marker = 0 produced sigma_Y2 <= 0."
      )
    }
  } else {
    # Solve for the marker scale that gives Var(Y) = 1
    quad_a <- var_zy +
      v_res_marker_direction * (1 - r2_y_marker) / r2_y_marker
    quad_b <- 2 * cov_base_zy
    quad_c <- var_base - 1
    
    marker_scale <- solve_positive_quadratic(
      quad_a,
      quad_b,
      quad_c,
      label = "Outcome calibration"
    )
    
    beta_M <- marker_scale * dY
    
    # Set residual variance to achieve the target marker partial R-squared
    sigma_Y2 <- marker_scale^2 * v_res_marker_direction *
      (1 - r2_y_marker) / r2_y_marker
  }
  
  # Center the outcome at zero
  beta0 <- -mean(
    theta * a_mean_x +
      x_beta +
      as.vector(mu_M %*% beta_M)
  )
  
  # Check the marker partial R-squared conditional on A and X
  marker_signal_resid <- as.numeric(
    t(beta_M) %*% Sigma_M_given_AX %*% beta_M
  )
  
  partial_r2_check <- if (r2_y_marker == 0) {
    0
  } else {
    marker_signal_resid / (marker_signal_resid + sigma_Y2)
  }
  
  # Calculate all components needed to check the marginal variance of Y
  cov_A_Mbeta <- cov(a_mean_x, as.vector(mu_M %*% beta_M)) +
    as.numeric(t(gamma_M) %*% Sigma %*% beta_M)
  
  cov_xbeta_Mbeta <- cov(
    x_beta,
    as.vector(mu_M %*% beta_M)
  )
  
  var_Mbeta <- var(as.vector(mu_M %*% beta_M)) +
    as.numeric(t(beta_M) %*% Sigma %*% beta_M)
  
  var_Y_check <- theta^2 * var_A +
    var_x_beta +
    var_Mbeta +
    2 * theta * cov_A_xbeta +
    2 * theta * cov_A_Mbeta +
    2 * cov_xbeta_Mbeta +
    sigma_Y2
  
  list(
    beta0 = beta0,
    beta_X = beta_X,
    beta_M = beta_M,
    sigma_Y2 = sigma_Y2,
    theta = theta,
    var_Y_check = var_Y_check,
    partial_r2_check = partial_r2_check,
    Sigma_M_given_AX = Sigma_M_given_AX
  )
}

# 4. Simulate complete data

# Generate X, M, A, and Y before introducing phase-2 marker missingness
simulate_complete_data <- function(
    N,
    calibration,
    exposure_parameters,
    outcome_parameters,
    marker_error = c("mvn", "empirical")
) {
  marker_error <- match.arg(marker_error)
  
  # Resample phase-1 covariate profiles from the analytic MIDUS sample
  draw <- sample.int(
    nrow(calibration$x_empirical),
    size = N,
    replace = TRUE
  )
  
  xdat <- calibration$x_empirical[draw, , drop = FALSE] %>%
    mutate(
      sex = factor(sex, levels = calibration$sex_levels),
      race_eth = factor(race_eth, levels = calibration$race_levels)
    )
  
  X <- make_x_matrix(xdat, calibration)
  X_no <- X[, calibration$x_no_intercept, drop = FALSE]
  mu_M <- X %*% calibration$B_marker
  
  # Generate marker residuals from a multivariate normal distribution or
  # by resampling the empirical residual vectors
  if (marker_error == "mvn") {
    U <- mvtnorm::rmvnorm(
      N,
      mean = rep(0, length(calibration$marker_names)),
      sigma = calibration$Sigma_marker
    )
  } else {
    residual_draw <- sample.int(
      nrow(calibration$marker_residuals),
      size = N,
      replace = TRUE
    )
    
    U <- calibration$marker_residuals[
      residual_draw,
      ,
      drop = FALSE
    ]
  }
  
  M <- mu_M + U
  colnames(M) <- calibration$marker_names
  
  # Generate the continuous exposure
  A <- exposure_parameters$gamma0 +
    as.vector(X_no %*% exposure_parameters$gamma_X) +
    as.vector(M %*% exposure_parameters$gamma_M) +
    rnorm(
      N,
      mean = 0,
      sd = sqrt(exposure_parameters$sigma_A2)
    )
  
  # Generate the continuous outcome with true exposure effect theta
  Y <- outcome_parameters$beta0 +
    outcome_parameters$theta * A +
    as.vector(X_no %*% outcome_parameters$beta_X) +
    as.vector(M %*% outcome_parameters$beta_M) +
    rnorm(
      N,
      mean = 0,
      sd = sqrt(outcome_parameters$sigma_Y2)
    )
  
  data.frame(
    id = seq_len(N),
    A = A,
    Y = Y,
    age_STD = xdat$age_STD,
    sex = xdat$sex,
    race_eth = xdat$race_eth,
    M,
    check.names = FALSE
  )
}

# 5. Calibrate the phase-2 selection intercept

# Find the logistic-model intercept that produces the target phase-2 fraction
calibrate_selection_intercept <- function(
    calibration,
    exposure_parameters,
    outcome_parameters,
    target_fraction,
    n_cal = 50000L,
    seed = 1L,
    marker_error = "mvn"
) {
  set.seed(seed)
  
  # Use a large reference sample to approximate the expected selection fraction
  ref <- simulate_complete_data(
    N = n_cal,
    calibration = calibration,
    exposure_parameters = exposure_parameters,
    outcome_parameters = outcome_parameters,
    marker_error = marker_error
  )
  
  X <- make_x_matrix(
    ref[, c("age_STD", "sex", "race_eth")],
    calibration
  )
  
  X_no <- X[, calibration$x_no_intercept, drop = FALSE]
  
  # Construct the selection predictor without its intercept
  lp_no_intercept <- calibration$eta_A * ref$A +
    calibration$eta_Y * ref$Y +
    as.vector(X_no %*% calibration$eta_X)
  
  # Solve for an intercept whose mean selection probability equals the target
  objective <- function(eta0) {
    mean(plogis(eta0 + lp_no_intercept)) - target_fraction
  }
  
  eta0 <- uniroot(
    objective,
    interval = c(-30, 30),
    tol = 1e-12
  )$root
  
  list(
    eta0 = eta0,
    expected_fraction_check = mean(plogis(eta0 + lp_no_intercept))
  )
}

# 6. Build the parameters for one DGM

# Calibrate all fixed parameters once for a specified simulation scenario
build_dgm_parameters <- function(
    calibration,
    r2_a_marker,
    r2_y_marker,
    theta,
    phase2_fraction,
    selection = c("mar", "mcar"),
    selection_calibration_n = 50000L,
    selection_seed = 1L,
    marker_error = "mvn"
) {
  selection <- match.arg(selection)
  
  exposure_parameters <- calibrate_exposure_parameters(
    calibration,
    r2_a_marker = r2_a_marker
  )
  
  outcome_parameters <- calibrate_outcome_parameters(
    calibration,
    exposure_parameters = exposure_parameters,
    r2_y_marker = r2_y_marker,
    theta = theta
  )
  
  if (selection == "mar") {
    selection_parameters <- calibrate_selection_intercept(
      calibration = calibration,
      exposure_parameters = exposure_parameters,
      outcome_parameters = outcome_parameters,
      target_fraction = phase2_fraction,
      n_cal = selection_calibration_n,
      seed = selection_seed,
      marker_error = marker_error
    )
  } else {
    # MCAR uses a constant selection probability and requires no intercept
    selection_parameters <- list(
      eta0 = NA_real_,
      expected_fraction_check = phase2_fraction
    )
  }
  
  list(
    r2_a_marker = r2_a_marker,
    r2_y_marker = r2_y_marker,
    theta = theta,
    phase2_fraction = phase2_fraction,
    selection = selection,
    marker_error = marker_error,
    exposure = exposure_parameters,
    outcome = outcome_parameters,
    selection_parameters = selection_parameters
  )
}

# 7. Generate one observed two-phase dataset

# Generate complete data, sample phase 2, and remove unobserved RNA markers
generate_two_phase_data <- function(N, calibration, dgm_parameters) {
  full <- simulate_complete_data(
    N = N,
    calibration = calibration,
    exposure_parameters = dgm_parameters$exposure,
    outcome_parameters = dgm_parameters$outcome,
    marker_error = dgm_parameters$marker_error
  )
  
  X <- make_x_matrix(
    full[, c("age_STD", "sex", "race_eth")],
    calibration
  )
  
  X_no <- X[, calibration$x_no_intercept, drop = FALSE]
  
  if (dgm_parameters$selection == "mar") {
    # Under MAR, selection depends on fully observed A, Y, and X
    lp <- dgm_parameters$selection_parameters$eta0 +
      calibration$eta_A * full$A +
      calibration$eta_Y * full$Y +
      as.vector(X_no %*% calibration$eta_X)
    
    pi_true <- plogis(lp)
  } else {
    # Under MCAR, everyone has the same phase-2 selection probability
    pi_true <- rep(dgm_parameters$phase2_fraction, N)
  }
  
  # Draw the phase-2 indicator from the true selection probabilities
  S <- rbinom(N, size = 1, prob = pi_true)
  
  observed <- full
  observed$phase2 <- S
  observed$pi_true <- pi_true
  
  # RNA markers are observed only for participants selected into phase 2
  observed[S == 0, calibration$marker_names] <- NA_real_
  
  # Retain the complete markers in full for simulation validation
  full$phase2 <- S
  full$pi_true <- pi_true
  
  list(
    full = full,
    observed = observed
  )
}

# 8. Validate the DGM in a large simulated sample

# Compare the realized partial R-squared values, exposure effect, moments,
# phase-2 fraction, and selection-probability distribution with their targets
check_dgm <- function(
    calibration,
    dgm_parameters,
    N_check = 100000L,
    seed = 12345L
) {
  set.seed(seed)
  
  dat <- generate_two_phase_data(
    N = N_check,
    calibration = calibration,
    dgm_parameters = dgm_parameters
  )$full
  
  marker_rhs <- paste(calibration$marker_names, collapse = " + ")
  
  # Estimate the marker partial R-squared in the exposure model
  fit_A_reduced <- lm(
    A ~ age_STD + sex + race_eth,
    data = dat
  )
  
  fit_A_full <- lm(
    as.formula(
      paste("A ~ age_STD + sex + race_eth +", marker_rhs)
    ),
    data = dat
  )
  
  r2_A_reduced <- summary(fit_A_reduced)$r.squared
  r2_A_full <- summary(fit_A_full)$r.squared
  
  partial_r2_A <- (r2_A_full - r2_A_reduced) /
    (1 - r2_A_reduced)
  
  # Estimate the marker partial R-squared in the outcome model
  fit_Y_reduced <- lm(
    Y ~ A + age_STD + sex + race_eth,
    data = dat
  )
  
  fit_Y_full <- lm(
    as.formula(
      paste("Y ~ A + age_STD + sex + race_eth +", marker_rhs)
    ),
    data = dat
  )
  
  r2_Y_reduced <- summary(fit_Y_reduced)$r.squared
  r2_Y_full <- summary(fit_Y_full)$r.squared
  
  partial_r2_Y <- (r2_Y_full - r2_Y_reduced) /
    (1 - r2_Y_reduced)
  
  # Recover theta from the correctly specified full-data outcome model
  theta_hat_full <- coef(fit_Y_full)["A"]
  
  data.frame(
    target_r2_A_M_given_X = dgm_parameters$r2_a_marker,
    observed_r2_A_M_given_X = partial_r2_A,
    target_r2_Y_M_given_AX = dgm_parameters$r2_y_marker,
    observed_r2_Y_M_given_AX = partial_r2_Y,
    target_theta = dgm_parameters$theta,
    observed_theta_full_model = theta_hat_full,
    mean_A = mean(dat$A),
    sd_A = sd(dat$A),
    mean_Y = mean(dat$Y),
    sd_Y = sd(dat$Y),
    target_phase2_fraction = dgm_parameters$phase2_fraction,
    realized_phase2_fraction = mean(dat$phase2),
    min_pi = min(dat$pi_true),
    p01_pi = unname(quantile(dat$pi_true, 0.01)),
    median_pi = median(dat$pi_true),
    p99_pi = unname(quantile(dat$pi_true, 0.99)),
    max_pi = max(dat$pi_true)
  )
}
