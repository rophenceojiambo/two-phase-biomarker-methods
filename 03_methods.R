# 03_methods.R
#
# Implements the six analysis methods compared in the MIDUS two-phase
# simulation:
# 1. Naive phase-1 analysis that omits the phase-2 RNA markers
# 2. Complete-case analysis restricted to the phase-2 sample
# 3. Fully conditional specification multiple imputation using mice
# 4. Joint-model multiple imputation using jomo
# 5. Inverse probability weighting with stacked sandwich variance
# 6. Augmented inverse probability weighting with stacked sandwich variance
#
# The AIPW estimator models phase-2 selection and the conditional distribution
# of the RNA markers. Its stacked sandwich variance includes uncertainty from
# the selection model, marker regression, marker residual covariance, and
# target outcome regression.
#
# Each fitting function returns the exposure-effect estimate, standard error,
# confidence interval, p-value, and method status. Sourcing this script defines
# the functions but does not fit any models.

# Load shared settings
source("00_config.R")

library(dplyr)
library(mice)
library(jomo)

# 1. Common helper functions

# Construct the outcome model containing all RNA markers
full_formula <- function(marker_names) {
  as.formula(
    paste(
      "Y ~ A + age_STD + sex + race_eth +",
      paste(marker_names, collapse = " + ")
    )
  )
}

# Construct the phase-1 outcome model that omits the RNA markers
naive_formula <- function() {
  Y ~ A + age_STD + sex + race_eth
}

# Construct the design matrix for the phase-2 selection model
selection_matrix <- function(dat, include_y = TRUE) {
  if (include_y) {
    model.matrix(~ A + Y + age_STD + sex + race_eth, data = dat)
  } else {
    model.matrix(~ A + age_STD + sex + race_eth, data = dat)
  }
}

# Construct the predictor matrix for the conditional marker model
marker_predictor_matrix <- function(dat, include_y = TRUE) {
  if (include_y) {
    model.matrix(~ A + Y + age_STD + sex + race_eth, data = dat)
  } else {
    model.matrix(~ A + age_STD + sex + race_eth, data = dat)
  }
}

# Construct the phase-1 component of the target outcome design matrix
target_w_matrix <- function(dat) {
  model.matrix(~ A + age_STD + sex + race_eth, data = dat)
}

# Extract the exposure coefficient and its inference from a linear model
extract_lm_A <- function(fit, method_name) {
  coefficients <- summary(fit)$coefficients
  
  if (!("A" %in% rownames(coefficients))) {
    stop("A coefficient not found.")
  }
  
  estimate <- coefficients["A", "Estimate"]
  se <- coefficients["A", "Std. Error"]
  df <- df.residual(fit)
  statistic <- estimate / se
  p_value <- 2 * pt(-abs(statistic), df = df)
  critical_value <- qt(0.975, df = df)
  
  data.frame(
    method = method_name,
    estimate = estimate,
    se = se,
    df = df,
    p_value = p_value,
    conf_low = estimate - critical_value * se,
    conf_high = estimate + critical_value * se,
    status = "ok",
    message = NA_character_
  )
}

# Return a standardized result when a method cannot be fitted
failed_result <- function(method_name, msg) {
  data.frame(
    method = method_name,
    estimate = NA_real_,
    se = NA_real_,
    df = NA_real_,
    p_value = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    status = "failed",
    message = as.character(msg)
  )
}

# Allow one method to fail without stopping the remaining methods
safe_run <- function(method_name, expr) {
  tryCatch(
    expr,
    error = function(e) {
      failed_result(method_name, conditionMessage(e))
    }
  )
}

# Pool exposure estimates and variances using Rubin's rules
pool_A_from_lm_list <- function(fits, N_complete) {
  estimates <- vapply(
    fits,
    function(fit) coef(fit)["A"],
    numeric(1)
  )
  
  variances <- vapply(
    fits,
    function(fit) vcov(fit)["A", "A"],
    numeric(1)
  )
  
  n_parameters <- length(coef(fits[[1]]))
  
  pooled <- mice::pool.scalar(
    Q = estimates,
    U = variances,
    n = N_complete,
    k = n_parameters,
    rule = "rubin1987"
  )
  
  estimate <- pooled$qbar
  se <- sqrt(pooled$t)
  df <- pooled$df
  statistic <- estimate / se
  p_value <- 2 * pt(-abs(statistic), df = df)
  critical_value <- qt(0.975, df = df)
  
  c(
    estimate = estimate,
    se = se,
    df = df,
    p_value = p_value,
    conf_low = estimate - critical_value * se,
    conf_high = estimate + critical_value * se
  )
}

# 2. Naive phase-1 analysis

# Fit Y | A, X to the full phase-1 sample without adjusting for RNA markers
fit_naive <- function(dat) {
  fit <- lm(naive_formula(), data = dat)
  extract_lm_A(fit, "Naive")
}

# 3. Complete-case analysis

# Fit the complete outcome model using only participants selected into phase 2
fit_cca <- function(dat, marker_names) {
  complete_cases <- dat[dat$phase2 == 1, , drop = FALSE]
  
  if (nrow(complete_cases) <= length(marker_names) + 10L) {
    stop("Too few Phase-2 observations for complete-case model.")
  }
  
  fit <- lm(
    full_formula(marker_names),
    data = complete_cases
  )
  
  extract_lm_A(fit, "CCA")
}

# 4. Fully conditional specification multiple imputation

# Impute the continuous RNA markers using chained Bayesian linear regressions
fit_fcs_mi <- function(
    dat,
    marker_names,
    nimp = nimp_primary,
    maxit = mice_maxit_primary,
    method = mice_method_primary
) {
  imp_data <- dat %>%
    select(
      A, Y, age_STD, sex, race_eth,
      all_of(marker_names)
    )
  
  # Obtain the default method vector and predictor matrix
  initialization <- mice::mice(
    imp_data,
    maxit = 0,
    printFlag = FALSE
  )
  
  methods <- initialization$method
  predictors <- initialization$predictorMatrix
  
  # Impute only the RNA markers
  methods[] <- ""
  methods[marker_names] <- method
  
  # Use A, Y, X, and the other markers to impute each RNA marker
  predictors[,] <- 0
  candidate_predictors <- c(
    "A", "Y", "age_STD", "sex", "race_eth",
    marker_names
  )
  
  for (marker in marker_names) {
    predictors[
      marker,
      setdiff(candidate_predictors, marker)
    ] <- 1
  }
  
  imputed <- mice::mice(
    imp_data,
    m = nimp,
    maxit = maxit,
    method = methods,
    predictorMatrix = predictors,
    printFlag = FALSE
  )
  
  # Fit the complete outcome model in each imputed dataset
  fits <- lapply(seq_len(nimp), function(j) {
    completed <- mice::complete(imputed, action = j)
    lm(full_formula(marker_names), data = completed)
  })
  
  pooled <- pool_A_from_lm_list(
    fits = fits,
    N_complete = nrow(dat)
  )
  
  log_message <- if (is.null(imputed$loggedEvents)) {
    NA_character_
  } else {
    paste("mice logged events:", nrow(imputed$loggedEvents))
  }
  
  data.frame(
    method = "FCS-MI",
    estimate = pooled["estimate"],
    se = pooled["se"],
    df = pooled["df"],
    p_value = pooled["p_value"],
    conf_low = pooled["conf_low"],
    conf_high = pooled["conf_high"],
    status = "ok",
    message = log_message
  )
}

# 5. Joint-model multiple imputation

# Jointly impute the RNA markers using a multivariate normal model
fit_joint_mi <- function(
    dat,
    marker_names,
    nimp = nimp_primary,
    nburn = jomo_nburn_primary,
    nbetween = jomo_nbetween_primary
) {
  # RNA markers are the jointly imputed outcomes
  Yimp <- as.data.frame(dat[, marker_names, drop = FALSE])
  
  # A, Y, and X are fully observed predictors in the imputation model
  Ximp <- model.matrix(
    ~ A + Y + age_STD + sex + race_eth,
    data = dat
  )
  
  imputed_long <- jomo::jomo(
    Y = Yimp,
    X = Ximp,
    nburn = as.integer(nburn),
    nbetween = as.integer(nbetween),
    nimp = as.integer(nimp),
    output = 0
  )
  
  if (!("Imputation" %in% names(imputed_long))) {
    stop("jomo output does not contain the Imputation index.")
  }
  
  fits <- vector("list", nimp)
  
  # Reconstruct and analyze each imputed dataset
  for (j in seq_len(nimp)) {
    imputed_markers <- imputed_long[
      imputed_long$Imputation == j,
      marker_names,
      drop = FALSE
    ]
    
    if (nrow(imputed_markers) != nrow(dat)) {
      stop(
        "Unexpected jomo output size for imputation ", j,
        ": expected ", nrow(dat),
        ", got ", nrow(imputed_markers), "."
      )
    }
    
    completed <- dat %>%
      select(A, Y, age_STD, sex, race_eth)
    
    completed[marker_names] <- imputed_markers
    
    fits[[j]] <- lm(
      full_formula(marker_names),
      data = completed
    )
  }
  
  pooled <- pool_A_from_lm_list(
    fits = fits,
    N_complete = nrow(dat)
  )
  
  data.frame(
    method = "JM-MI",
    estimate = pooled["estimate"],
    se = pooled["se"],
    df = pooled["df"],
    p_value = pooled["p_value"],
    conf_low = pooled["conf_low"],
    conf_high = pooled["conf_high"],
    status = "ok",
    message = NA_character_
  )
}

# 6. Inverse probability weighting

# Fit the full outcome model in phase 2 using inverse selection probabilities
# and estimate variance with stacked selection and outcome equations
fit_ipw <- function(
    dat,
    marker_names,
    prob_floor = probability_floor
) {
  S <- dat$phase2
  N <- nrow(dat)
  
  # Fit the phase-2 selection model S | A, Y, X
  Z <- selection_matrix(dat, include_y = TRUE)
  
  selection_fit <- glm.fit(
    x = Z,
    y = S,
    family = binomial()
  )
  
  alpha <- selection_fit$coefficients
  
  if (any(!is.finite(alpha))) {
    stop("Non-finite coefficient in IPW selection model.")
  }
  
  pi_hat <- as.vector(plogis(Z %*% alpha))
  
  # Stop rather than silently truncate unstable probabilities
  if (any(pi_hat < prob_floor | pi_hat > 1 - prob_floor)) {
    stop("Estimated Phase-2 probabilities too close to 0/1 for stable IPW.")
  }
  
  observed <- which(S == 1)
  
  # Construct the complete outcome design for phase-2 participants
  C_obs <- model.matrix(
    full_formula(marker_names),
    data = dat[observed, , drop = FALSE]
  )
  
  Y_obs <- dat$Y[observed]
  pi_obs <- pi_hat[observed]
  weights <- 1 / pi_obs
  
  # Solve the weighted least-squares estimating equation
  XtWX <- crossprod(C_obs, C_obs * weights)
  XtWY <- crossprod(C_obs, Y_obs * weights)
  
  beta <- as.vector(solve(XtWX, XtWY))
  names(beta) <- colnames(C_obs)
  
  # Stack the selection and weighted outcome estimating equations
  q <- ncol(Z)
  d <- ncol(C_obs)
  
  psi_alpha <- Z * (S - pi_hat)
  
  psi_beta <- matrix(0, nrow = N, ncol = d)
  colnames(psi_beta) <- colnames(C_obs)
  
  residuals_obs <- Y_obs - as.vector(C_obs %*% beta)
  outcome_scores_obs <- C_obs * residuals_obs
  
  psi_beta[observed, ] <- outcome_scores_obs / pi_obs
  
  # Construct the Jacobian blocks
  J_aa <- -crossprod(
    Z,
    Z * (pi_hat * (1 - pi_hat))
  )
  
  J_bb <- -crossprod(
    C_obs,
    C_obs * (1 / pi_obs)
  )
  
  J_ba <- -crossprod(
    outcome_scores_obs,
    Z[observed, , drop = FALSE] *
      ((1 - pi_obs) / pi_obs)
  )
  
  J <- rbind(
    cbind(J_aa, matrix(0, q, d)),
    cbind(J_ba, J_bb)
  )
  
  # Calculate the stacked sandwich covariance matrix
  scores <- cbind(psi_alpha, psi_beta)
  meat <- crossprod(scores)
  
  J_inverse <- solve(J)
  stacked_variance <- J_inverse %*% meat %*% t(J_inverse)
  
  beta_indices <- q + seq_len(d)
  beta_variance <- stacked_variance[
    beta_indices,
    beta_indices,
    drop = FALSE
  ]
  
  A_index <- match("A", names(beta))
  
  if (is.na(A_index)) {
    stop("A coefficient not found in IPW outcome design.")
  }
  
  estimate <- beta[A_index]
  se <- sqrt(beta_variance[A_index, A_index])
  statistic <- estimate / se
  p_value <- 2 * pnorm(-abs(statistic))
  critical_value <- qnorm(0.975)
  
  data.frame(
    method = "IPW",
    estimate = estimate,
    se = se,
    df = Inf,
    p_value = p_value,
    conf_low = estimate - critical_value * se,
    conf_high = estimate + critical_value * se,
    status = "ok",
    message = NA_character_,
    weight_min = min(weights),
    weight_p99 = unname(quantile(weights, 0.99)),
    weight_max = max(weights),
    weight_cv = sd(weights) / mean(weights),
    weight_ess = sum(weights)^2 / sum(weights^2)
  )
}

# 7. AIPW helper

# Return lower-triangular index pairs for the unique covariance parameters
vech_pairs <- function(k) {
  pairs <- vector("list", k * (k + 1) / 2)
  index <- 1L
  
  for (column in seq_len(k)) {
    for (row in column:k) {
      pairs[[index]] <- c(row, column)
      index <- index + 1L
    }
  }
  
  pairs
}

# 8. Augmented inverse probability weighting

# Fit the full outcome model using the observed complete-data score plus its
# conditional expectation under the marker model
fit_aipw <- function(
    dat,
    marker_names,
    selection_include_y = TRUE,
    marker_include_y = TRUE,
    method_label = "AIPW",
    prob_floor = probability_floor
) {
  S <- dat$phase2
  N <- nrow(dat)
  observed <- which(S == 1)
  
  if (length(observed) <= length(marker_names) + 10L) {
    stop("Too few Phase-2 observations for AIPW.")
  }
  
  # A. Fit the phase-2 selection model
  
  Zs <- selection_matrix(
    dat,
    include_y = selection_include_y
  )
  
  selection_fit <- glm.fit(
    x = Zs,
    y = S,
    family = binomial()
  )
  
  alpha <- selection_fit$coefficients
  
  if (any(!is.finite(alpha))) {
    stop("Non-finite coefficient in AIPW selection model.")
  }
  
  pi_hat <- as.vector(plogis(Zs %*% alpha))
  
  if (any(pi_hat < prob_floor | pi_hat > 1 - prob_floor)) {
    stop("Estimated Phase-2 probabilities too close to 0/1 for stable AIPW.")
  }
  
  inverse_probability_factor <- S / pi_hat
  
  # B. Fit the joint marker model M | A, Y, X
  
  Zm <- marker_predictor_matrix(
    dat,
    include_y = marker_include_y
  )
  
  M_obs <- as.matrix(
    dat[observed, marker_names, drop = FALSE]
  )
  
  Zm_obs <- Zm[observed, , drop = FALSE]
  
  # Estimate the conditional marker means
  Bhat <- qr.solve(Zm_obs, M_obs)
  marker_means <- Zm %*% Bhat
  colnames(marker_means) <- marker_names
  
  marker_residuals <- M_obs - Zm_obs %*% Bhat
  
  # Use the MLE covariance so the covariance score sums exactly to zero
  Sigma_hat <- crossprod(marker_residuals) / length(observed)
  
  minimum_eigenvalue <- min(
    eigen(
      Sigma_hat,
      symmetric = TRUE,
      only.values = TRUE
    )$values
  )
  
  if (minimum_eigenvalue <= matrix_tolerance) {
    stop("Estimated AIPW marker covariance is singular/near-singular.")
  }
  
  # C. Solve the AIPW target estimating equation
  
  W <- target_w_matrix(dat)
  r <- ncol(W)
  k <- length(marker_names)
  d <- r + k
  
  # Cbar contains W and the conditional marker means for everyone
  Cbar <- cbind(W, marker_means)
  
  # Replace conditional means with observed markers for phase-2 participants
  Cobs <- Cbar
  Cobs[observed, (r + 1):d] <- M_obs
  
  # Add the marker covariance to the lower-right design block
  Vsig <- matrix(0, nrow = d, ncol = d)
  Vsig[(r + 1):d, (r + 1):d] <- Sigma_hat
  
  sqrt_factor <- sqrt(inverse_probability_factor)
  
  A_sum <- crossprod(Cbar) +
    N * Vsig +
    crossprod(Cobs * sqrt_factor) -
    crossprod(Cbar * sqrt_factor) -
    sum(inverse_probability_factor) * Vsig
  
  b_sum <- colSums(Cbar * dat$Y) +
    colSums(
      (Cobs - Cbar) *
        (inverse_probability_factor * dat$Y)
    )
  
  beta <- as.vector(solve(A_sum, b_sum))
  names(beta) <- colnames(Cbar)
  
  # D. Construct the target AIPW estimating equations
  
  beta_w <- beta[seq_len(r)]
  beta_m <- beta[(r + 1):d]
  
  residual_bar <- dat$Y -
    as.vector(W %*% beta_w) -
    as.vector(marker_means %*% beta_m)
  
  expected_score_top <- W * residual_bar
  sigma_beta_m <- as.vector(Sigma_hat %*% beta_m)
  
  expected_score_bottom <- sweep(
    marker_means * residual_bar,
    MARGIN = 2,
    STATS = sigma_beta_m,
    FUN = "-"
  )
  
  expected_scores <- cbind(
    expected_score_top,
    expected_score_bottom
  )
  
  complete_residual <- dat$Y - as.vector(Cobs %*% beta)
  complete_scores <- Cobs * complete_residual
  
  psi_beta <- expected_scores +
    (complete_scores - expected_scores) *
    inverse_probability_factor
  
  # E. Construct the nuisance-model estimating equations
  
  # Selection-model score
  psi_alpha <- Zs * (S - pi_hat)
  
  qm <- ncol(Zm)
  
  # Marker-regression scores
  residuals_all <- matrix(0, nrow = N, ncol = k)
  residuals_all[observed, ] <- marker_residuals
  
  psi_B_blocks <- lapply(seq_len(k), function(j) {
    Zm * (S * residuals_all[, j])
  })
  
  psi_B <- do.call(cbind, psi_B_blocks)
  
  # Unique marker-covariance scores
  pairs <- vech_pairs(k)
  n_covariance_parameters <- length(pairs)
  
  psi_Sigma <- matrix(
    0,
    nrow = N,
    ncol = n_covariance_parameters
  )
  
  for (j in seq_along(pairs)) {
    row <- pairs[[j]][1]
    column <- pairs[[j]][2]
    
    psi_Sigma[, j] <- S * (
      residuals_all[, row] * residuals_all[, column] -
        Sigma_hat[row, column]
    )
  }
  
  # F. Construct the stacked Jacobian
  
  qs <- ncol(Zs)
  n_marker_coefficients <- qm * k
  
  # Selection, marker-regression, and marker-covariance blocks
  J_aa <- -crossprod(
    Zs,
    Zs * (pi_hat * (1 - pi_hat))
  )
  
  J_BB <- -kronecker(
    diag(k),
    crossprod(Zm_obs)
  )
  
  J_SS <- -length(observed) *
    diag(n_covariance_parameters)
  
  # Derivative of the target score with respect to selection coefficients
  dh_dalpha <- Zs * (
    -S * (1 - pi_hat) / pi_hat
  )
  
  J_beta_alpha <- crossprod(
    complete_scores - expected_scores,
    dh_dalpha
  )
  
  # Derivative of the target score with respect to marker coefficients
  J_beta_B <- matrix(
    0,
    nrow = d,
    ncol = n_marker_coefficients
  )
  
  one_minus_factor <- 1 - inverse_probability_factor
  
  for (j in seq_len(k)) {
    derivative_column <- matrix(0, nrow = N, ncol = d)
    
    # Derivative of the upper score block with respect to marker mean j
    derivative_column[, seq_len(r)] <- -W * beta_m[j]
    
    # Derivative of the lower score block with respect to marker mean j
    bottom <- -marker_means * beta_m[j]
    bottom[, j] <- bottom[, j] + residual_bar
    
    derivative_column[, (r + 1):d] <- bottom
    
    coefficient_indices <- ((j - 1) * qm + 1):(j * qm)
    
    J_beta_B[, coefficient_indices] <- crossprod(
      derivative_column * one_minus_factor,
      Zm
    )
  }
  
  # Derivative of the target score with respect to covariance parameters
  covariance_derivative <- matrix(
    0,
    nrow = d,
    ncol = n_covariance_parameters
  )
  
  for (j in seq_along(pairs)) {
    row <- pairs[[j]][1]
    column <- pairs[[j]][2]
    derivative_bottom <- rep(0, k)
    
    if (row == column) {
      derivative_bottom[row] <- -beta_m[row]
    } else {
      derivative_bottom[row] <- -beta_m[column]
      derivative_bottom[column] <- -beta_m[row]
    }
    
    covariance_derivative[
      (r + 1):d,
      j
    ] <- derivative_bottom
  }
  
  J_beta_Sigma <- sum(one_minus_factor) *
    covariance_derivative
  
  J_beta_beta <- -A_sum
  
  # Assemble the full block Jacobian
  total_parameters <- qs +
    n_marker_coefficients +
    n_covariance_parameters +
    d
  
  J <- matrix(
    0,
    nrow = total_parameters,
    ncol = total_parameters
  )
  
  i_alpha <- seq_len(qs)
  i_B <- qs + seq_len(n_marker_coefficients)
  i_Sigma <- qs +
    n_marker_coefficients +
    seq_len(n_covariance_parameters)
  i_beta <- qs +
    n_marker_coefficients +
    n_covariance_parameters +
    seq_len(d)
  
  J[i_alpha, i_alpha] <- J_aa
  J[i_B, i_B] <- J_BB
  J[i_Sigma, i_Sigma] <- J_SS
  
  # The covariance-score derivative with respect to B is zero at the OLS root
  J[i_beta, i_alpha] <- J_beta_alpha
  J[i_beta, i_B] <- J_beta_B
  J[i_beta, i_Sigma] <- J_beta_Sigma
  J[i_beta, i_beta] <- J_beta_beta
  
  # Calculate the stacked sandwich covariance matrix
  scores <- cbind(
    psi_alpha,
    psi_B,
    psi_Sigma,
    psi_beta
  )
  
  meat <- crossprod(scores)
  J_inverse <- solve(J)
  stacked_variance <- J_inverse %*% meat %*% t(J_inverse)
  
  beta_variance <- stacked_variance[
    i_beta,
    i_beta,
    drop = FALSE
  ]
  
  A_index <- match("A", names(beta))
  
  if (is.na(A_index)) {
    stop("A coefficient not found in AIPW target design.")
  }
  
  estimate <- beta[A_index]
  se <- sqrt(beta_variance[A_index, A_index])
  statistic <- estimate / se
  p_value <- 2 * pnorm(-abs(statistic))
  critical_value <- qnorm(0.975)
  
  # Report inverse-probability diagnostics among phase-2 participants
  observed_weights <- 1 / pi_hat[observed]
  
  data.frame(
    method = method_label,
    estimate = estimate,
    se = se,
    df = Inf,
    p_value = p_value,
    conf_low = estimate - critical_value * se,
    conf_high = estimate + critical_value * se,
    status = "ok",
    message = NA_character_,
    weight_min = min(observed_weights),
    weight_p99 = unname(quantile(observed_weights, 0.99)),
    weight_max = max(observed_weights),
    weight_cv = sd(observed_weights) / mean(observed_weights),
    weight_ess = sum(observed_weights)^2 / sum(observed_weights^2)
  )
}

# 9. Run all six primary methods

# Fit every method and retain a failure record when an individual method errors
run_all_methods <- function(
    dat,
    marker_names,
    nimp = nimp_primary,
    mice_maxit = mice_maxit_primary,
    jomo_nburn = jomo_nburn_primary,
    jomo_nbetween = jomo_nbetween_primary
) {
  results <- list(
    safe_run(
      "Naive",
      fit_naive(dat)
    ),
    safe_run(
      "CCA",
      fit_cca(dat, marker_names)
    ),
    safe_run(
      "FCS-MI",
      fit_fcs_mi(
        dat,
        marker_names,
        nimp = nimp,
        maxit = mice_maxit,
        method = mice_method_primary
      )
    ),
    safe_run(
      "JM-MI",
      fit_joint_mi(
        dat,
        marker_names,
        nimp = nimp,
        nburn = jomo_nburn,
        nbetween = jomo_nbetween
      )
    ),
    safe_run(
      "IPW",
      fit_ipw(dat, marker_names)
    ),
    safe_run(
      "AIPW",
      fit_aipw(
        dat,
        marker_names,
        selection_include_y = TRUE,
        marker_include_y = TRUE,
        method_label = "AIPW"
      )
    )
  )
  
  bind_rows(results)
}
