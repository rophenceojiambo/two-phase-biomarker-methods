################################################################################
# verify_sandwich_bootstrap.R
#
#
# Purpose:
#   Check the IPW and AIPW sandwich variance calculations in the stress
#   scenario.
#
# This script compares, on the SAME simulated datasets:
#
#   1. Current hand-derived analytic stacked-sandwich SE
#   2. Numerical stacked-sandwich SE
#      - same estimating equations
#      - Jacobian obtained by central finite differences
#   3. Nonparametric bootstrap SE
#      - resamples Phase-1 participants
#      - refits selection / marker / target models in every bootstrap sample
#
# Interpretation:
#
#   A. Analytic SE ~= Numerical SE
#      -> the hand-derived Jacobian is implemented consistently with the
#         stacked estimating equations.
#
#   B. Both sandwich SEs < Bootstrap SE, and bootstrap is near the empirical
#      SE from the 2,000-repetition validation
#      -> likely a genuine finite-sample limitation of the asymptotic sandwich.
#
#   C. Analytic SE differs materially from Numerical SE
#      -> inspect/fix the analytic Jacobian before the production run.
#
# Default diagnostic:
#   Stress scenario only
#   20 independently simulated datasets
#   250 bootstrap resamples per dataset
#
# This does NOT run FCS-MI or JM-MI.
################################################################################

source("00_config.R")
source("02_dgm.R")
source("03_methods.R")

library(dplyr)
library(readr)

# ==============================================================================
# 1. SETTINGS
# ==============================================================================

verification_nsim <- 20L
bootstrap_B <- 250L

verification_dir <- file.path(
  simulation_dir,
  "sandwich_bootstrap_verification"
)

dir.create(
  verification_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

if (!file.exists(calibration_rds)) {
  stop("Calibration file not found. Run 01_calibrate_midus.R first.")
}

calibration <- readRDS(calibration_rds)

# The stress setting is used because the variance discrepancy was clearest there.
scenario <- list(
  name = "Stress",
  N = 500L,
  phase2_fraction = 0.25,
  r2_a_marker = 0.10,
  r2_y_marker = 0.10,
  theta = 0.15
)

# ==============================================================================
# 2. GENERIC CENTRAL-DIFFERENCE JACOBIAN
#
# Returns d U(theta) / d theta' for the SUMMED estimating equations.
# ==============================================================================

central_jacobian <- function(
  fn,
  par,
  rel_step = 1e-5
) {

  f0 <- fn(par)

  p <- length(par)
  q <- length(f0)

  J <- matrix(
    NA_real_,
    nrow = q,
    ncol = p
  )

  for (j in seq_len(p)) {

    h <- rel_step * max(
      1,
      abs(par[j])
    )

    par_plus <- par
    par_minus <- par

    par_plus[j] <- par_plus[j] + h
    par_minus[j] <- par_minus[j] - h

    f_plus <- fn(par_plus)
    f_minus <- fn(par_minus)

    J[, j] <- (
      f_plus - f_minus
    ) / (2 * h)
  }

  J
}

# ==============================================================================
# 3. IPW: RECONSTRUCT FIT + NUMERICAL STACKED SANDWICH
# ==============================================================================

ipw_numerical_sandwich <- function(
  dat,
  marker_names
) {

  S <- dat$phase2
  N <- nrow(dat)

  Z <- selection_matrix(
    dat,
    include_y = TRUE
  )

  sel_fit <- glm.fit(
    x = Z,
    y = S,
    family = binomial()
  )

  alpha_hat <- sel_fit$coefficients

  if (any(!is.finite(alpha_hat))) {
    stop("Non-finite coefficient in IPW selection model.")
  }

  pi_hat <- as.vector(
    plogis(
      Z %*% alpha_hat
    )
  )

  obs <- which(S == 1)

  C_obs <- model.matrix(
    full_formula(marker_names),
    data = dat[obs, , drop = FALSE]
  )

  Y_obs <- dat$Y[obs]

  beta_hat <- solve(
    crossprod(
      C_obs,
      C_obs * (1 / pi_hat[obs])
    ),
    crossprod(
      C_obs,
      Y_obs / pi_hat[obs]
    )
  )

  beta_hat <- as.vector(beta_hat)
  names(beta_hat) <- colnames(C_obs)

  q <- ncol(Z)
  d <- ncol(C_obs)

  par_hat <- c(
    alpha_hat,
    beta_hat
  )

  # ---------------------------------------------------------------------------
  # Unit-level stacked estimating equations
  # ---------------------------------------------------------------------------

  psi_matrix <- function(par) {

    alpha <- par[
      seq_len(q)
    ]

    beta <- par[
      q + seq_len(d)
    ]

    pi <- as.vector(
      plogis(
        Z %*% alpha
      )
    )

    psi_alpha <- Z * (
      S - pi
    )

    psi_beta <- matrix(
      0,
      nrow = N,
      ncol = d
    )

    resid_obs <- Y_obs -
      as.vector(
        C_obs %*% beta
      )

    psi_beta[obs, ] <- (
      C_obs * resid_obs
    ) / pi[obs]

    cbind(
      psi_alpha,
      psi_beta
    )
  }

  summed_ee <- function(par) {
    colSums(
      psi_matrix(par)
    )
  }

  root_check <- max(
    abs(
      summed_ee(par_hat)
    )
  )

  Psi_hat <- psi_matrix(
    par_hat
  )

  meat <- crossprod(
    Psi_hat
  )

  J_num <- central_jacobian(
    summed_ee,
    par_hat
  )

  J_inv <- solve(
    J_num
  )

  V_num <- J_inv %*%
    meat %*%
    t(J_inv)

  beta_indices <- q +
    seq_len(d)

  V_beta_num <- V_num[
    beta_indices,
    beta_indices,
    drop = FALSE
  ]

  A_index <- match(
    "A",
    names(beta_hat)
  )

  list(
    estimate = beta_hat[A_index],
    numerical_se = sqrt(
      V_beta_num[
        A_index,
        A_index
      ]
    ),
    root_check = root_check,
    parameter_count = length(par_hat)
  )
}

# ==============================================================================
# 4. AIPW: HELPERS TO PACK / UNPACK SIGMA
# ==============================================================================

sigma_to_vech <- function(
  Sigma,
  pairs
) {

  vapply(
    pairs,
    function(ab) {
      Sigma[
        ab[1],
        ab[2]
      ]
    },
    numeric(1)
  )
}

vech_to_sigma <- function(
  values,
  pairs,
  k
) {

  Sigma <- matrix(
    0,
    nrow = k,
    ncol = k
  )

  for (j in seq_along(pairs)) {

    a <- pairs[[j]][1]
    b <- pairs[[j]][2]

    Sigma[a, b] <- values[j]
    Sigma[b, a] <- values[j]
  }

  Sigma
}

# ==============================================================================
# 5. AIPW: RECONSTRUCT FIT + NUMERICAL STACKED SANDWICH
# ==============================================================================

aipw_numerical_sandwich <- function(
  dat,
  marker_names
) {

  S <- dat$phase2
  N <- nrow(dat)
  obs <- which(S == 1)

  # ---------------------------------------------------------------------------
  # Selection nuisance model
  # ---------------------------------------------------------------------------

  Zs <- selection_matrix(
    dat,
    include_y = TRUE
  )

  sel_fit <- glm.fit(
    x = Zs,
    y = S,
    family = binomial()
  )

  alpha_hat <- sel_fit$coefficients

  if (any(!is.finite(alpha_hat))) {
    stop("Non-finite coefficient in AIPW selection model.")
  }

  pi_hat <- as.vector(
    plogis(
      Zs %*% alpha_hat
    )
  )

  h_hat <- S / pi_hat

  # ---------------------------------------------------------------------------
  # Marker nuisance model
  # ---------------------------------------------------------------------------

  Zm <- marker_predictor_matrix(
    dat,
    include_y = TRUE
  )

  M_obs <- as.matrix(
    dat[
      obs,
      marker_names,
      drop = FALSE
    ]
  )

  Zm_obs <- Zm[
    obs,
    ,
    drop = FALSE
  ]

  Bhat <- qr.solve(
    Zm_obs,
    M_obs
  )

  mu_hat <- Zm %*% Bhat
  colnames(mu_hat) <- marker_names

  E_obs <- M_obs -
    Zm_obs %*% Bhat

  Sigma_hat <- crossprod(
    E_obs
  ) / length(obs)

  # ---------------------------------------------------------------------------
  # Target parameter
  # ---------------------------------------------------------------------------

  W <- target_w_matrix(
    dat
  )

  r <- ncol(W)
  k <- length(marker_names)
  d <- r + k

  Cbar_hat <- cbind(
    W,
    mu_hat
  )

  Cobs_hat <- Cbar_hat

  Cobs_hat[
    obs,
    (r + 1):d
  ] <- M_obs

  Vsig_hat <- matrix(
    0,
    nrow = d,
    ncol = d
  )

  Vsig_hat[
    (r + 1):d,
    (r + 1):d
  ] <- Sigma_hat

  sqrt_h <- sqrt(
    h_hat
  )

  A_sum <- crossprod(
    Cbar_hat
  ) +
    N * Vsig_hat +
    crossprod(
      Cobs_hat * sqrt_h
    ) -
    crossprod(
      Cbar_hat * sqrt_h
    ) -
    sum(h_hat) * Vsig_hat

  b_sum <- colSums(
    Cbar_hat * dat$Y
  ) +
    colSums(
      (
        Cobs_hat - Cbar_hat
      ) *
        (
          h_hat * dat$Y
        )
    )

  beta_hat <- solve(
    A_sum,
    b_sum
  )

  beta_hat <- as.vector(
    beta_hat
  )

  names(beta_hat) <- colnames(
    Cbar_hat
  )

  # ---------------------------------------------------------------------------
  # Parameter-vector layout
  # ---------------------------------------------------------------------------

  qs <- ncol(Zs)
  qm <- ncol(Zm)
  pB <- qm * k

  pairs <- vech_pairs(k)
  hs <- length(pairs)

  sigma_hat_vech <- sigma_to_vech(
    Sigma_hat,
    pairs
  )

  par_hat <- c(
    alpha_hat,
    as.vector(Bhat),
    sigma_hat_vech,
    beta_hat
  )

  i_alpha <- seq_len(qs)

  i_B <- qs +
    seq_len(pB)

  i_Sigma <- qs +
    pB +
    seq_len(hs)

  i_beta <- qs +
    pB +
    hs +
    seq_len(d)

  # ---------------------------------------------------------------------------
  # Unit-level full stacked equations at arbitrary parameter vector
  # ---------------------------------------------------------------------------

  psi_matrix <- function(par) {

    alpha <- par[
      i_alpha
    ]

    B_vec <- par[
      i_B
    ]

    sigma_vec <- par[
      i_Sigma
    ]

    beta <- par[
      i_beta
    ]

    B <- matrix(
      B_vec,
      nrow = qm,
      ncol = k
    )

    Sigma <- vech_to_sigma(
      sigma_vec,
      pairs = pairs,
      k = k
    )

    pi <- as.vector(
      plogis(
        Zs %*% alpha
      )
    )

    h <- S / pi

    mu <- Zm %*% B
    colnames(mu) <- marker_names

    E_all <- matrix(
      0,
      nrow = N,
      ncol = k
    )

    E_all[
      obs,
    ] <- M_obs -
      Zm_obs %*% B

    # ----------------------------
    # Selection score
    # ----------------------------

    psi_alpha <- Zs * (
      S - pi
    )

    # ----------------------------
    # Marker regression score
    # ----------------------------

    psi_B_blocks <- lapply(
      seq_len(k),
      function(j) {
        Zm * (
          S *
            E_all[, j]
        )
      }
    )

    psi_B <- do.call(
      cbind,
      psi_B_blocks
    )

    # ----------------------------
    # Marker covariance moments
    # ----------------------------

    psi_Sigma <- matrix(
      0,
      nrow = N,
      ncol = hs
    )

    for (j in seq_along(pairs)) {

      a <- pairs[[j]][1]
      b <- pairs[[j]][2]

      psi_Sigma[, j] <- S * (
        E_all[, a] *
          E_all[, b] -
          Sigma[a, b]
      )
    }

    # ----------------------------
    # AIPW target equation
    # ----------------------------

    Cbar <- cbind(
      W,
      mu
    )

    Cobs <- Cbar

    Cobs[
      obs,
      (r + 1):d
    ] <- M_obs

    beta_w <- beta[
      seq_len(r)
    ]

    beta_m <- beta[
      (r + 1):d
    ]

    resid_bar <- dat$Y -
      as.vector(
        W %*% beta_w
      ) -
      as.vector(
        mu %*% beta_m
      )

    m_top <- W *
      resid_bar

    sigma_beta_m <- as.vector(
      Sigma %*%
        beta_m
    )

    m_bottom <- sweep(
      mu * resid_bar,
      MARGIN = 2,
      STATS = sigma_beta_m,
      FUN = "-"
    )

    m_psi <- cbind(
      m_top,
      m_bottom
    )

    resid_complete <- dat$Y -
      as.vector(
        Cobs %*% beta
      )

    u_psi <- Cobs *
      resid_complete

    psi_beta <- m_psi +
      (
        u_psi -
          m_psi
      ) * h

    cbind(
      psi_alpha,
      psi_B,
      psi_Sigma,
      psi_beta
    )
  }

  summed_ee <- function(par) {
    colSums(
      psi_matrix(par)
    )
  }

  root_check <- max(
    abs(
      summed_ee(
        par_hat
      )
    )
  )

  Psi_hat <- psi_matrix(
    par_hat
  )

  meat <- crossprod(
    Psi_hat
  )

  J_num <- central_jacobian(
    summed_ee,
    par_hat
  )

  J_inv <- solve(
    J_num
  )

  V_num <- J_inv %*%
    meat %*%
    t(J_inv)

  V_beta_num <- V_num[
    i_beta,
    i_beta,
    drop = FALSE
  ]

  A_index <- match(
    "A",
    names(beta_hat)
  )

  list(
    estimate = beta_hat[A_index],
    numerical_se = sqrt(
      V_beta_num[
        A_index,
        A_index
      ]
    ),
    root_check = root_check,
    parameter_count = length(par_hat)
  )
}

# ==============================================================================
# 6. NONPARAMETRIC BOOTSTRAP SE FOR ONE METHOD
# ==============================================================================

bootstrap_method_se <- function(
  dat,
  marker_names,
  method = c(
    "IPW",
    "AIPW"
  ),
  B = bootstrap_B
) {

  method <- match.arg(
    method
  )

  N <- nrow(dat)

  estimates <- rep(
    NA_real_,
    B
  )

  for (b in seq_len(B)) {

    idx <- sample.int(
      N,
      size = N,
      replace = TRUE
    )

    db <- dat[
      idx,
      ,
      drop = FALSE
    ]

    # Keep factor levels even if a bootstrap resample omits a category.
    db$sex <- factor(
      db$sex,
      levels = levels(dat$sex)
    )

    db$race_eth <- factor(
      db$race_eth,
      levels = levels(dat$race_eth)
    )

    fit_b <- tryCatch(
      {
        if (method == "IPW") {

          fit_ipw(
            dat = db,
            marker_names = marker_names
          )

        } else {

          fit_aipw(
            dat = db,
            marker_names = marker_names,
            selection_include_y = TRUE,
            marker_include_y = TRUE,
            method_label = "AIPW"
          )
        }
      },
      error = function(e) {
        NULL
      }
    )

    if (
      !is.null(fit_b) &&
        nrow(fit_b) == 1 &&
        fit_b$status == "ok" &&
        is.finite(fit_b$estimate)
    ) {
      estimates[b] <- fit_b$estimate
    }
  }

  successful <- is.finite(
    estimates
  )

  n_success <- sum(
    successful
  )

  if (
    n_success <
      max(
        50L,
        ceiling(
          0.80 * B
        )
      )
  ) {
    return(
      list(
        bootstrap_se = NA_real_,
        bootstrap_success = n_success,
        bootstrap_failure = B - n_success
      )
    )
  }

  list(
    bootstrap_se = sd(
      estimates[
        successful
      ]
    ),
    bootstrap_success = n_success,
    bootstrap_failure = B - n_success
  )
}

# ==============================================================================
# 7. BUILD THE STRESS DGM
# ==============================================================================

dgm <- build_dgm_parameters(
  calibration = calibration,
  r2_a_marker = scenario$r2_a_marker,
  r2_y_marker = scenario$r2_y_marker,
  theta = scenario$theta,
  phase2_fraction = scenario$phase2_fraction,
  selection = "mar",
  selection_calibration_n = 50000L,
  selection_seed = 41001L,
  marker_error = "mvn"
)

# ==============================================================================
# 8. RUN VERIFICATION DATASETS
# ==============================================================================

RNGkind(
  "L'Ecuyer-CMRG"
)

set.seed(
  rng_seed_master + 9000L
)

results <- vector(
  "list",
  verification_nsim * 2L
)

z <- 1L

for (i in seq_len(verification_nsim)) {

  cat(
    "\nVerification dataset ",
    i,
    " / ",
    verification_nsim,
    "\n",
    sep = ""
  )

  dat_i <- generate_two_phase_data(
    N = scenario$N,
    calibration = calibration,
    dgm_parameters = dgm
  )$observed

  # ---------------------------------------------------------------------------
  # Current analytic fits
  # ---------------------------------------------------------------------------

  current_ipw <- fit_ipw(
    dat = dat_i,
    marker_names = calibration$marker_names
  )

  current_aipw <- fit_aipw(
    dat = dat_i,
    marker_names = calibration$marker_names,
    selection_include_y = TRUE,
    marker_include_y = TRUE,
    method_label = "AIPW"
  )

  # ---------------------------------------------------------------------------
  # Numerical stacked sandwiches
  # ---------------------------------------------------------------------------

  cat(
    "  Numerical IPW sandwich...\n"
  )

  num_ipw <- ipw_numerical_sandwich(
    dat = dat_i,
    marker_names = calibration$marker_names
  )

  cat(
    "  Numerical AIPW sandwich...\n"
  )

  num_aipw <- aipw_numerical_sandwich(
    dat = dat_i,
    marker_names = calibration$marker_names
  )

  # ---------------------------------------------------------------------------
  # Bootstrap
  # ---------------------------------------------------------------------------

  cat(
    "  IPW bootstrap (",
    bootstrap_B,
    ")...\n",
    sep = ""
  )

  boot_ipw <- bootstrap_method_se(
    dat = dat_i,
    marker_names = calibration$marker_names,
    method = "IPW",
    B = bootstrap_B
  )

  cat(
    "  AIPW bootstrap (",
    bootstrap_B,
    ")...\n",
    sep = ""
  )

  boot_aipw <- bootstrap_method_se(
    dat = dat_i,
    marker_names = calibration$marker_names,
    method = "AIPW",
    B = bootstrap_B
  )

  # ---------------------------------------------------------------------------
  # Store
  # ---------------------------------------------------------------------------

  results[[z]] <- data.frame(
    repetition = i,
    method = "IPW",
    theta_true = scenario$theta,
    n_phase2 = sum(
      dat_i$phase2
    ),
    analytic_estimate =
      current_ipw$estimate,
    analytic_se =
      current_ipw$se,
    numerical_estimate =
      num_ipw$estimate,
    numerical_se =
      num_ipw$numerical_se,
    bootstrap_se =
      boot_ipw$bootstrap_se,
    bootstrap_success =
      boot_ipw$bootstrap_success,
    bootstrap_failure =
      boot_ipw$bootstrap_failure,
    numerical_root_check =
      num_ipw$root_check,
    parameter_count =
      num_ipw$parameter_count
  )

  z <- z + 1L

  results[[z]] <- data.frame(
    repetition = i,
    method = "AIPW",
    theta_true = scenario$theta,
    n_phase2 = sum(
      dat_i$phase2
    ),
    analytic_estimate =
      current_aipw$estimate,
    analytic_se =
      current_aipw$se,
    numerical_estimate =
      num_aipw$estimate,
    numerical_se =
      num_aipw$numerical_se,
    bootstrap_se =
      boot_aipw$bootstrap_se,
    bootstrap_success =
      boot_aipw$bootstrap_success,
    bootstrap_failure =
      boot_aipw$bootstrap_failure,
    numerical_root_check =
      num_aipw$root_check,
    parameter_count =
      num_aipw$parameter_count
  )

  z <- z + 1L

  # Checkpoint after every dataset.
  partial <- bind_rows(
    results[
      seq_len(z - 1L)
    ]
  )

  saveRDS(
    partial,
    file.path(
      verification_dir,
      "sandwich_bootstrap_verification_checkpoint.rds"
    )
  )
}

verification_results <- bind_rows(
  results
)

write_csv(
  verification_results,
  file.path(
    verification_dir,
    "sandwich_bootstrap_verification_results.csv"
  )
)

saveRDS(
  verification_results,
  file.path(
    verification_dir,
    "sandwich_bootstrap_verification_results.rds"
  )
)

# ==============================================================================
# 9. SUMMARY
# ==============================================================================

verification_summary <- verification_results %>%
  group_by(
    method
  ) %>%
  summarise(
    n_datasets = n(),

    mean_n_phase2 =
      mean(
        n_phase2
      ),

    max_abs_estimate_difference =
      max(
        abs(
          analytic_estimate -
            numerical_estimate
        )
      ),

    mean_analytic_se =
      mean(
        analytic_se
      ),

    mean_numerical_se =
      mean(
        numerical_se
      ),

    mean_bootstrap_se =
      mean(
        bootstrap_se,
        na.rm = TRUE
      ),

    mean_analytic_to_numerical_ratio =
      mean(
        analytic_se /
          numerical_se
      ),

    min_analytic_to_numerical_ratio =
      min(
        analytic_se /
          numerical_se
      ),

    max_analytic_to_numerical_ratio =
      max(
        analytic_se /
          numerical_se
      ),

    mean_bootstrap_to_analytic_ratio =
      mean(
        bootstrap_se /
          analytic_se,
        na.rm = TRUE
      ),

    mean_bootstrap_to_numerical_ratio =
      mean(
        bootstrap_se /
          numerical_se,
        na.rm = TRUE
      ),

    max_numerical_root_check =
      max(
        numerical_root_check
      ),

    bootstrap_failures =
      sum(
        bootstrap_failure
      ),

    .groups = "drop"
  )

# ==============================================================================
# 10. ADD THE 2,000-REPETITION EMPIRICAL-SE BENCHMARK IF AVAILABLE
# ==============================================================================

validation_summary_file <- file.path(
  simulation_dir,
  "variance_validation_2000",
  "variance_validation_performance_summary.csv"
)

if (
  file.exists(
    validation_summary_file
  )
) {

  validation_benchmark <- read_csv(
    validation_summary_file,
    show_col_types = FALSE
  ) %>%
    filter(
      scenario == "Stress",
      method %in% c(
        "IPW",
        "AIPW"
      )
    ) %>%
    select(
      method,
      validation_empse = empse,
      validation_modse = modse,
      validation_se_ratio = se_ratio,
      validation_coverage = coverage
    )

  verification_summary <- verification_summary %>%
    left_join(
      validation_benchmark,
      by = "method"
    ) %>%
    mutate(
      bootstrap_to_validation_empse =
        mean_bootstrap_se /
          validation_empse,
      numerical_to_validation_empse =
        mean_numerical_se /
          validation_empse,
      analytic_to_validation_empse =
        mean_analytic_se /
          validation_empse
    )
}

write_csv(
  verification_summary,
  file.path(
    verification_dir,
    "sandwich_bootstrap_verification_summary.csv"
  )
)

# ==============================================================================
# 11. PRINT
# ==============================================================================

cat(
  "\n\n================ VERIFICATION SUMMARY ================\n"
)

print(
  verification_summary,
  n = Inf,
  width = Inf
)

