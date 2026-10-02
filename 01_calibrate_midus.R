# 01_calibrate_midus.R
#
# Uses the analytic MIDUS data to estimate the empirical quantities needed to
# generate simulated datasets. These quantities describe:
# - the distribution of fully observed phase-1 covariates;
# - the relationships between the RNA markers and phase-1 covariates;
# - the residual covariance among the RNA markers;
# - the directions of marker effects in the exposure and outcome models; and
# - the slopes governing selection into the phase-2 sample.
#
# The resulting calibration object is used by the remaining simulation scripts.
# Standardization places variables on a common scale but does not reproduce
# their observed associations or the correlation among the RNA markers.

# Load shared paths, variables, and numerical settings
source("00_config.R")

library(dplyr)
library(readr)

# Align model coefficients with a required set of terms
align_coef <- function(coef_vector, required_names) {
  aligned <- setNames(rep(0, length(required_names)), required_names)
  common_names <- intersect(names(coef_vector), required_names)
  aligned[common_names] <- coef_vector[common_names]
  aligned
}

# 1. Read and check the analytic data

# Stop if the analytic dataset cannot be found
if (!file.exists(analytic_rds)) {
  stop("Analytic RDS not found: ", analytic_rds)
}

df <- readRDS(analytic_rds)

# Confirm that all variables needed for calibration are available
required_variables <- c(
  exposure_name, calibration_outcome, "age_STD",
  "sex", "race_eth", "phase2", marker_names
)

missing_variables <- setdiff(required_variables, names(df))

if (length(missing_variables) > 0L) {
  stop(
    "The analytic dataset is missing: ",
    paste(missing_variables, collapse = ", ")
  )
}

# Retain the observed coding of sex and race/ethnicity
df <- df %>%
  mutate(
    sex = droplevels(factor(sex)),
    race_eth = droplevels(factor(race_eth))
  )

sex_levels <- levels(df$sex)
race_levels <- levels(df$race_eth)

# Phase-1 variables should be observed for everyone in the analytic sample
phase1_vars <- c(
  exposure_name, calibration_outcome, "age_STD",
  "sex", "race_eth", "phase2"
)

if (anyNA(df[, phase1_vars])) {
  stop("Phase-1 variables contain missing values. Recheck analytic sample code.")
}

# Check that the eight RNA markers are either all observed or all missing
n_marker_observed <- rowSums(!is.na(df[, marker_names, drop = FALSE]))

if (!all(n_marker_observed %in% c(0, length(marker_names)))) {
  stop(
    "RNA markers are not blockwise observed/missing. ",
    "Some participants have only a subset of the 8 markers."
  )
}

# Confirm that phase2 agrees with joint availability of the RNA markers
if (!all(df$phase2 == as.integer(n_marker_observed == length(marker_names)))) {
  stop("phase2 does not agree with joint availability of all 8 RNA markers.")
}

# Restrict marker-based calibration models to the phase-2 sample
df_p2 <- df %>% filter(phase2 == 1)

# 2. Construct the covariate model matrices

# Create identical phase-1 covariate matrices for the full and phase-2 samples
X_all <- model.matrix(x_formula, data = df)
X_p2 <- model.matrix(x_formula, data = df_p2)

x_columns <- colnames(X_all)
x_no_intercept <- setdiff(x_columns, "(Intercept)")

# Factor coding should produce the same columns in both samples
if (!identical(colnames(X_all), colnames(X_p2))) {
  stop("X model-matrix columns differ between Phase 1 and Phase 2.")
}

# 3. Estimate the marker model M | X

# Fit the multivariate linear model M = XB + U in the phase-2 sample
M_p2 <- as.matrix(df_p2[, marker_names, drop = FALSE])

B_marker <- qr.solve(X_p2, M_p2)
marker_residuals <- M_p2 - X_p2 %*% B_marker

# Use the residual degrees of freedom to estimate the marker covariance
df_resid_marker <- nrow(M_p2) - ncol(X_p2)

if (df_resid_marker <= 0L) {
  stop("Not enough Phase-2 observations to estimate marker covariance.")
}

Sigma_marker <- crossprod(marker_residuals) / df_resid_marker

# Add a small ridge only if the covariance matrix is nearly singular
marker_eigenvalues <- eigen(
  Sigma_marker,
  symmetric = TRUE,
  only.values = TRUE
)$values

minimum_eigenvalue <- min(marker_eigenvalues)

if (minimum_eigenvalue <= matrix_tolerance) {
  ridge <- abs(minimum_eigenvalue) + 1e-8
  Sigma_marker <- Sigma_marker + diag(ridge, ncol(Sigma_marker))
  
  warning(
    "A very small ridge was added to Sigma_marker because its minimum ",
    "eigenvalue was ", signif(minimum_eigenvalue, 4), "."
  )
}

# 4. Estimate exposure-model coefficient directions

# Fit A | X, M only to determine empirical coefficient directions
exposure_predictors <- c("age_STD", "sex", "race_eth", marker_names)

exposure_formula <- as.formula(
  paste(exposure_name, "~", paste(exposure_predictors, collapse = " + "))
)

fit_exposure <- lm(exposure_formula, data = df_p2)
coef_exposure <- coef(fit_exposure)

# Retain the empirical X coefficients and normalize the marker coefficients
gamma_X_empirical <- align_coef(coef_exposure, x_no_intercept)
gamma_M_raw <- align_coef(coef_exposure, marker_names)

gamma_M_norm <- sqrt(sum(gamma_M_raw^2))

if (gamma_M_norm < 1e-12) {
  stop("Empirical exposure-marker coefficient direction is essentially zero.")
}

gamma_M_direction <- gamma_M_raw / gamma_M_norm

# 5. Estimate outcome-model coefficient directions

# Fit Y | A, X, M only to determine X coefficients and marker directions
outcome_predictors <- c(
  exposure_name, "age_STD", "sex", "race_eth", marker_names
)

outcome_formula <- as.formula(
  paste(
    calibration_outcome,
    "~",
    paste(outcome_predictors, collapse = " + ")
  )
)

fit_outcome <- lm(outcome_formula, data = df_p2)
coef_outcome <- coef(fit_outcome)

# Retain the empirical X coefficients and normalize the marker coefficients
beta_X_empirical <- align_coef(coef_outcome, x_no_intercept)
beta_M_raw <- align_coef(coef_outcome, marker_names)

beta_M_norm <- sqrt(sum(beta_M_raw^2))

if (beta_M_norm < 1e-12) {
  stop("Empirical outcome-marker coefficient direction is essentially zero.")
}

beta_M_direction <- beta_M_raw / beta_M_norm

# 6. Estimate phase-2 selection slopes

# Model selection into phase 2 using the observed exposure, outcome, and X
selection_predictors <- c(
  exposure_name, calibration_outcome, "age_STD", "sex", "race_eth"
)

selection_formula <- as.formula(
  paste("phase2 ~", paste(selection_predictors, collapse = " + "))
)

fit_selection <- glm(
  selection_formula,
  data = df,
  family = binomial()
)

coef_selection <- coef(fit_selection)

# Extract slopes for the exposure, outcome, and phase-1 covariates
eta_A <- unname(coef_selection[exposure_name])
eta_Y <- unname(coef_selection[calibration_outcome])
eta_X <- align_coef(coef_selection, x_no_intercept)

if (any(!is.finite(c(eta_A, eta_Y, eta_X)))) {
  stop("Non-finite coefficient in empirical Phase-2 selection model.")
}

# 7. Save the calibration object

# Store the quantities needed to generate and diagnose simulated datasets
calibration <- list(
  exposure_name = exposure_name,
  calibration_outcome = calibration_outcome,
  marker_names = marker_names,
  x_formula = x_formula,
  x_columns = x_columns,
  x_no_intercept = x_no_intercept,
  sex_levels = sex_levels,
  race_levels = race_levels,
  
  # Empirical phase-1 covariate distribution used for resampling
  x_empirical = df %>%
    select(age_STD, sex, race_eth),
  
  # Marker model parameters
  B_marker = B_marker,
  Sigma_marker = Sigma_marker,
  marker_residuals = marker_residuals,
  
  # Exposure-model coefficients and marker direction
  gamma_X_empirical = gamma_X_empirical,
  gamma_M_direction = gamma_M_direction,
  
  # Outcome-model coefficients and marker direction
  beta_X_empirical = beta_X_empirical,
  beta_M_direction = beta_M_direction,
  
  # Phase-2 selection slopes
  eta_A = eta_A,
  eta_Y = eta_Y,
  eta_X = eta_X,
  
  # Observed sample characteristics
  n_phase1 = nrow(df),
  n_phase2 = sum(df$phase2 == 1),
  observed_phase2_fraction = mean(df$phase2 == 1),
  
  # Fitted models retained for diagnostic review
  fitted = list(
    exposure = fit_exposure,
    outcome = fit_outcome,
    selection = fit_selection
  )
)

saveRDS(calibration, calibration_rds)

# 8. Write diagnostic outputs

# Summarize the sample sizes, sampling fraction, and marker covariance
calibration_summary <- data.frame(
  quantity = c(
    "Phase-1 N",
    "Phase-2 N",
    "Observed Phase-2 fraction",
    "Minimum eigenvalue of Sigma_M"
  ),
  value = c(
    calibration$n_phase1,
    calibration$n_phase2,
    calibration$observed_phase2_fraction,
    min(eigen(Sigma_marker, symmetric = TRUE, only.values = TRUE)$values)
  )
)

write_csv(
  calibration_summary,
  file.path(calibration_dir, "calibration_summary.csv")
)

# Save the residual covariance matrix among the RNA markers
write_csv(
  as.data.frame(Sigma_marker) %>%
    mutate(marker = marker_names, .before = 1),
  file.path(calibration_dir, "marker_residual_covariance.csv")
)

# Save the normalized marker directions used in the exposure and outcome DGMs
write_csv(
  data.frame(
    marker = marker_names,
    gamma_M_direction = as.numeric(gamma_M_direction),
    beta_M_direction = as.numeric(beta_M_direction)
  ),
  file.path(calibration_dir, "marker_coefficient_directions.csv")
)

# Print the calibration location, summary, and fitted models
cat("\nCalibration saved to:\n", calibration_rds, "\n\n")
print(calibration_summary)

cat("\nExposure model used for calibration:\n")
print(summary(fit_exposure))

cat("\nOutcome model used for calibration:\n")
print(summary(fit_outcome))

cat("\nPhase-2 selection model used for calibration:\n")
print(summary(fit_selection))
