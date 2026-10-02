################################################################################
# 11_run_aipw_robustness_canary.R
#
# Small implementation canary for the AIPW double-robustness experiment.
#
# Default:
#   - 20 exact primary datasets from the MIDUS-like setting
#   - 20 exact primary datasets from the stress setting
#   - each dataset analyzed under all four AIPW working-model specifications
#   - 40 generated datasets x 4 AIPW specifications = 160 AIPW fits
#
# The exact primary RNG streams are reused so the "Both correct" AIPW estimate
# and SE can be checked against the already-validated primary simulation output.
################################################################################

source("10_aipw_robustness_helpers.R")

library(dplyr)
library(readr)

# ------------------------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------------------------

robustness_root <- file.path(
  results_dir,
  "robustness",
  "aipw_double_robustness"
)

canary_dir <- file.path(
  robustness_root,
  "canary"
)

dir.create(
  canary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Number of repetitions PER generating setting.
n_canary <- as.integer(
  Sys.getenv(
    "AIPW_ROBUSTNESS_CANARY_N",
    unset = "20"
  )
)

if (
  is.na(n_canary) ||
    n_canary < 1L ||
    n_canary > nsim_primary
) {
  stop(
    "AIPW_ROBUSTNESS_CANARY_N must be between 1 and ",
    nsim_primary,
    "."
  )
}

# ------------------------------------------------------------------------------
# 2. Read calibration
# ------------------------------------------------------------------------------

if (!file.exists(calibration_rds)) {
  stop(
    "Calibration object not found:\n",
    calibration_rds
  )
}

calibration <- readRDS(
  calibration_rds
)

# ------------------------------------------------------------------------------
# 3. Run both generating settings
# ------------------------------------------------------------------------------

result_list <- list()
diag_list <- list()
result_index <- 1L
diag_index <- 1L

for (
  setting_index in seq_len(
    nrow(
      aipw_robustness_settings
    )
  )
) {

  ss <- aipw_robustness_settings[
    setting_index,
    ,
    drop = FALSE
  ]

  cache_file <- robustness_cache_filename(
    phase2_fraction =
      ss$phase2_fraction,
    r2_a_marker =
      ss$r2_a_marker,
    r2_y_marker =
      ss$r2_y_marker,
    theta =
      ss$theta
  )

  if (!file.exists(cache_file)) {
    stop(
      "Required primary DGM cache is missing:\n",
      cache_file,
      "\nRun 04_build_dgm_cache.R first."
    )
  }

  dgm <- readRDS(
    cache_file
  )

  # Existing validated primary estimate file, used only for exact comparison.
  primary_estimate_file <- file.path(
    combined_dir,
    "scenarios",
    sprintf(
      "scenario_%03d_estimates.rds",
      ss$primary_scenario_id
    )
  )

  primary_results <- NULL

  if (file.exists(primary_estimate_file)) {

    primary_results <- readRDS(
      primary_estimate_file
    ) %>%
      filter(
        method == "AIPW",
        repetition <= n_canary
      ) %>%
      select(
        repetition,
        primary_estimate = estimate,
        primary_se = se,
        primary_status = status
      )
  }

  cat(
    "\n============================================================\n",
    "AIPW robustness canary\n",
    "Setting: ", ss$setting, "\n",
    "Primary scenario ID: ", ss$primary_scenario_id, "\n",
    "N: ", ss$N, "\n",
    "Phase-2 fraction: ", ss$phase2_fraction, "\n",
    "R2 A~M|X: ", ss$r2_a_marker, "\n",
    "R2 Y~M|A,X: ", ss$r2_y_marker, "\n",
    "Theta: ", ss$theta, "\n",
    "Canary repetitions: ", n_canary, "\n",
    "============================================================\n",
    sep = ""
  )

  for (
    repetition in seq_len(
      n_canary
    )
  ) {

    states <- get_primary_repetition_states(
      scenario_id =
        ss$primary_scenario_id,
      repetition =
        repetition
    )

    # ------------------------------------------------------------------------
    # Generate exactly the same dataset used in the corresponding primary run.
    # ------------------------------------------------------------------------

    .Random.seed <- states$data_state

    sim <- generate_two_phase_data(
      N = ss$N,
      calibration = calibration,
      dgm_parameters = dgm
    )

    dat <- sim$observed

    # ------------------------------------------------------------------------
    # Fit all four AIPW specifications to this ONE fixed dataset.
    # ------------------------------------------------------------------------

    .Random.seed <- states$method_state

    ans <- run_aipw_robustness_specs(
      dat = dat,
      marker_names =
        calibration$marker_names
    ) %>%
      mutate(
        setting_id =
          ss$setting_id,
        setting =
          ss$setting,
        primary_scenario_id =
          ss$primary_scenario_id,
        repetition =
          repetition,
        N =
          ss$N,
        target_phase2_fraction =
          ss$phase2_fraction,
        realized_phase2_fraction =
          mean(
            dat$phase2
          ),
        n_phase2 =
          sum(
            dat$phase2
          ),
        r2_a_marker =
          ss$r2_a_marker,
        r2_y_marker =
          ss$r2_y_marker,
        theta_true =
          ss$theta,
        min_pi_true =
          min(
            dat$pi_true
          ),
        p01_pi_true =
          unname(
            quantile(
              dat$pi_true,
              0.01
            )
          ),
        median_pi_true =
          median(
            dat$pi_true
          ),
        p99_pi_true =
          unname(
            quantile(
              dat$pi_true,
              0.99
            )
          ),
        max_pi_true =
          max(
            dat$pi_true
          ),
        data_rng_state =
          state_to_string(
            states$data_state
          )
      )

    # Exact primary comparison for the both-correct specification.
    if (!is.null(primary_results)) {

      pp <- primary_results %>%
        filter(
          repetition ==
            !!repetition
        )

      if (nrow(pp) == 1L) {

        ans <- ans %>%
          mutate(
            primary_estimate =
              ifelse(
                specification ==
                  "Both correct",
                pp$primary_estimate,
                NA_real_
              ),
            primary_se =
              ifelse(
                specification ==
                  "Both correct",
                pp$primary_se,
                NA_real_
              ),
            primary_status =
              ifelse(
                specification ==
                  "Both correct",
                pp$primary_status,
                NA_character_
              )
          )

      } else {

        ans <- ans %>%
          mutate(
            primary_estimate =
              NA_real_,
            primary_se =
              NA_real_,
            primary_status =
              NA_character_
          )
      }

    } else {

      ans <- ans %>%
        mutate(
          primary_estimate =
            NA_real_,
          primary_se =
            NA_real_,
          primary_status =
            NA_character_
        )
    }

    ans <- ans %>%
      mutate(
        estimate_difference_from_primary =
          ifelse(
            specification ==
              "Both correct",
            estimate -
              primary_estimate,
            NA_real_
          ),
        se_difference_from_primary =
          ifelse(
            specification ==
              "Both correct",
            se -
              primary_se,
            NA_real_
          )
      )

    result_list[[result_index]] <- ans
    result_index <- result_index + 1L

    # ------------------------------------------------------------------------
    # Quantify how much omitting Y changes the two nuisance components.
    # ------------------------------------------------------------------------

    dd <- tryCatch(
      compute_misspecification_diagnostics(
        dat = dat,
        marker_names =
          calibration$marker_names
      ),
      error = function(e) {
        data.frame(
          mean_abs_pi_difference = NA_real_,
          max_abs_pi_difference = NA_real_,
          rmse_pi_difference = NA_real_,
          cor_pi_correct_misspecified = NA_real_,
          mean_abs_marker_mean_difference = NA_real_,
          max_abs_marker_mean_difference = NA_real_,
          rmse_marker_mean_difference = NA_real_,
          relative_frobenius_sigma_difference = NA_real_,
          pi_correct_min = NA_real_,
          pi_correct_max = NA_real_,
          pi_misspecified_min = NA_real_,
          pi_misspecified_max = NA_real_,
          diagnostic_message = conditionMessage(e)
        )
      }
    )

    if (!("diagnostic_message" %in% names(dd))) {
      dd$diagnostic_message <- NA_character_
    }

    dd <- dd %>%
      mutate(
        setting_id =
          ss$setting_id,
        setting =
          ss$setting,
        primary_scenario_id =
          ss$primary_scenario_id,
        repetition =
          repetition,
        N =
          ss$N,
        target_phase2_fraction =
          ss$phase2_fraction,
        realized_phase2_fraction =
          mean(
            dat$phase2
          ),
        n_phase2 =
          sum(
            dat$phase2
          ),
        theta_true =
          ss$theta
      ) %>%
      relocate(
        setting_id,
        setting,
        primary_scenario_id,
        repetition,
        N,
        target_phase2_fraction,
        realized_phase2_fraction,
        n_phase2,
        theta_true
      )

    diag_list[[diag_index]] <- dd
    diag_index <- diag_index + 1L

    cat(
      "Setting ",
      ss$setting,
      ": completed repetition ",
      repetition,
      " / ",
      n_canary,
      "\n",
      sep = ""
    )
  }
}

# ------------------------------------------------------------------------------
# 4. Save canary outputs
# ------------------------------------------------------------------------------

canary_results <- bind_rows(
  result_list
)

canary_diagnostics <- bind_rows(
  diag_list
)

write_csv(
  canary_results,
  file.path(
    canary_dir,
    "aipw_robustness_canary_results.csv"
  )
)

saveRDS(
  canary_results,
  file.path(
    canary_dir,
    "aipw_robustness_canary_results.rds"
  ),
  compress = "xz"
)

write_csv(
  canary_diagnostics,
  file.path(
    canary_dir,
    "aipw_robustness_canary_misspecification_diagnostics.csv"
  )
)

cat(
  "\nAIPW robustness canary complete.\n",
  "Results:\n  ",
  file.path(
    canary_dir,
    "aipw_robustness_canary_results.csv"
  ),
  "\nDiagnostics:\n  ",
  file.path(
    canary_dir,
    "aipw_robustness_canary_misspecification_diagnostics.csv"
  ),
  "\n",
  sep = ""
)
