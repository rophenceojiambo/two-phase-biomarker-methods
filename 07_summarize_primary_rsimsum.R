# 07_summarize_primary_rsimsum.R
#
# Summarizes the primary MIDUS simulation using the ADEMP framework. Run this
# script after 06_combine_primary_results.R has created one estimate file for
# each of the 324 scenarios.
#
# The rsimsum package is the authoritative source for:
# - average estimate and bias;
# - empirical and model-based standard errors;
# - mean squared error and relative precision;
# - relative error in the model-based standard error;
# - confidence-interval coverage and bias-eliminated coverage; and
# - Monte Carlo standard errors for these performance measures.
#
# This script also calculates:
# - the ratio of model-based to empirical standard error;
# - the null rejection rate when theta = 0;
# - method failure rates;
# - method runtimes;
# - realized phase-2 sample diagnostics; and
# - IPW and AIPW weight diagnostics.
#
# Scenario-level rsimsum objects are retained so individual results can be
# inspected without repeating the full summarization.

# Load shared paths, design settings, and reporting thresholds
source("00_config.R")

library(dplyr)
library(tidyr)
library(readr)
library(rsimsum)

# Define input and output directories
scenario_dir <- file.path(combined_dir, "scenarios")
rsimsum_object_dir <- file.path(summary_dir, "rsimsum_objects")

dir.create(
  rsimsum_object_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# Construct the path to a combined scenario estimate file
scenario_estimate_file <- function(scenario_id) {
  file.path(
    scenario_dir,
    sprintf("scenario_%03d_estimates.rds", scenario_id)
  )
}

# Add the DGM settings to a scenario-level summary
add_scenario_metadata <- function(data, scenario_number, scenario_metadata) {
  data %>%
    mutate(
      scenario_id = scenario_number,
      N = scenario_metadata$N[[1]],
      phase2_fraction = scenario_metadata$phase2_fraction[[1]],
      r2_a_marker = scenario_metadata$r2_a_marker[[1]],
      r2_y_marker = scenario_metadata$r2_y_marker[[1]],
      theta_true = scenario_metadata$theta[[1]]
    ) %>%
    relocate(
      scenario_id,
      N,
      phase2_fraction,
      r2_a_marker,
      r2_y_marker,
      theta_true
    )
}

# Request the ADEMP measures needed for reporting
stats_to_keep <- c(
  "nsim",
  "thetamean",
  "bias",
  "empse",
  "mse",
  "relprec",
  "modelse",
  "relerror",
  "cover",
  "becover",
  "power"
)

# Initialize scenario-level result storage
rsimsum_tidy_list <- vector("list", nrow(primary_grid))
type1_list <- list()
failure_list <- vector("list", nrow(primary_grid))
runtime_list <- vector("list", nrow(primary_grid))
phase2_list <- vector("list", nrow(primary_grid))
weight_list <- vector("list", nrow(primary_grid))

type1_index <- 1L

# Summarize each scenario separately
for (s in seq_len(nrow(primary_grid))) {
  scenario_id <- primary_grid$scenario_id[s]
  
  cat(
    "Summarising scenario ",
    scenario_id,
    " / ",
    nrow(primary_grid),
    " with rsimsum\n",
    sep = ""
  )
  
  scenario_file <- scenario_estimate_file(scenario_id)
  
  if (!file.exists(scenario_file)) {
    stop(
      "Scenario estimate file missing:\n",
      scenario_file,
      "\nRun 06_combine_primary_results.R first."
    )
  }
  
  scenario_data <- readRDS(scenario_file)
  
  scenario_metadata <- primary_grid[
    primary_grid$scenario_id == scenario_id,
    ,
    drop = FALSE
  ]
  
  theta_true <- scenario_metadata$theta[[1]]
  
  # Use method-specific confidence intervals from the analysis functions.
  # This retains t-based Rubin intervals for MI and normal intervals for
  # methods with infinite degrees of freedom.
  rsimsum_object <- rsimsum::simsum(
    data = scenario_data,
    estvarname = "estimate",
    true = theta_true,
    se = "se",
    methodvar = "method",
    ref = "CCA",
    ci.limits = c("conf_low", "conf_high"),
    dropbig = FALSE,
    x = FALSE,
    control = list(
      mcse = TRUE,
      level = nominal_level,
      na.rm = TRUE
    )
  )
  
  # Retain the complete rsimsum object for scenario-specific diagnostics
  saveRDS(
    rsimsum_object,
    file.path(
      rsimsum_object_dir,
      sprintf("scenario_%03d_rsimsum.rds", scenario_id)
    )
  )
  
  rsimsum_summary <- summary(
    rsimsum_object,
    stats = stats_to_keep
  )
  
  rsimsum_tidy_list[[s]] <- generics::tidy(rsimsum_summary) %>%
    add_scenario_metadata(
      scenario_number = scenario_id,
      scenario_metadata = scenario_metadata
    )
  
  # Calculate null rejection rates only for scenarios with theta = 0
  if (isTRUE(all.equal(theta_true, 0))) {
    type1_rows <- scenario_data %>%
      group_by(method) %>%
      summarise(
        nsim_total = n(),
        nsim_success = sum(
          status == "ok" & is.finite(p_value)
        ),
        null_rejection_rate = mean(
          p_value[
            status == "ok" & is.finite(p_value)
          ] < alpha_level
        ),
        mcse_null_rejection = sqrt(
          null_rejection_rate *
            (1 - null_rejection_rate) /
            nsim_success
        ),
        .groups = "drop"
      ) %>%
      add_scenario_metadata(
        scenario_number = scenario_id,
        scenario_metadata = scenario_metadata
      )
    
    type1_list[[type1_index]] <- type1_rows
    type1_index <- type1_index + 1L
  }
  
  # Calculate the proportion of repetitions in which each method failed
  failure_list[[s]] <- scenario_data %>%
    group_by(method) %>%
    summarise(
      nsim_total = n(),
      n_failed = sum(
        status != "ok" |
          !is.finite(estimate) |
          !is.finite(se)
      ),
      failure_rate = n_failed / nsim_total,
      mcse_failure = sqrt(
        failure_rate *
          (1 - failure_rate) /
          nsim_total
      ),
      .groups = "drop"
    ) %>%
    add_scenario_metadata(
      scenario_number = scenario_id,
      scenario_metadata = scenario_metadata
    )
  
  # Summarize the elapsed runtime of each method
  runtime_list[[s]] <- scenario_data %>%
    group_by(method) %>%
    summarise(
      mean_runtime_seconds = mean(
        elapsed_seconds,
        na.rm = TRUE
      ),
      median_runtime_seconds = median(
        elapsed_seconds,
        na.rm = TRUE
      ),
      p95_runtime_seconds = unname(
        quantile(
          elapsed_seconds,
          0.95,
          na.rm = TRUE
        )
      ),
      total_cpu_hours = sum(
        elapsed_seconds,
        na.rm = TRUE
      ) / 3600,
      .groups = "drop"
    ) %>%
    add_scenario_metadata(
      scenario_number = scenario_id,
      scenario_metadata = scenario_metadata
    )
  
  # Sampling quantities are repeated for all six methods, so retain one row
  # per repetition before calculating phase-2 diagnostics
  phase2_list[[s]] <- scenario_data %>%
    distinct(
      repetition,
      realized_phase2_fraction,
      n_phase2,
      min_pi_true,
      p01_pi_true,
      median_pi_true,
      p99_pi_true,
      max_pi_true
    ) %>%
    summarise(
      mean_n_phase2 = mean(n_phase2),
      sd_n_phase2 = sd(n_phase2),
      min_n_phase2 = min(n_phase2),
      max_n_phase2 = max(n_phase2),
      mean_realized_phase2_fraction = mean(
        realized_phase2_fraction
      ),
      mean_min_pi_true = mean(min_pi_true),
      mean_p01_pi_true = mean(p01_pi_true),
      mean_median_pi_true = mean(median_pi_true),
      mean_p99_pi_true = mean(p99_pi_true),
      mean_max_pi_true = mean(max_pi_true)
    ) %>%
    add_scenario_metadata(
      scenario_number = scenario_id,
      scenario_metadata = scenario_metadata
    )
  
  # Summarize inverse-probability weights for successful IPW and AIPW fits
  weight_list[[s]] <- scenario_data %>%
    filter(
      method %in% c("IPW", "AIPW"),
      status == "ok"
    ) %>%
    group_by(method) %>%
    summarise(
      mean_weight_min = mean(
        weight_min,
        na.rm = TRUE
      ),
      mean_weight_p99 = mean(
        weight_p99,
        na.rm = TRUE
      ),
      mean_weight_max = mean(
        weight_max,
        na.rm = TRUE
      ),
      mean_weight_cv = mean(
        weight_cv,
        na.rm = TRUE
      ),
      p95_weight_cv = unname(
        quantile(
          weight_cv,
          0.95,
          na.rm = TRUE
        )
      ),
      mean_weight_ess = mean(
        weight_ess,
        na.rm = TRUE
      ),
      min_weight_ess = min(
        weight_ess,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    add_scenario_metadata(
      scenario_number = scenario_id,
      scenario_metadata = scenario_metadata
    )
}

# 1. Save the tidy rsimsum output

rsimsum_tidy <- bind_rows(rsimsum_tidy_list)

write_csv(
  rsimsum_tidy,
  file.path(summary_dir, "primary_rsimsum_tidy.csv")
)

saveRDS(
  rsimsum_tidy,
  file.path(summary_dir, "primary_rsimsum_tidy.rds"),
  compress = "xz"
)

# 2. Create the scenario-by-method performance table

performance <- rsimsum_tidy %>%
  mutate(stat = as.character(stat))

performance_wide <- performance %>%
  select(
    scenario_id,
    N,
    phase2_fraction,
    r2_a_marker,
    r2_y_marker,
    theta_true,
    method,
    stat,
    est,
    mcse
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = c(est, mcse),
    names_sep = "__"
  )

# Convert rsimsum's percentage error in model SE to ModelSE / EmpSE
if ("est__relerror" %in% names(performance_wide)) {
  performance_wide <- performance_wide %>%
    mutate(
      se_ratio = 1 + est__relerror / 100,
      mcse_se_ratio = mcse__relerror / 100
    )
}

write_csv(
  performance_wide,
  file.path(summary_dir, "primary_performance_wide.csv")
)

# 3. Save the null rejection-rate summary

type1_summary <- if (length(type1_list) == 0L) {
  tibble::tibble()
} else {
  bind_rows(type1_list)
}

write_csv(
  type1_summary,
  file.path(summary_dir, "primary_type1_error.csv")
)

# 4. Save method failure rates

failure_summary <- bind_rows(failure_list)

write_csv(
  failure_summary,
  file.path(summary_dir, "primary_failure_rates.csv")
)

# 5. Save method runtime summaries

runtime_summary <- bind_rows(runtime_list)

write_csv(
  runtime_summary,
  file.path(summary_dir, "primary_runtime.csv")
)

# 6. Save phase-2 sampling diagnostics

phase2_summary <- bind_rows(phase2_list)

write_csv(
  phase2_summary,
  file.path(summary_dir, "primary_phase2_diagnostics.csv")
)

# 7. Save IPW and AIPW weight diagnostics

weight_summary <- bind_rows(weight_list)

write_csv(
  weight_summary,
  file.path(summary_dir, "primary_weight_diagnostics.csv")
)

cat(
  "\nPrimary rsimsum analysis complete.\n",
  "Main files:\n",
  file.path(summary_dir, "primary_rsimsum_tidy.csv"),
  "\n",
  file.path(summary_dir, "primary_performance_wide.csv"),
  "\n",
  file.path(summary_dir, "primary_type1_error.csv"),
  "\n",
  file.path(summary_dir, "primary_failure_rates.csv"),
  "\n",
  sep = ""
)
