################################################################################
# 44_run_midus_real_world_application.R
#
# Real-world MIDUS application of the six analysis methods used in the final
# two-phase simulation framework.
#
# Outcomes:
#   - standardized GrimAge2
#   - standardized DunedinPACE
#
# Exposure:
#   - standardized everyday discrimination
#
# Phase-1 covariates:
#   - standardized age
#   - sex
#   - race/ethnicity
#
# Phase-2 covariates:
#   - eight standardized RNA leukocyte-marker transcripts
#
# Methods:
#   1. Naive
#   2. CCA
#   3. FCS-MI
#   4. JM-MI
#   5. IPW
#   6. AIPW
#
# Analysis specification:
#   - BMI is not used.
#   - DNAm plate and RNA plate are not used.
#   - The script uses the validated 03_methods.R implementation, so the real-data
#     application matches the final simulation estimators.
#   - Run this script from the repository root, where 00_config.R and
#     03_methods.R are located.
################################################################################

# ==============================================================================
# 1. PATHS
# ==============================================================================

# Use an explicitly supplied authorized MIDUS file when available. Otherwise,
# fall back to the conventional ignored path under data/.
analytic_file <- Sys.getenv(
  "SIM_ANALYTIC_RDS",
  unset = file.path(
    "data",
    "MIDUS_discrimination_analysis.rds"
  )
)

if (!file.exists(analytic_file)) {
  stop(
    "Analytic dataset not found:\n",
    analytic_file,
    "\nProvide an authorized file via SIM_ANALYTIC_RDS or run the ",
    "MIDUS analytic-sample construction script first."
  )
}

required_project_files <- c(
  "00_config.R",
  "03_methods.R"
)

missing_project_files <- required_project_files[
  !file.exists(required_project_files)
]

if (length(missing_project_files) > 0L) {
  stop(
    "Run this script from the repository root.\n",
    "Missing project file(s): ",
    paste(missing_project_files, collapse = ", ")
  )
}

Sys.setenv(
  SIM_ANALYTIC_RDS = normalizePath(
    analytic_file,
    winslash = "/",
    mustWork = TRUE
  )
)

# ==============================================================================
# 2. LOAD THE VALIDATED PRIMARY METHOD IMPLEMENTATION
# ==============================================================================

source("03_methods.R")

library(dplyr)
library(readr)
library(tidyr)

real_world_dir <- file.path(
  results_dir,
  "real_world_application"
)

dir.create(
  real_world_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ==============================================================================
# 3. READ AND VALIDATE THE REAL-WORLD ANALYTIC DATA
# ==============================================================================

df <- readRDS(analytic_rds)

real_marker_names <- c(
  "log2_cd19_STD",
  "log2_cd3d_STD",
  "log2_cd3e_STD",
  "log2_cd4_STD",
  "log2_cd8a_STD",
  "log2_cd14_STD",
  "log2_fcgr3a_STD",
  "log2_ncam1_STD"
)

required_columns <- c(
  "id",
  "discrimination_STD",
  "grimage2_STD",
  "dunedinpace_STD",
  "age_STD",
  "sex",
  "race_eth",
  real_marker_names,
  "phase2"
)

missing_columns <- setdiff(
  required_columns,
  names(df)
)

if (length(missing_columns) > 0L) {
  stop(
    "MIDUS_discrimination_analysis.rds is missing required variable(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

df <- df %>%
  mutate(
    sex = droplevels(factor(sex)),
    race_eth = factor(
      race_eth,
      levels = c(
        "non-Hispanic White",
        "non-Hispanic Black",
        "Other"
      )
    ),
    phase2 = as.integer(phase2)
  )

if (anyNA(df$phase2) || !all(df$phase2 %in% c(0L, 1L))) {
  stop("phase2 must be a complete binary indicator coded 0/1.")
}

phase1_required <- c(
  "discrimination_STD",
  "grimage2_STD",
  "dunedinpace_STD",
  "age_STD",
  "sex",
  "race_eth"
)

phase1_missing <- vapply(
  df[phase1_required],
  function(x) sum(is.na(x)),
  integer(1)
)

if (any(phase1_missing > 0L)) {
  stop(
    "Unexpected missingness remains in Phase-1 variables:\n",
    paste(
      names(phase1_missing),
      phase1_missing,
      sep = " = ",
      collapse = "\n"
    )
  )
}

n_marker_observed <- rowSums(
  !is.na(df[real_marker_names])
)

marker_complete <- as.integer(
  n_marker_observed == length(real_marker_names)
)

phase2_marker_mismatch <- sum(
  df$phase2 != marker_complete
)

if (phase2_marker_mismatch > 0L) {
  stop(
    "The phase2 indicator disagrees with complete observation of all eight ",
    "RNA markers for ",
    phase2_marker_mismatch,
    " participant(s). Recheck the analytic-sample construction."
  )
}

n_partial_marker_rows <- sum(
  n_marker_observed > 0L &
    n_marker_observed < length(real_marker_names)
)

# ==============================================================================
# 4. SAMPLE AND MISSINGNESS DIAGNOSTICS
# ==============================================================================

sample_summary <- tibble::tibble(
  N_phase1 = nrow(df),
  N_phase2 = sum(df$phase2 == 1L),
  phase2_fraction = mean(df$phase2),
  N_not_phase2 = sum(df$phase2 == 0L),
  N_partial_marker_rows = n_partial_marker_rows
)

marker_missingness <- tibble::tibble(
  marker = real_marker_names,
  n_missing = vapply(
    df[real_marker_names],
    function(x) sum(is.na(x)),
    integer(1)
  )
) %>%
  mutate(
    percent_missing = 100 * n_missing / nrow(df)
  )

write_csv(
  sample_summary,
  file.path(
    real_world_dir,
    "midus_real_world_sample_summary.csv"
  )
)

write_csv(
  marker_missingness,
  file.path(
    real_world_dir,
    "midus_real_world_marker_missingness.csv"
  )
)

cat(
  "\n============================================================\n",
  "MIDUS REAL-WORLD APPLICATION\n",
  "============================================================\n",
  "Phase-1 N: ", sample_summary$N_phase1, "\n",
  "Phase-2 N: ", sample_summary$N_phase2, "\n",
  "Phase-2 fraction: ",
  round(sample_summary$phase2_fraction, 4),
  "\n",
  "Rows with partially observed RNA marker blocks: ",
  sample_summary$N_partial_marker_rows,
  "\n",
  "Primary FCS method: ", mice_method_primary, "\n",
  "Imputations: ", nimp_primary, "\n",
  "MICE max iterations: ", mice_maxit_primary, "\n",
  "JOMO burn-in: ", jomo_nburn_primary, "\n",
  "JOMO between imputations: ", jomo_nbetween_primary, "\n",
  "============================================================\n\n",
  sep = ""
)

# ==============================================================================
# 5. PREPARE ONE OUTCOME-SPECIFIC ANALYSIS DATASET
# ==============================================================================

prepare_analysis_data <- function(data, outcome_variable) {

  out <- data %>%
    transmute(
      id = id,
      A = discrimination_STD,
      Y = .data[[outcome_variable]],
      age_STD = age_STD,
      sex = sex,
      race_eth = race_eth,
      across(
        all_of(real_marker_names)
      ),
      phase2 = phase2
    )

  if (
    anyNA(out$A) ||
      anyNA(out$Y) ||
      anyNA(out$age_STD) ||
      anyNA(out$sex) ||
      anyNA(out$race_eth)
  ) {
    stop(
      "Unexpected Phase-1 missingness after preparing outcome: ",
      outcome_variable
    )
  }

  out
}

# ==============================================================================
# 6. RUN ALL SIX METHODS FOR BOTH OUTCOMES
# ==============================================================================

outcome_map <- tibble::tribble(
  ~outcome_variable, ~outcome_label, ~analysis_seed,
  "grimage2_STD",    "GrimAge2",      20260901L,
  "dunedinpace_STD", "DunedinPACE",   20260902L
)

results_list <- vector(
  "list",
  nrow(outcome_map)
)

for (i in seq_len(nrow(outcome_map))) {

  outcome_variable <- outcome_map$outcome_variable[[i]]
  outcome_label <- outcome_map$outcome_label[[i]]
  analysis_seed <- outcome_map$analysis_seed[[i]]

  cat(
    "\n------------------------------------------------------------\n",
    "Outcome: ", outcome_label, "\n",
    "------------------------------------------------------------\n",
    sep = ""
  )

  analysis_data <- prepare_analysis_data(
    data = df,
    outcome_variable = outcome_variable
  )

  set.seed(analysis_seed)

  start_time <- proc.time()[["elapsed"]]

  method_results <- run_all_methods(
    dat = analysis_data,
    marker_names = real_marker_names,
    nimp = nimp_primary,
    mice_maxit = mice_maxit_primary,
    jomo_nburn = jomo_nburn_primary,
    jomo_nbetween = jomo_nbetween_primary
  )

  elapsed_seconds <- (
    proc.time()[["elapsed"]] -
      start_time
  )

  method_results <- method_results %>%
    mutate(
      outcome = outcome_label,
      outcome_variable = outcome_variable,
      exposure = "Everyday discrimination",
      exposure_variable = "discrimination_STD",
      N_phase1 = nrow(analysis_data),
      N_phase2 = sum(analysis_data$phase2 == 1L),
      phase2_fraction = mean(analysis_data$phase2),
      nimp = nimp_primary,
      mice_method = mice_method_primary,
      mice_maxit = mice_maxit_primary,
      jomo_nburn = jomo_nburn_primary,
      jomo_nbetween = jomo_nbetween_primary,
      analysis_seed = analysis_seed,
      total_all_methods_elapsed_seconds = elapsed_seconds
    ) %>%
    relocate(
      outcome,
      outcome_variable,
      method,
      estimate,
      se,
      conf_low,
      conf_high,
      p_value,
      df,
      status,
      message
    )

  results_list[[i]] <- method_results

  print(
    method_results %>%
      select(
        outcome,
        method,
        estimate,
        se,
        conf_low,
        conf_high,
        p_value,
        status
      )
  )

  cat(
    "Total elapsed seconds for ",
    outcome_label,
    ": ",
    round(elapsed_seconds, 2),
    "\n",
    sep = ""
  )
}

real_world_results <- bind_rows(
  results_list
)

# ==============================================================================
# 7. VALIDATE RESULTS
# ==============================================================================

expected_methods <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

method_qc <- real_world_results %>%
  count(
    outcome,
    method,
    name = "n_rows"
  )

missing_method_combinations <- tidyr::expand_grid(
  outcome = outcome_map$outcome_label,
  method = expected_methods
) %>%
  anti_join(
    method_qc,
    by = c("outcome", "method")
  )

duplicate_method_combinations <- method_qc %>%
  filter(n_rows != 1L)

if (nrow(missing_method_combinations) > 0L) {
  stop(
    "One or more outcome-method combinations are missing."
  )
}

if (nrow(duplicate_method_combinations) > 0L) {
  stop(
    "One or more outcome-method combinations are duplicated."
  )
}

result_qc <- real_world_results %>%
  transmute(
    outcome,
    method,
    status,
    finite_estimate = is.finite(estimate),
    finite_se = is.finite(se) & se > 0,
    finite_ci = (
      is.finite(conf_low) &
        is.finite(conf_high)
    ),
    passed = (
      status == "ok" &
        is.finite(estimate) &
        is.finite(se) &
        se > 0 &
        is.finite(conf_low) &
        is.finite(conf_high)
    )
  )

write_csv(
  result_qc,
  file.path(
    real_world_dir,
    "midus_real_world_result_qc.csv"
  )
)

# ==============================================================================
# 8. SAVE PRIMARY REAL-WORLD RESULTS
# ==============================================================================

method_order <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

real_world_results <- real_world_results %>%
  mutate(
    method = factor(
      method,
      levels = method_order
    )
  ) %>%
  arrange(
    outcome,
    method
  )

write_csv(
  real_world_results,
  file.path(
    real_world_dir,
    "midus_real_world_method_results.csv"
  )
)

saveRDS(
  real_world_results,
  file.path(
    real_world_dir,
    "midus_real_world_method_results.rds"
  )
)

# ==============================================================================
# 9. MANUSCRIPT-READY COMPACT TABLE
# ==============================================================================

manuscript_table <- real_world_results %>%
  transmute(
    outcome,
    method = as.character(method),
    estimate,
    standard_error = se,
    ci_95_low = conf_low,
    ci_95_high = conf_high,
    p_value,
    N_phase1,
    N_phase2,
    phase2_fraction,
    status
  )

write_csv(
  manuscript_table,
  file.path(
    real_world_dir,
    "midus_real_world_manuscript_table.csv"
  )
)

# ==============================================================================
# 10. WEIGHT DIAGNOSTICS FOR IPW AND AIPW
# ==============================================================================

weight_columns <- intersect(
  c(
    "weight_min",
    "weight_p99",
    "weight_max",
    "weight_cv",
    "weight_ess"
  ),
  names(real_world_results)
)

if (length(weight_columns) > 0L) {

  weight_diagnostics <- real_world_results %>%
    filter(
      as.character(method) %in% c("IPW", "AIPW")
    ) %>%
    select(
      outcome,
      method,
      all_of(weight_columns)
    )

  write_csv(
    weight_diagnostics,
    file.path(
      real_world_dir,
      "midus_real_world_weight_diagnostics.csv"
    )
  )
}

# ==============================================================================
# 11. FINAL CONSOLE OUTPUT
# ==============================================================================

cat(
  "\n============================================================\n",
  "REAL-WORLD APPLICATION COMPLETE\n",
  "============================================================\n",
  "Output directory:\n  ",
  real_world_dir,
  "\n\n",
  "Primary result file:\n  ",
  file.path(
    real_world_dir,
    "midus_real_world_method_results.csv"
  ),
  "\n\n",
  "All outcome-method fits passed QC: ",
  all(result_qc$passed),
  "\n",
  "============================================================\n\n",
  sep = ""
)

print(manuscript_table)
