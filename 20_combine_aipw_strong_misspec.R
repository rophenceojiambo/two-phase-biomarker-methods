################################################################################
# 20_combine_aipw_strong_misspec.R
#
# Validates and combines the two strong-misspecification production files.
#
# Each file must contain all 2,000 repetitions, all four AIPW specifications,
# one diagnostic row per repetition, and exact primary AIPW reproduction.
################################################################################

source("16_aipw_strong_misspec_helpers.R")

library(dplyr)
library(readr)

# Set the production and combined-result directories.
sensitivity_root <- file.path(
  results_dir,
  "sensitivity",
  "aipw_strong_misspecification"
)

production_dir <- file.path(
  sensitivity_root,
  "production"
)

combined_dir_strong <- file.path(
  sensitivity_root,
  "combined"
)

dir.create(
  combined_dir_strong,
  recursive = TRUE,
  showWarnings = FALSE
)

n_settings <- nrow(aipw_strong_settings)
expected_repetitions <- seq_len(nsim_primary)

expected_specifications <- sort(
  aipw_strong_specs$specification
)

expected_specification_string <- paste(
  expected_specifications,
  collapse = " | "
)

exact_tolerance <- 1e-10

production_objects <- vector("list", n_settings)
manifest_rows <- vector("list", n_settings)

# Validate each production file before combining results.
for (i in seq_len(n_settings)) {
  setting <- aipw_strong_settings[
    i,
    ,
    drop = FALSE
  ]
  
  setting_slug <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    setting$setting
  )
  
  setting_dir <- file.path(
    production_dir,
    sprintf(
      "setting_%02d_%s",
      setting$setting_id,
      setting_slug
    )
  )
  
  setting_file <- file.path(
    setting_dir,
    sprintf(
      "aipw_strong_setting_%02d_results.rds",
      setting$setting_id
    )
  )
  
  if (!file.exists(setting_file)) {
    manifest_rows[[i]] <- data.frame(
      setting_id = setting$setting_id,
      setting = setting$setting,
      file_found = FALSE,
      repetitions_found = NA_integer_,
      result_rows = NA_integer_,
      diagnostic_rows = NA_integer_,
      n_failures = NA_integer_,
      both_correct_primary_exact = FALSE,
      complete = FALSE,
      problem = "Production file missing"
    )
    
    next
  }
  
  production_object <- readRDS(setting_file)
  
  if (
    !is.list(production_object) ||
    !all(c("results", "diagnostics") %in%
         names(production_object))
  ) {
    manifest_rows[[i]] <- data.frame(
      setting_id = setting$setting_id,
      setting = setting$setting,
      file_found = TRUE,
      repetitions_found = NA_integer_,
      result_rows = NA_integer_,
      diagnostic_rows = NA_integer_,
      n_failures = NA_integer_,
      both_correct_primary_exact = FALSE,
      complete = FALSE,
      problem =
        "Production object does not contain results and diagnostics"
    )
    
    next
  }
  
  setting_results <- production_object$results
  setting_diagnostics <- production_object$diagnostics
  
  if (
    !is.data.frame(setting_results) ||
    !is.data.frame(setting_diagnostics)
  ) {
    manifest_rows[[i]] <- data.frame(
      setting_id = setting$setting_id,
      setting = setting$setting,
      file_found = TRUE,
      repetitions_found = NA_integer_,
      result_rows = if (
        is.data.frame(setting_results)
      ) {
        nrow(setting_results)
      } else {
        NA_integer_
      },
      diagnostic_rows = if (
        is.data.frame(setting_diagnostics)
      ) {
        nrow(setting_diagnostics)
      } else {
        NA_integer_
      },
      n_failures = NA_integer_,
      both_correct_primary_exact = FALSE,
      complete = FALSE,
      problem = "Results or diagnostics are not data frames"
    )
    
    next
  }
  
  required_result_columns <- c(
    "setting_id",
    "setting",
    "repetition",
    "specification",
    "status",
    "estimate",
    "se",
    "estimate_difference_from_primary",
    "se_difference_from_primary"
  )
  
  required_diagnostic_columns <- c(
    "setting_id",
    "setting",
    "repetition"
  )
  
  missing_result_columns <- setdiff(
    required_result_columns,
    names(setting_results)
  )
  
  missing_diagnostic_columns <- setdiff(
    required_diagnostic_columns,
    names(setting_diagnostics)
  )
  
  problems <- character()
  
  if (length(missing_result_columns) > 0L) {
    problems <- c(
      problems,
      paste(
        "missing result columns:",
        paste(missing_result_columns, collapse = ", ")
      )
    )
  }
  
  if (length(missing_diagnostic_columns) > 0L) {
    problems <- c(
      problems,
      paste(
        "missing diagnostic columns:",
        paste(missing_diagnostic_columns, collapse = ", ")
      )
    )
  }
  
  if (
    length(missing_result_columns) > 0L ||
    length(missing_diagnostic_columns) > 0L
  ) {
    manifest_rows[[i]] <- data.frame(
      setting_id = setting$setting_id,
      setting = setting$setting,
      file_found = TRUE,
      repetitions_found = NA_integer_,
      result_rows = nrow(setting_results),
      diagnostic_rows = nrow(setting_diagnostics),
      n_failures = NA_integer_,
      both_correct_primary_exact = FALSE,
      complete = FALSE,
      problem = paste(problems, collapse = "; ")
    )
    
    next
  }
  
  result_repetitions <- sort(
    unique(setting_results$repetition)
  )
  
  diagnostic_repetitions <- sort(
    unique(setting_diagnostics$repetition)
  )
  
  # Check one result per repetition and specification.
  duplicate_result_keys <- setting_results %>%
    count(
      setting_id,
      repetition,
      specification,
      name = "n"
    ) %>%
    filter(.data$n != 1L)
  
  specification_sets <- setting_results %>%
    group_by(repetition) %>%
    summarise(
      n_specifications = n_distinct(specification),
      specifications = paste(
        sort(unique(specification)),
        collapse = " | "
      ),
      .groups = "drop"
    )
  
  incorrect_specification_sets <- specification_sets %>%
    filter(
      .data$n_specifications !=
        length(expected_specifications) |
        .data$specifications !=
        expected_specification_string
    )
  
  # Check one diagnostic row per repetition.
  duplicate_diagnostic_keys <- setting_diagnostics %>%
    count(setting_id, repetition, name = "n") %>%
    filter(.data$n != 1L)
  
  result_setting_correct <- all(
    !is.na(setting_results$setting_id) &
      setting_results$setting_id == setting$setting_id &
      !is.na(setting_results$setting) &
      setting_results$setting == setting$setting
  )
  
  diagnostic_setting_correct <- all(
    !is.na(setting_diagnostics$setting_id) &
      setting_diagnostics$setting_id ==
      setting$setting_id &
      !is.na(setting_diagnostics$setting) &
      setting_diagnostics$setting ==
      setting$setting
  )
  
  # Confirm exact agreement with the primary AIPW results.
  both_correct <- setting_results %>%
    filter(.data$specification == "Both correct")
  
  primary_exact <- (
    nrow(both_correct) == nsim_primary &&
      setequal(
        both_correct$repetition,
        expected_repetitions
      ) &&
      all(
        is.finite(
          both_correct$estimate_difference_from_primary
        )
      ) &&
      all(
        is.finite(
          both_correct$se_difference_from_primary
        )
      ) &&
      max(
        abs(
          both_correct$estimate_difference_from_primary
        )
      ) < exact_tolerance &&
      max(
        abs(
          both_correct$se_difference_from_primary
        )
      ) < exact_tolerance
  )
  
  failures <- setting_results %>%
    filter(
      is.na(status) |
        status != "ok" |
        !is.finite(estimate) |
        !is.finite(se)
    )
  
  if (
    !setequal(
      result_repetitions,
      expected_repetitions
    )
  ) {
    problems <- c(
      problems,
      "result repetition set incomplete"
    )
  }
  
  if (
    !setequal(
      diagnostic_repetitions,
      expected_repetitions
    )
  ) {
    problems <- c(
      problems,
      "diagnostic repetition set incomplete"
    )
  }
  
  if (
    nrow(setting_results) !=
    nsim_primary * length(expected_specifications)
  ) {
    problems <- c(
      problems,
      "unexpected result-row count"
    )
  }
  
  if (nrow(setting_diagnostics) != nsim_primary) {
    problems <- c(
      problems,
      "unexpected diagnostic-row count"
    )
  }
  
  if (nrow(duplicate_result_keys) > 0L) {
    problems <- c(
      problems,
      "duplicate result keys"
    )
  }
  
  if (nrow(incorrect_specification_sets) > 0L) {
    problems <- c(
      problems,
      "incorrect specification set"
    )
  }
  
  if (nrow(duplicate_diagnostic_keys) > 0L) {
    problems <- c(
      problems,
      "duplicate diagnostic keys"
    )
  }
  
  if (!result_setting_correct) {
    problems <- c(
      problems,
      "incorrect setting in results"
    )
  }
  
  if (!diagnostic_setting_correct) {
    problems <- c(
      problems,
      "incorrect setting in diagnostics"
    )
  }
  
  if (!primary_exact) {
    problems <- c(
      problems,
      "both-correct does not reproduce primary AIPW"
    )
  }
  
  manifest_rows[[i]] <- data.frame(
    setting_id = setting$setting_id,
    setting = setting$setting,
    file_found = TRUE,
    repetitions_found = length(result_repetitions),
    result_rows = nrow(setting_results),
    diagnostic_rows = nrow(setting_diagnostics),
    n_failures = nrow(failures),
    both_correct_primary_exact = primary_exact,
    complete = length(problems) == 0L,
    problem = paste(problems, collapse = "; ")
  )
  
  production_objects[[i]] <- production_object
}

# Save the completion report before stopping for missing or invalid tasks.
completion_manifest <- bind_rows(manifest_rows)

write_csv(
  completion_manifest,
  file.path(
    combined_dir_strong,
    "aipw_strong_completion_manifest.csv"
  )
)

if (any(!completion_manifest$complete)) {
  print(completion_manifest, n = Inf)
  
  stop(
    "Strong AIPW sensitivity production is incomplete. ",
    "See aipw_strong_completion_manifest.csv."
  )
}

# Combine results only after both settings pass validation.
all_results <- bind_rows(
  lapply(production_objects, function(x) x$results)
)

all_diagnostics <- bind_rows(
  lapply(production_objects, function(x) x$diagnostics)
)

failures <- all_results %>%
  filter(
    is.na(status) |
      status != "ok" |
      !is.finite(estimate) |
      !is.finite(se)
  )

# Save failed fits separately for review.
write_csv(
  failures,
  file.path(
    combined_dir_strong,
    "aipw_strong_failures.csv"
  )
)

# Save the combined estimator results.
saveRDS(
  all_results,
  file.path(
    combined_dir_strong,
    "aipw_strong_estimates.rds"
  ),
  compress = "xz"
)

write_csv(
  all_results,
  file.path(
    combined_dir_strong,
    "aipw_strong_estimates.csv"
  )
)

# Save the combined misspecification diagnostics.
saveRDS(
  all_diagnostics,
  file.path(
    combined_dir_strong,
    "aipw_strong_misspecification_diagnostics.rds"
  ),
  compress = "xz"
)

write_csv(
  all_diagnostics,
  file.path(
    combined_dir_strong,
    "aipw_strong_misspecification_diagnostics.csv"
  )
)

cat(
  "\nStrong AIPW sensitivity combination complete.\n",
  "Settings: ", n_settings, "\n",
  "Result rows: ", nrow(all_results), "\n",
  "Diagnostic rows: ", nrow(all_diagnostics), "\n",
  "Failures: ", nrow(failures), "\n",
  sep = ""
)

