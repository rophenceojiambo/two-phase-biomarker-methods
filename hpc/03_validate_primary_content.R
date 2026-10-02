################################################################################
# hpc/03_validate_primary_content.R
#
# Read-only content-level QC for the completed MIDUS two-phase primary
# simulation. This script DOES NOT combine, overwrite, or rerun simulations.
#
# It checks all 6,480 final chunk files for:
#   - readability and expected object structure
#   - correct task/scenario/chunk metadata
#   - exactly 100 repetitions per chunk
#   - exactly six methods per repetition, once each
#   - correct scenario parameters
#   - valid RNG-state records
#   - finite/valid numerical outputs among rows marked status == "ok"
#   - method-level failures/non-finite estimates (reported separately)
#
# Method failures are NOT treated as structural corruption because failure rate
# is itself a simulation performance measure. They are counted and reported.
################################################################################

source("00_config.R")

library(dplyr)
library(readr)

qc_dir <- file.path(results_dir, "qc")
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

expected_methods <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

required_result_columns <- c(
  "task_id",
  "scenario_id",
  "chunk_id",
  "repetition",
  "method",
  "estimate",
  "se",
  "df",
  "p_value",
  "conf_low",
  "conf_high",
  "status",
  "message",
  "elapsed_seconds",
  "N",
  "target_phase2_fraction",
  "realized_phase2_fraction",
  "n_phase2",
  "r2_a_marker",
  "r2_y_marker",
  "theta_true",
  "min_pi_true",
  "p01_pi_true",
  "median_pi_true",
  "p99_pi_true",
  "max_pi_true"
)

required_rng_columns <- c(
  "scenario_id",
  "repetition",
  "data_rng_state",
  "method_rng_state"
)

same_num <- function(x, expected, tol = 1e-12) {
  length(x) > 0L &&
    all(!is.na(x)) &&
    all(abs(as.numeric(x) - as.numeric(expected)) <= tol)
}

same_int <- function(x, expected) {
  length(x) > 0L &&
    all(!is.na(x)) &&
    all(as.integer(x) == as.integer(expected))
}

same_chr <- function(x, expected) {
  length(x) > 0L &&
    all(!is.na(x)) &&
    all(as.character(x) == as.character(expected))
}

append_problem <- function(problems, condition, text) {
  if (!isTRUE(condition)) c(problems, text) else problems
}

qc_rows <- vector("list", n_primary_array_tasks)

method_summary <- data.frame(
  method = expected_methods,
  rows = integer(length(expected_methods)),
  status_failed = integer(length(expected_methods)),
  nonfinite_estimate_or_se = integer(length(expected_methods)),
  invalid_ok_rows = integer(length(expected_methods)),
  stringsAsFactors = FALSE
)

issue_sample <- list()
issue_sample_n <- 0L
issue_sample_limit <- 500L

cat(
  "Starting content-level QC of ",
  n_primary_array_tasks,
  " chunk files...\n",
  sep = ""
)

for (task_id in seq_len(n_primary_array_tasks)) {

  scenario_id <- ((task_id - 1L) %/% chunks_per_scenario) + 1L
  chunk_id <- ((task_id - 1L) %% chunks_per_scenario) + 1L

  rep_start <- ((chunk_id - 1L) * reps_per_chunk) + 1L
  rep_end <- min(chunk_id * reps_per_chunk, nsim_primary)
  expected_reps <- seq.int(rep_start, rep_end)
  expected_nrep <- length(expected_reps)
  expected_rows <- expected_nrep * length(expected_methods)

  scenario_row <- primary_grid[
    primary_grid$scenario_id == scenario_id,
    ,
    drop = FALSE
  ]

  chunk_file <- file.path(
    chunk_dir,
    sprintf("scenario_%03d", scenario_id),
    sprintf("scenario_%03d_chunk_%02d.rds", scenario_id, chunk_id)
  )

  problems <- character()

  file_exists <- file.exists(chunk_file)
  read_ok <- FALSE
  object_structure_ok <- FALSE
  metadata_ok <- FALSE
  results_columns_ok <- FALSE
  result_rows_ok <- FALSE
  repetition_set_ok <- FALSE
  method_set_ok <- FALSE
  repetition_method_keys_ok <- FALSE
  identifiers_ok <- FALSE
  scenario_parameters_ok <- FALSE
  rng_columns_ok <- FALSE
  rng_rows_ok <- FALSE
  rng_repetition_set_ok <- FALSE
  rng_identifiers_ok <- FALSE
  rng_states_ok <- FALSE
  numerical_ok_rows_valid <- FALSE

  result_rows <- NA_integer_
  rng_rows <- NA_integer_
  method_issue_rows <- NA_integer_
  invalid_ok_rows <- NA_integer_

  if (!file_exists) {
    problems <- c(problems, "final chunk file missing")
  } else {

    obj <- tryCatch(
      readRDS(chunk_file),
      error = function(e) e
    )

    read_ok <- !inherits(obj, "error")

    if (!read_ok) {
      problems <- c(
        problems,
        paste0("readRDS failed: ", conditionMessage(obj))
      )
    } else {

      object_structure_ok <- (
        is.list(obj) &&
          all(c("metadata", "results", "rng_states") %in% names(obj)) &&
          is.list(obj$metadata) &&
          is.data.frame(obj$results) &&
          is.data.frame(obj$rng_states)
      )

      problems <- append_problem(
        problems,
        object_structure_ok,
        "object does not contain metadata/results/rng_states in expected form"
      )

      if (object_structure_ok) {

        md <- obj$metadata
        results <- obj$results
        rng_states <- obj$rng_states

        result_rows <- nrow(results)
        rng_rows <- nrow(rng_states)

        # ----------------------------------------------------------------------
        # Metadata checks
        # ----------------------------------------------------------------------

        metadata_ok <- (
          same_int(md$task_id, task_id) &&
            same_int(md$scenario_id, scenario_id) &&
            same_int(md$chunk_id, chunk_id) &&
            same_int(md$rep_start, rep_start) &&
            same_int(md$rep_end, rep_end) &&
            same_int(md$N, scenario_row$N) &&
            same_num(md$phase2_fraction, scenario_row$phase2_fraction) &&
            same_num(md$r2_a_marker, scenario_row$r2_a_marker) &&
            same_num(md$r2_y_marker, scenario_row$r2_y_marker) &&
            same_num(md$theta, scenario_row$theta) &&
            same_int(md$nimp, nimp_primary) &&
            same_chr(md$mice_method, mice_method_primary) &&
            same_int(md$mice_maxit, mice_maxit_primary) &&
            same_int(md$jomo_nburn, jomo_nburn_primary) &&
            same_int(md$jomo_nbetween, jomo_nbetween_primary)
        )

        problems <- append_problem(
          problems,
          metadata_ok,
          "metadata does not match expected production settings"
        )

        # ----------------------------------------------------------------------
        # Result-table structural checks
        # ----------------------------------------------------------------------

        results_columns_ok <- all(
          required_result_columns %in% names(results)
        )

        problems <- append_problem(
          problems,
          results_columns_ok,
          paste0(
            "missing required result columns: ",
            paste(
              setdiff(required_result_columns, names(results)),
              collapse = ", "
            )
          )
        )

        if (results_columns_ok) {

          result_rows_ok <- nrow(results) == expected_rows

          repetition_set_ok <- (
            length(unique(results$repetition)) == expected_nrep &&
              setequal(as.integer(results$repetition), expected_reps)
          )

          method_set_ok <- setequal(
            unique(as.character(results$method)),
            expected_methods
          )

          keys <- paste(
            results$repetition,
            results$method,
            sep = "|"
          )

          method_counts <- table(
            factor(
              as.character(results$method),
              levels = expected_methods
            )
          )

          rep_counts <- table(results$repetition)

          repetition_method_keys_ok <- (
            length(keys) == expected_rows &&
              anyDuplicated(keys) == 0L &&
              all(as.integer(method_counts) == expected_nrep) &&
              length(rep_counts) == expected_nrep &&
              all(as.integer(rep_counts) == length(expected_methods))
          )

          identifiers_ok <- (
            same_int(results$task_id, task_id) &&
              same_int(results$scenario_id, scenario_id) &&
              same_int(results$chunk_id, chunk_id)
          )

          scenario_parameters_ok <- (
            same_int(results$N, scenario_row$N) &&
              same_num(
                results$target_phase2_fraction,
                scenario_row$phase2_fraction
              ) &&
              same_num(results$r2_a_marker, scenario_row$r2_a_marker) &&
              same_num(results$r2_y_marker, scenario_row$r2_y_marker) &&
              same_num(results$theta_true, scenario_row$theta)
          )

          problems <- append_problem(
            problems,
            result_rows_ok,
            paste0(
              "expected ", expected_rows,
              " result rows, found ", nrow(results)
            )
          )

          problems <- append_problem(
            problems,
            repetition_set_ok,
            "result repetition set is incomplete or incorrect"
          )

          problems <- append_problem(
            problems,
            method_set_ok,
            "result method set differs from the six expected methods"
          )

          problems <- append_problem(
            problems,
            repetition_method_keys_ok,
            "duplicate/missing repetition-method combinations"
          )

          problems <- append_problem(
            problems,
            identifiers_ok,
            "task/scenario/chunk identifiers in results are inconsistent"
          )

          problems <- append_problem(
            problems,
            scenario_parameters_ok,
            "scenario parameters in results do not match primary_grid"
          )

          # --------------------------------------------------------------------
          # Method-level and numerical-output checks
          # --------------------------------------------------------------------

          bad_status <- is.na(results$status) |
            as.character(results$status) != "ok"

          nonfinite_est_se <- (
            !is.finite(results$estimate) |
              !is.finite(results$se)
          )

          ok_rows <- !bad_status

          invalid_ok <- ok_rows & (
            !is.finite(results$estimate) |
              !is.finite(results$se) |
              results$se <= 0 |
              !is.finite(results$p_value) |
              results$p_value < 0 |
              results$p_value > 1 |
              !is.finite(results$conf_low) |
              !is.finite(results$conf_high) |
              results$conf_low > results$conf_high |
              !is.finite(results$elapsed_seconds) |
              results$elapsed_seconds < 0 |
              !is.finite(results$realized_phase2_fraction) |
              results$realized_phase2_fraction < 0 |
              results$realized_phase2_fraction > 1 |
              !is.finite(results$n_phase2) |
              results$n_phase2 < 0 |
              results$n_phase2 > results$N |
              abs(
                results$realized_phase2_fraction -
                  results$n_phase2 / results$N
              ) > 1e-12 |
              !is.finite(results$min_pi_true) |
              !is.finite(results$p01_pi_true) |
              !is.finite(results$median_pi_true) |
              !is.finite(results$p99_pi_true) |
              !is.finite(results$max_pi_true) |
              results$min_pi_true < 0 |
              results$max_pi_true > 1 |
              results$min_pi_true > results$p01_pi_true |
              results$p01_pi_true > results$median_pi_true |
              results$median_pi_true > results$p99_pi_true |
              results$p99_pi_true > results$max_pi_true
          )

          method_issue <- bad_status | nonfinite_est_se

          method_issue_rows <- sum(method_issue)
          invalid_ok_rows <- sum(invalid_ok)
          numerical_ok_rows_valid <- invalid_ok_rows == 0L

          problems <- append_problem(
            problems,
            numerical_ok_rows_valid,
            paste0(
              invalid_ok_rows,
              " row(s) marked status=ok have invalid numerical values"
            )
          )

          # Aggregate method-specific counts.
          for (m in seq_along(expected_methods)) {
            idx <- as.character(results$method) == expected_methods[m]

            method_summary$rows[m] <- method_summary$rows[m] + sum(idx)
            method_summary$status_failed[m] <-
              method_summary$status_failed[m] + sum(idx & bad_status)
            method_summary$nonfinite_estimate_or_se[m] <-
              method_summary$nonfinite_estimate_or_se[m] +
              sum(idx & nonfinite_est_se)
            method_summary$invalid_ok_rows[m] <-
              method_summary$invalid_ok_rows[m] + sum(idx & invalid_ok)
          }

          # Save a bounded sample of problematic method rows for inspection.
          if (any(method_issue | invalid_ok) &&
              issue_sample_n < issue_sample_limit) {

            issue_idx <- which(method_issue | invalid_ok)
            keep_n <- min(
              length(issue_idx),
              issue_sample_limit - issue_sample_n
            )

            keep <- issue_idx[seq_len(keep_n)]

            sample_cols <- intersect(
              c(
                "task_id",
                "scenario_id",
                "chunk_id",
                "repetition",
                "method",
                "estimate",
                "se",
                "p_value",
                "conf_low",
                "conf_high",
                "status",
                "message",
                "N",
                "target_phase2_fraction",
                "r2_a_marker",
                "r2_y_marker",
                "theta_true"
              ),
              names(results)
            )

            tmp <- results[keep, sample_cols, drop = FALSE]
            tmp$qc_invalid_ok_row <- invalid_ok[keep]

            issue_sample[[length(issue_sample) + 1L]] <- tmp
            issue_sample_n <- issue_sample_n + nrow(tmp)
          }
        }

        # ----------------------------------------------------------------------
        # RNG-state checks
        # ----------------------------------------------------------------------

        rng_columns_ok <- all(
          required_rng_columns %in% names(rng_states)
        )

        problems <- append_problem(
          problems,
          rng_columns_ok,
          paste0(
            "missing required RNG columns: ",
            paste(
              setdiff(required_rng_columns, names(rng_states)),
              collapse = ", "
            )
          )
        )

        if (rng_columns_ok) {

          rng_rows_ok <- nrow(rng_states) == expected_nrep

          rng_repetition_set_ok <- (
            length(unique(rng_states$repetition)) == expected_nrep &&
              setequal(as.integer(rng_states$repetition), expected_reps)
          )

          rng_identifiers_ok <- same_int(
            rng_states$scenario_id,
            scenario_id
          )

          data_state <- as.character(rng_states$data_rng_state)
          method_state <- as.character(rng_states$method_rng_state)

          rng_states_ok <- (
            all(!is.na(data_state)) &&
              all(!is.na(method_state)) &&
              all(nzchar(data_state)) &&
              all(nzchar(method_state)) &&
              anyDuplicated(data_state) == 0L &&
              anyDuplicated(method_state) == 0L &&
              all(data_state != method_state)
          )

          problems <- append_problem(
            problems,
            rng_rows_ok,
            paste0(
              "expected ", expected_nrep,
              " RNG rows, found ", nrow(rng_states)
            )
          )

          problems <- append_problem(
            problems,
            rng_repetition_set_ok,
            "RNG repetition set is incomplete or incorrect"
          )

          problems <- append_problem(
            problems,
            rng_identifiers_ok,
            "RNG scenario identifiers are inconsistent"
          )

          problems <- append_problem(
            problems,
            rng_states_ok,
            "RNG state strings are missing, duplicated within chunk, or paired incorrectly"
          )
        }
      }
    }
  }

  structural_ok <- (
    file_exists &&
      read_ok &&
      object_structure_ok &&
      metadata_ok &&
      results_columns_ok &&
      result_rows_ok &&
      repetition_set_ok &&
      method_set_ok &&
      repetition_method_keys_ok &&
      identifiers_ok &&
      scenario_parameters_ok &&
      rng_columns_ok &&
      rng_rows_ok &&
      rng_repetition_set_ok &&
      rng_identifiers_ok &&
      rng_states_ok
  )

  qc_rows[[task_id]] <- data.frame(
    task_id = task_id,
    scenario_id = scenario_id,
    chunk_id = chunk_id,
    rep_start = rep_start,
    rep_end = rep_end,
    file = chunk_file,
    file_exists = file_exists,
    read_ok = read_ok,
    object_structure_ok = object_structure_ok,
    metadata_ok = metadata_ok,
    results_columns_ok = results_columns_ok,
    result_rows = result_rows,
    result_rows_ok = result_rows_ok,
    repetition_set_ok = repetition_set_ok,
    method_set_ok = method_set_ok,
    repetition_method_keys_ok = repetition_method_keys_ok,
    identifiers_ok = identifiers_ok,
    scenario_parameters_ok = scenario_parameters_ok,
    rng_rows = rng_rows,
    rng_columns_ok = rng_columns_ok,
    rng_rows_ok = rng_rows_ok,
    rng_repetition_set_ok = rng_repetition_set_ok,
    rng_identifiers_ok = rng_identifiers_ok,
    rng_states_ok = rng_states_ok,
    method_issue_rows = method_issue_rows,
    invalid_ok_rows = invalid_ok_rows,
    numerical_ok_rows_valid = numerical_ok_rows_valid,
    structural_ok = structural_ok,
    problem = paste(unique(problems), collapse = "; "),
    stringsAsFactors = FALSE
  )

  if (task_id %% 250L == 0L || task_id == n_primary_array_tasks) {
    cat(
      "Checked ",
      task_id,
      " / ",
      n_primary_array_tasks,
      " files.\n",
      sep = ""
    )
  }
}

qc_manifest <- bind_rows(qc_rows)

# ------------------------------------------------------------------------------
# Scenario-level summaries
# ------------------------------------------------------------------------------

scenario_qc <- qc_manifest %>%
  group_by(scenario_id) %>%
  summarise(
    chunks_found = sum(file_exists),
    chunks_structurally_ok = sum(structural_ok),
    result_rows = sum(result_rows, na.rm = TRUE),
    rng_rows = sum(rng_rows, na.rm = TRUE),
    method_issue_rows = sum(method_issue_rows, na.rm = TRUE),
    invalid_ok_rows = sum(invalid_ok_rows, na.rm = TRUE),
    scenario_structural_ok = (
      n() == chunks_per_scenario &&
        chunks_found == chunks_per_scenario &&
        chunks_structurally_ok == chunks_per_scenario &&
        result_rows ==
          nsim_primary * length(expected_methods) &&
        rng_rows == nsim_primary
    ),
    .groups = "drop"
  )

method_summary <- method_summary %>%
  mutate(
    expected_rows = nsim_primary * nrow(primary_grid),
    failure_rate = status_failed / rows
  )

issues_sample_df <- if (length(issue_sample) == 0L) {
  tibble::tibble()
} else {
  bind_rows(issue_sample)
}

# ------------------------------------------------------------------------------
# Global checks
# ------------------------------------------------------------------------------

expected_total_result_rows <-
  nrow(primary_grid) *
  nsim_primary *
  length(expected_methods)

expected_total_rng_rows <-
  nrow(primary_grid) *
  nsim_primary

actual_total_result_rows <- sum(
  qc_manifest$result_rows,
  na.rm = TRUE
)

actual_total_rng_rows <- sum(
  qc_manifest$rng_rows,
  na.rm = TRUE
)

n_structural_bad <- sum(!qc_manifest$structural_ok)
n_scenario_bad <- sum(!scenario_qc$scenario_structural_ok)
n_invalid_ok <- sum(qc_manifest$invalid_ok_rows, na.rm = TRUE)
n_method_issues <- sum(qc_manifest$method_issue_rows, na.rm = TRUE)

structural_pass <- (
  nrow(qc_manifest) == n_primary_array_tasks &&
    n_structural_bad == 0L &&
    n_scenario_bad == 0L &&
    actual_total_result_rows == expected_total_result_rows &&
    actual_total_rng_rows == expected_total_rng_rows
)

numeric_pass <- n_invalid_ok == 0L

# ------------------------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------------------------

write_csv(
  qc_manifest,
  file.path(qc_dir, "primary_content_qc_chunk_manifest.csv")
)

write_csv(
  scenario_qc,
  file.path(qc_dir, "primary_content_qc_scenario_summary.csv")
)

write_csv(
  method_summary,
  file.path(qc_dir, "primary_content_qc_method_summary.csv")
)

write_csv(
  issues_sample_df,
  file.path(qc_dir, "primary_content_qc_issue_sample.csv")
)

summary_lines <- c(
  "MIDUS TWO-PHASE PRIMARY SIMULATION — CONTENT-LEVEL QC",
  paste0("QC completed: ", Sys.time()),
  "",
  paste0("Expected chunk files: ", n_primary_array_tasks),
  paste0("Chunk files checked: ", nrow(qc_manifest)),
  paste0("Structurally invalid chunks: ", n_structural_bad),
  paste0("Structurally invalid scenarios: ", n_scenario_bad),
  "",
  paste0(
    "Expected result rows: ",
    format(expected_total_result_rows, big.mark = ",")
  ),
  paste0(
    "Observed result rows: ",
    format(actual_total_result_rows, big.mark = ",")
  ),
  paste0(
    "Expected RNG rows: ",
    format(expected_total_rng_rows, big.mark = ",")
  ),
  paste0(
    "Observed RNG rows: ",
    format(actual_total_rng_rows, big.mark = ",")
  ),
  "",
  paste0(
    "Method-level issue rows (status failure or non-finite estimate/SE): ",
    format(n_method_issues, big.mark = ",")
  ),
  paste0(
    "Invalid numerical rows among status='ok': ",
    format(n_invalid_ok, big.mark = ",")
  ),
  "",
  paste0(
    "STRUCTURAL CONTENT QC: ",
    if (structural_pass) "PASS" else "FAIL"
  ),
  paste0(
    "NUMERICAL VALIDITY OF status='ok' ROWS: ",
    if (numeric_pass) "PASS" else "FAIL"
  ),
  "",
  "Note: method-level failures are reported but do not by themselves imply",
  "structural corruption because method failure rate is a simulation performance",
  "measure. Inspect primary_content_qc_method_summary.csv before combining."
)

writeLines(
  summary_lines,
  file.path(qc_dir, "primary_content_qc_summary.txt")
)

cat(
  "\n============================================================\n",
  paste(summary_lines, collapse = "\n"),
  "\n============================================================\n",
  sep = ""
)

if (!structural_pass || !numeric_pass) {
  quit(save = "no", status = 2L)
}

cat(
  "\nCONTENT-LEVEL QC PASSED.\n",
  "It is safe to proceed to 06_combine_primary_results.R, ",
  "after reviewing the method-level issue summary.\n",
  sep = ""
)
