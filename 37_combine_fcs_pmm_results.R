################################################################################
# 37_combine_fcs_pmm_results.R
#
# Validates and combines all FCS-PMM production chunks. Production repetitions
# 1-2 must exactly reproduce the completed FCS-PMM canary.
################################################################################

source("00_config.R")
source("34_fcs_pmm_sensitivity_helpers.R")

library(dplyr)
library(readr)


# ------------------------------------------------------------------------------
# 1. Settings, paths, and gates
# ------------------------------------------------------------------------------

nsim_fcs_pmm <- read_fcs_pmm_positive_integer_env(
  "SIM_FCS_PMM_NSIM",
  2000L
)

reps_per_fcs_pmm_chunk <- read_fcs_pmm_positive_integer_env(
  "SIM_FCS_PMM_REPS_PER_CHUNK",
  10L
)

if (nsim_fcs_pmm != 2000L || reps_per_fcs_pmm_chunk != 10L) {
  stop(
    "The prespecified FCS-PMM combination requires ",
    "SIM_FCS_PMM_NSIM=2000 and SIM_FCS_PMM_REPS_PER_CHUNK=10."
  )
}

chunks_per_fcs_pmm_scenario <- as.integer(
  ceiling(nsim_fcs_pmm / reps_per_fcs_pmm_chunk)
)

fcs_pmm_results_dir <- file.path(results_dir, "fcs_pmm_sensitivity")
fcs_pmm_canary_dir <- file.path(fcs_pmm_results_dir, "canary")
fcs_pmm_production_dir <- file.path(fcs_pmm_results_dir, "production")
fcs_pmm_chunk_dir <- file.path(fcs_pmm_production_dir, "chunks")
fcs_pmm_combined_dir <- file.path(fcs_pmm_production_dir, "combined")
fcs_pmm_scenario_dir <- file.path(fcs_pmm_combined_dir, "scenarios")
fcs_pmm_rng_dir <- file.path(fcs_pmm_combined_dir, "rng_seeds")
fcs_pmm_dataset_qc_dir <- file.path(fcs_pmm_combined_dir, "dataset_qc")

invisible(
  lapply(
    c(
      fcs_pmm_combined_dir,
      fcs_pmm_scenario_dir,
      fcs_pmm_rng_dir,
      fcs_pmm_dataset_qc_dir
    ),
    dir.create,
    showWarnings = FALSE,
    recursive = TRUE
  )
)

combine_pass_file <- file.path(
  fcs_pmm_combined_dir,
  "FCS_PMM_PRODUCTION_COMBINE_PASS.txt"
)

if (file.exists(combine_pass_file)) {
  unlink(combine_pass_file)
}

mcar_summary_pass_file <- file.path(
  results_dir,
  "phase2_mcar_sensitivity",
  "production",
  "summary",
  "PHASE2_MCAR_PRODUCTION_SUMMARY_PASS.txt"
)

canary_pass_file <- file.path(
  fcs_pmm_canary_dir,
  "FCS_PMM_CANARY_PASS.txt"
)

canary_qc_file <- file.path(
  fcs_pmm_canary_dir,
  "fcs_pmm_canary_qc.csv"
)

canary_defaults_file <- file.path(
  fcs_pmm_canary_dir,
  "fcs_pmm_default_settings.csv"
)

canary_results_file <- file.path(
  fcs_pmm_canary_dir,
  "fcs_pmm_canary_results.csv"
)

required_gate_files <- c(
  mcar_summary_pass_file,
  canary_pass_file,
  canary_qc_file,
  canary_defaults_file,
  canary_results_file
)

missing_gate_files <- required_gate_files[!file.exists(required_gate_files)]

if (length(missing_gate_files) > 0L) {
  stop(
    "Required FCS-PMM production gate files are missing:\n- ",
    paste(missing_gate_files, collapse = "\n- ")
  )
}

if (!any(grepl(
  "MCAR PRODUCTION SUMMARY: PASS",
  readLines(mcar_summary_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The Phase-2 MCAR production summary PASS marker is invalid.")
}

if (!any(grepl(
  "FCS-PMM CANARY: PASS",
  readLines(canary_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The FCS-PMM canary PASS marker is invalid.")
}

pmm_defaults <- fcs_pmm_default_settings()

canary_qc <- read_csv(canary_qc_file, show_col_types = FALSE)
canary_defaults <- read_csv(
  canary_defaults_file,
  show_col_types = FALSE
)

if (
  nrow(canary_qc) == 0L ||
    !fcs_pmm_all_true(canary_qc$passed) ||
    !fcs_pmm_all_true(canary_qc$overall_canary_pass)
) {
  stop("The FCS-PMM canary QC file does not record a full pass.")
}

if (!fcs_pmm_all_true(pmm_defaults$defaults_match_prespecified)) {
  stop("The current mice PMM defaults are not the prespecified defaults.")
}

if (!isTRUE(all.equal(
  pmm_defaults,
  canary_defaults,
  tolerance = 0,
  check.attributes = FALSE
))) {
  stop("The current mice version or PMM defaults differ from the canary.")
}


# ------------------------------------------------------------------------------
# 2. Combine and validate each scenario
# ------------------------------------------------------------------------------

manifest_list <- vector("list", nrow(fcs_pmm_sensitivity_grid))
combined_result_list <- vector("list", nrow(fcs_pmm_sensitivity_grid))
failure_list <- list()
failure_index <- 1L

for (s in seq_len(nrow(fcs_pmm_sensitivity_grid))) {

  row <- fcs_pmm_sensitivity_grid[s, , drop = FALSE]
  sensitivity_scenario_id <- row$sensitivity_scenario_id

  cat(
    "Combining FCS-PMM scenario ",
    sensitivity_scenario_id,
    " / ",
    nrow(fcs_pmm_sensitivity_grid),
    ": ",
    row$scenario_label,
    "\n",
    sep = ""
  )

  scenario_chunk_dir <- file.path(
    fcs_pmm_chunk_dir,
    sprintf(
      "scenario_%02d_%s",
      sensitivity_scenario_id,
      row$scenario_key
    )
  )

  expected_files <- file.path(
    scenario_chunk_dir,
    sprintf(
      "fcs_pmm_scenario_%02d_chunk_%03d.rds",
      sensitivity_scenario_id,
      seq_len(chunks_per_fcs_pmm_scenario)
    )
  )

  missing_files <- expected_files[!file.exists(expected_files)]

  if (length(missing_files) > 0L) {

    manifest_list[[s]] <- bind_cols(
      row,
      data.frame(
        chunks_expected = chunks_per_fcs_pmm_scenario,
        chunks_found =
          chunks_per_fcs_pmm_scenario - length(missing_files),
        repetitions_found = NA_integer_,
        result_rows = NA_integer_,
        complete = FALSE,
        problem = paste0("Missing ", length(missing_files), " chunk file(s)"),
        stringsAsFactors = FALSE
      )
    )

    next
  }

  chunk_objects <- lapply(
    expected_files,
    function(path) {
      tryCatch(readRDS(path), error = function(e) NULL)
    }
  )

  unreadable_chunks <- vapply(chunk_objects, is.null, logical(1))

  required_components <- c(
    "metadata",
    "results",
    "rng_states",
    "dataset_qc"
  )

  component_checks <- vapply(
    chunk_objects,
    function(x) {
      !is.null(x) && all(required_components %in% names(x))
    },
    logical(1)
  )

  if (any(unreadable_chunks) || !fcs_pmm_all_true(component_checks)) {
    manifest_list[[s]] <- bind_cols(
      row,
      data.frame(
        chunks_expected = chunks_per_fcs_pmm_scenario,
        chunks_found = chunks_per_fcs_pmm_scenario,
        repetitions_found = NA_integer_,
        result_rows = NA_integer_,
        complete = FALSE,
        problem = "Unreadable chunk or missing chunk-object component",
        stringsAsFactors = FALSE
      )
    )

    next
  }

  results <- bind_rows(lapply(chunk_objects, function(x) x$results))
  rng_seeds <- bind_rows(lapply(chunk_objects, function(x) x$rng_states))
  dataset_qc <- bind_rows(lapply(chunk_objects, function(x) x$dataset_qc))

  required_result_columns <- c(
    "sensitivity_scenario_id",
    "primary_scenario_id",
    "repetition",
    "method",
    "status",
    "estimate",
    "se",
    "mice_method",
    "nimp",
    "mice_maxit",
    "mice_version",
    "selection",
    "marker_error"
  )

  required_dataset_columns <- c(
    "sensitivity_scenario_id",
    "repetition",
    "marker_blockwise",
    "phase2_matches_markers",
    "finite_valid_probabilities",
    "sufficient_phase2_n"
  )

  required_rng_columns <- c(
    "sensitivity_scenario_id",
    "repetition",
    "data_rng_state",
    "method_rng_state"
  )

  missing_bound_columns <- c(
    setdiff(required_result_columns, names(results)),
    setdiff(required_dataset_columns, names(dataset_qc)),
    setdiff(required_rng_columns, names(rng_seeds))
  )

  if (length(missing_bound_columns) > 0L) {
    manifest_list[[s]] <- bind_cols(
      row,
      data.frame(
        chunks_expected = chunks_per_fcs_pmm_scenario,
        chunks_found = chunks_per_fcs_pmm_scenario,
        repetitions_found = NA_integer_,
        result_rows = nrow(results),
        complete = FALSE,
        problem = paste0(
          "Missing bound column(s): ",
          paste(unique(missing_bound_columns), collapse = ", ")
        ),
        stringsAsFactors = FALSE
      )
    )

    next
  }

  repetitions_found <- sort(unique(results$repetition))
  expected_repetitions <- seq_len(nsim_fcs_pmm)

  duplicate_keys <- results %>%
    count(sensitivity_scenario_id, repetition, method, name = "n_rows") %>%
    filter(n_rows != 1L)

  dataset_duplicate_keys <- dataset_qc %>%
    count(sensitivity_scenario_id, repetition, name = "n_rows") %>%
    filter(n_rows != 1L)

  rng_duplicate_keys <- rng_seeds %>%
    count(sensitivity_scenario_id, repetition, name = "n_rows") %>%
    filter(n_rows != 1L)

  metadata_matches <- vapply(
    seq_along(chunk_objects),
    function(j) {
      x <- chunk_objects[[j]]
      expected_rep_start <-
        (j - 1L) * reps_per_fcs_pmm_chunk + 1L
      expected_rep_end <- min(
        j * reps_per_fcs_pmm_chunk,
        nsim_fcs_pmm
      )

      identical(
        x$metadata$sensitivity_scenario_id,
        sensitivity_scenario_id
      ) &&
        identical(x$metadata$primary_scenario_id, row$primary_scenario_id) &&
        identical(x$metadata$chunk_id, j) &&
        identical(x$metadata$rep_start, expected_rep_start) &&
        identical(x$metadata$rep_end, expected_rep_end) &&
        identical(x$metadata$nsim_fcs_pmm, nsim_fcs_pmm) &&
        identical(
          x$metadata$reps_per_fcs_pmm_chunk,
          reps_per_fcs_pmm_chunk
        ) &&
        identical(x$metadata$selection, "mar") &&
        identical(x$metadata$marker_error, "mvn") &&
        identical(x$metadata$mice_method, "pmm") &&
        identical(x$metadata$nimp, nimp_primary) &&
        identical(x$metadata$mice_maxit, mice_maxit_primary) &&
        identical(x$metadata$pmm_donors, pmm_defaults$donors) &&
        identical(x$metadata$pmm_matchtype, pmm_defaults$matchtype) &&
        isTRUE(all.equal(
          x$metadata$pmm_ridge,
          pmm_defaults$ridge,
          tolerance = 0
        )) &&
        identical(x$metadata$mice_version, pmm_defaults$mice_version)
    },
    logical(1)
  )

  complete <-
    setequal(repetitions_found, expected_repetitions) &&
    setequal(unique(dataset_qc$repetition), expected_repetitions) &&
    setequal(unique(rng_seeds$repetition), expected_repetitions) &&
    nrow(results) == nsim_fcs_pmm &&
    nrow(dataset_qc) == nsim_fcs_pmm &&
    nrow(rng_seeds) == nsim_fcs_pmm &&
    nrow(duplicate_keys) == 0L &&
    nrow(dataset_duplicate_keys) == 0L &&
    nrow(rng_duplicate_keys) == 0L &&
    identical(unique(results$method), fcs_pmm_label) &&
    fcs_pmm_all_true(results$mice_method == "pmm") &&
    fcs_pmm_all_true(results$nimp == nimp_primary) &&
    fcs_pmm_all_true(results$mice_maxit == mice_maxit_primary) &&
    fcs_pmm_all_true(results$mice_version == pmm_defaults$mice_version) &&
    fcs_pmm_all_true(results$selection == "mar") &&
    fcs_pmm_all_true(results$marker_error == "mvn") &&
    fcs_pmm_all_true(dataset_qc$marker_blockwise) &&
    fcs_pmm_all_true(dataset_qc$phase2_matches_markers) &&
    fcs_pmm_all_true(dataset_qc$finite_valid_probabilities) &&
    fcs_pmm_all_true(dataset_qc$sufficient_phase2_n) &&
    fcs_pmm_all_true(metadata_matches)

  problems <- character()

  if (!setequal(repetitions_found, expected_repetitions)) {
    problems <- c(problems, "repetition set incomplete")
  }

  if (
    !setequal(unique(dataset_qc$repetition), expected_repetitions) ||
      !setequal(unique(rng_seeds$repetition), expected_repetitions)
  ) {
    problems <- c(problems, "dataset-QC or RNG repetition set incomplete")
  }

  if (nrow(results) != nsim_fcs_pmm) {
    problems <- c(
      problems,
      paste0("expected ", nsim_fcs_pmm, " rows, found ", nrow(results))
    )
  }

  if (
    nrow(dataset_qc) != nsim_fcs_pmm ||
      nrow(rng_seeds) != nsim_fcs_pmm
  ) {
    problems <- c(problems, "unexpected dataset-QC or RNG row count")
  }

  if (
    nrow(duplicate_keys) > 0L ||
      nrow(dataset_duplicate_keys) > 0L ||
      nrow(rng_duplicate_keys) > 0L
  ) {
    problems <- c(problems, "duplicate result, dataset-QC, or RNG keys")
  }

  if (!identical(unique(results$method), fcs_pmm_label)) {
    problems <- c(problems, "unexpected method label")
  }

  settings_match <-
    fcs_pmm_all_true(results$mice_method == "pmm") &&
    fcs_pmm_all_true(results$nimp == nimp_primary) &&
    fcs_pmm_all_true(results$mice_maxit == mice_maxit_primary) &&
    fcs_pmm_all_true(results$mice_version == pmm_defaults$mice_version) &&
    fcs_pmm_all_true(results$selection == "mar") &&
    fcs_pmm_all_true(results$marker_error == "mvn")

  if (!settings_match) {
    problems <- c(problems, "analysis or DGM settings mismatch")
  }

  if (!fcs_pmm_all_true(metadata_matches)) {
    problems <- c(problems, "chunk metadata mismatch")
  }

  if (
    !fcs_pmm_all_true(dataset_qc$marker_blockwise) ||
      !fcs_pmm_all_true(dataset_qc$phase2_matches_markers) ||
      !fcs_pmm_all_true(dataset_qc$finite_valid_probabilities) ||
      !fcs_pmm_all_true(dataset_qc$sufficient_phase2_n)
  ) {
    problems <- c(problems, "dataset structural QC failure")
  }

  if (complete) {

    fcs_pmm_atomic_save_rds(
      results,
      file.path(
        fcs_pmm_scenario_dir,
        sprintf(
          "fcs_pmm_scenario_%02d_estimates.rds",
          sensitivity_scenario_id
        )
      ),
      compress = "xz"
    )

    fcs_pmm_atomic_save_rds(
      rng_seeds,
      file.path(
        fcs_pmm_rng_dir,
        sprintf(
          "fcs_pmm_scenario_%02d_rng_seeds.rds",
          sensitivity_scenario_id
        )
      ),
      compress = "xz"
    )

    fcs_pmm_atomic_save_rds(
      dataset_qc,
      file.path(
        fcs_pmm_dataset_qc_dir,
        sprintf(
          "fcs_pmm_scenario_%02d_dataset_qc.rds",
          sensitivity_scenario_id
        )
      ),
      compress = "xz"
    )
  }

  failed_rows <- results %>%
    filter(
      is.na(status) |
        status != "ok" |
        !is.finite(estimate) |
        !is.finite(se) |
        se <= 0
    )

  if (nrow(failed_rows) > 0L) {
    failure_list[[failure_index]] <- failed_rows
    failure_index <- failure_index + 1L
  }

  combined_result_list[[s]] <- results

  manifest_list[[s]] <- bind_cols(
    row,
    data.frame(
      chunks_expected = chunks_per_fcs_pmm_scenario,
      chunks_found = chunks_per_fcs_pmm_scenario,
      repetitions_found = length(repetitions_found),
      result_rows = nrow(results),
      complete = complete,
      problem = if (length(problems) == 0L) {
        ""
      } else {
        paste(problems, collapse = "; ")
      },
      stringsAsFactors = FALSE
    )
  )
}

manifest <- bind_rows(manifest_list)

fcs_pmm_atomic_write_csv(
  manifest,
  file.path(fcs_pmm_combined_dir, "fcs_pmm_completion_manifest.csv")
)

failures <- if (length(failure_list) == 0L) {
  tibble::tibble(
    sensitivity_scenario_id = integer(),
    primary_scenario_id = integer(),
    repetition = integer(),
    method = character(),
    status = character(),
    message = character(),
    estimate = numeric(),
    se = numeric()
  )
} else {
  bind_rows(failure_list)
}

fcs_pmm_atomic_write_csv(
  failures,
  file.path(fcs_pmm_combined_dir, "fcs_pmm_method_failures.csv")
)

if (
  nrow(manifest) != nrow(fcs_pmm_sensitivity_grid) ||
    !fcs_pmm_all_true(manifest$complete)
) {
  cat(
    "\nFCS-PMM production is incomplete. See:\n",
    file.path(fcs_pmm_combined_dir, "fcs_pmm_completion_manifest.csv"),
    "\n",
    sep = ""
  )
  quit(save = "no", status = 2L)
}


# ------------------------------------------------------------------------------
# 3. Reproduce FCS-PMM canary results from production repetitions 1-2
# ------------------------------------------------------------------------------

canary_pmm <- read_csv(canary_results_file, show_col_types = FALSE) %>%
  filter(method == fcs_pmm_label) %>%
  arrange(primary_scenario_id, repetition)

production_pmm <- bind_rows(combined_result_list) %>%
  filter(
    primary_scenario_id %in% fcs_pmm_canary_grid$primary_scenario_id,
    repetition %in% 1:2
  ) %>%
  arrange(primary_scenario_id, repetition)

excluded_columns <- c(
  "elapsed_seconds",
  "task_id",
  "chunk_id",
  "message"
)

comparison_columns <- setdiff(
  intersect(names(canary_pmm), names(production_pmm)),
  excluded_columns
)

comparison_key <- c("primary_scenario_id", "repetition", "method")
comparison_columns <- unique(c(comparison_key, comparison_columns))

canary_reproduction_list <- vector(
  "list",
  nrow(fcs_pmm_canary_grid)
)

for (i in seq_len(nrow(fcs_pmm_canary_grid))) {

  primary_id <- fcs_pmm_canary_grid$primary_scenario_id[[i]]

  canary_subset <- canary_pmm %>%
    filter(primary_scenario_id == primary_id) %>%
    select(all_of(comparison_columns)) %>%
    arrange(repetition)

  production_subset <- production_pmm %>%
    filter(primary_scenario_id == primary_id) %>%
    select(all_of(comparison_columns)) %>%
    arrange(repetition)

  rownames(canary_subset) <- NULL
  rownames(production_subset) <- NULL

  rows_complete <-
    nrow(canary_subset) == 2L &&
    nrow(production_subset) == 2L

  # Canonicalize canary/production comparison-key storage types before
  # requiring exact equality.
  #
  # Canary results are read from CSV, so integer-valued identifiers may be
  # imported as doubles. Production results are read from RDS and retain
  # integer storage. Convert both sides to the same explicit representation
  # before using identical(); the scientific key values are unchanged.
  canary_keys_compare <- canary_subset %>%
    transmute(
      primary_scenario_id = as.integer(primary_scenario_id),
      repetition = as.integer(repetition),
      method = as.character(method)
    ) %>%
    as.data.frame()

  production_keys_compare <- production_subset %>%
    transmute(
      primary_scenario_id = as.integer(primary_scenario_id),
      repetition = as.integer(repetition),
      method = as.character(method)
    ) %>%
    as.data.frame()

  keys_complete <-
    rows_complete &&
    !any(duplicated(canary_keys_compare)) &&
    !any(duplicated(production_keys_compare)) &&
    identical(
      canary_keys_compare,
      production_keys_compare
    )

  comparison <- if (keys_complete) {
    all.equal(
      production_subset,
      canary_subset,
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  } else {
    "Expected row counts or comparison keys do not match."
  }

  numeric_columns <- comparison_columns[
    vapply(canary_subset, is.numeric, logical(1))
  ]

  numeric_differences <- if (keys_complete) unlist(
    lapply(
      numeric_columns,
      function(column_name) {
        difference <- abs(
          production_subset[[column_name]] - canary_subset[[column_name]]
        )
        difference[is.finite(difference)]
      }
    ),
    use.names = FALSE
  ) else numeric()

  reproduced <- rows_complete && keys_complete && isTRUE(comparison)

  canary_reproduction_list[[i]] <- data.frame(
    primary_scenario_id = primary_id,
    expected_rows = 2L,
    canary_rows = nrow(canary_subset),
    production_rows = nrow(production_subset),
    max_abs_numeric_difference = if (
      length(numeric_differences) == 0L
    ) {
      NA_real_
    } else {
      max(numeric_differences)
    },
    reproduced_within_tolerance = reproduced,
    comparison_message = if (reproduced) {
      ""
    } else {
      paste(comparison, collapse = "; ")
    },
    stringsAsFactors = FALSE
  )
}

canary_reproduction_qc <- bind_rows(canary_reproduction_list)

fcs_pmm_atomic_write_csv(
  canary_reproduction_qc,
  file.path(fcs_pmm_combined_dir, "fcs_pmm_canary_reproduction_qc.csv")
)


# ------------------------------------------------------------------------------
# 4. Final production-content gate
# ------------------------------------------------------------------------------

production_qc <- data.frame(
  check = c(
    "All four FCS-PMM production scenarios are complete",
    "Each scenario contains exactly 2,000 repetitions",
    "Each scenario contains exactly 2,000 FCS-PMM result rows",
    "All chunk and dataset structural checks pass",
    "Production repetitions reproduce the completed FCS-PMM canary"
  ),
  passed = c(
    nrow(manifest) == 4L && fcs_pmm_all_true(manifest$complete),
    fcs_pmm_all_true(manifest$repetitions_found == nsim_fcs_pmm),
    fcs_pmm_all_true(manifest$result_rows == nsim_fcs_pmm),
    fcs_pmm_all_true(manifest$problem == ""),
    nrow(canary_reproduction_qc) == 2L &&
      fcs_pmm_all_true(canary_reproduction_qc$canary_rows == 2L) &&
      fcs_pmm_all_true(canary_reproduction_qc$production_rows == 2L) &&
      fcs_pmm_all_true(canary_reproduction_qc$reproduced_within_tolerance)
  ),
  stringsAsFactors = FALSE
)

production_qc$overall_production_content_pass <-
  fcs_pmm_all_true(production_qc$passed)

fcs_pmm_atomic_write_csv(
  production_qc,
  file.path(fcs_pmm_combined_dir, "fcs_pmm_production_qc.csv")
)

if (!fcs_pmm_all_true(production_qc$passed)) {

  failed_checks <- production_qc$check[!production_qc$passed]

  stop(
    "FCS-PMM production combination FAILED. Failed checks:\n- ",
    paste(failed_checks, collapse = "\n- "),
    "\nSee outputs in: ",
    fcs_pmm_combined_dir
  )
}

fcs_pmm_atomic_write_lines(
  c(
    "FCS-PMM PRODUCTION COMBINATION: PASS",
    paste("Completed:", Sys.time()),
    paste("Production scenarios:", nrow(manifest)),
    paste("Repetitions per scenario:", nsim_fcs_pmm),
    paste("Method-level failures recorded:", nrow(failures)),
    "Canary reproduction: PASS"
  ),
  combine_pass_file
)

cat(
  "\nFCS-PMM production combination PASSED.\n",
  "Combined directory: ",
  fcs_pmm_combined_dir,
  "\n",
  sep = ""
)
