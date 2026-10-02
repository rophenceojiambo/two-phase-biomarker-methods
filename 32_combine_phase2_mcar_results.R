################################################################################
# 32_combine_phase2_mcar_results.R
#
# Validates and combines all Phase-2 MCAR production chunks.
# The script preserves repetition-level estimates, RNG seeds, and dataset QC.
# It also verifies that production repetitions 1-2 exactly reproduce the
# completed canary for primary scenarios 168 and 208, within numerical tolerance.
################################################################################

source("00_config.R")
source("28_phase2_mcar_sensitivity_helpers.R")

library(dplyr)
library(readr)


# ------------------------------------------------------------------------------
# 1. Settings and paths
# ------------------------------------------------------------------------------

nsim_mcar <- read_mcar_positive_integer_env(
  "SIM_MCAR_NSIM",
  2000L
)

reps_per_mcar_chunk <- read_mcar_positive_integer_env(
  "SIM_MCAR_REPS_PER_CHUNK",
  10L
)

if (nsim_mcar != 2000L || reps_per_mcar_chunk != 10L) {
  stop(
    "The prespecified MCAR production combination requires ",
    "SIM_MCAR_NSIM=2000 and SIM_MCAR_REPS_PER_CHUNK=10."
  )
}

chunks_per_mcar_scenario <- as.integer(
  ceiling(nsim_mcar / reps_per_mcar_chunk)
)

mcar_results_dir <- file.path(
  results_dir,
  "phase2_mcar_sensitivity"
)

mcar_canary_dir <- file.path(
  mcar_results_dir,
  "canary"
)

mcar_production_dir <- file.path(
  mcar_results_dir,
  "production"
)

mcar_cache_dir <- file.path(
  mcar_production_dir,
  "cache"
)

mcar_chunk_dir <- file.path(
  mcar_production_dir,
  "chunks"
)

mcar_combined_dir <- file.path(
  mcar_production_dir,
  "combined"
)

mcar_scenario_dir <- file.path(
  mcar_combined_dir,
  "scenarios"
)

mcar_rng_dir <- file.path(
  mcar_combined_dir,
  "rng_seeds"
)

mcar_dataset_qc_dir <- file.path(
  mcar_combined_dir,
  "dataset_qc"
)

combine_pass_file <- file.path(
  mcar_combined_dir,
  "PHASE2_MCAR_PRODUCTION_COMBINE_PASS.txt"
)

invisible(
  lapply(
    c(
      mcar_combined_dir,
      mcar_scenario_dir,
      mcar_rng_dir,
      mcar_dataset_qc_dir
    ),
    dir.create,
    showWarnings = FALSE,
    recursive = TRUE
  )
)

if (file.exists(combine_pass_file)) {
  unlink(combine_pass_file)
}

cache_pass_file <- file.path(
  mcar_cache_dir,
  "PHASE2_MCAR_PRODUCTION_CACHE_PASS.txt"
)

cache_qc_file <- file.path(
  mcar_cache_dir,
  "phase2_mcar_production_cache_qc.csv"
)

canary_pass_file <- file.path(
  mcar_canary_dir,
  "PHASE2_MCAR_CANARY_PASS.txt"
)

canary_qc_file <- file.path(
  mcar_canary_dir,
  "phase2_mcar_canary_qc.csv"
)

canary_results_file <- file.path(
  mcar_canary_dir,
  "phase2_mcar_canary_results.csv"
)

required_gate_files <- c(
  cache_pass_file,
  cache_qc_file,
  canary_pass_file,
  canary_qc_file,
  canary_results_file
)

missing_gate_files <- required_gate_files[
  !file.exists(required_gate_files)
]

if (length(missing_gate_files) > 0L) {
  stop(
    "Required MCAR production gate files are missing:\n- ",
    paste(missing_gate_files, collapse = "\n- ")
  )
}

if (!any(grepl(
  "MCAR PRODUCTION CACHE: PASS",
  readLines(cache_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The Phase-2 MCAR production cache PASS marker is invalid.")
}

if (!any(grepl(
  "MCAR CANARY: PASS",
  readLines(canary_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The Phase-2 MCAR canary PASS marker is invalid.")
}

cache_qc <- read_csv(cache_qc_file, show_col_types = FALSE)
canary_qc <- read_csv(canary_qc_file, show_col_types = FALSE)

if (
  nrow(cache_qc) == 0L ||
    !mcar_all_true(cache_qc$passed) ||
    !mcar_all_true(cache_qc$overall_cache_pass)
) {
  stop("The MCAR production cache QC file does not record a full pass.")
}

if (
  nrow(canary_qc) == 0L ||
    !mcar_all_true(canary_qc$passed) ||
    !mcar_all_true(canary_qc$overall_canary_pass)
) {
  stop("The MCAR canary QC file does not record a full pass.")
}


# ------------------------------------------------------------------------------
# 2. Combine and validate each scenario
# ------------------------------------------------------------------------------

expected_methods <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

manifest_list <- vector(
  "list",
  nrow(mcar_sensitivity_grid)
)

failure_list <- list()
failure_index <- 1L

combined_result_list <- vector(
  "list",
  nrow(mcar_sensitivity_grid)
)

combined_rng_list <- vector(
  "list",
  nrow(mcar_sensitivity_grid)
)

combined_dataset_qc_list <- vector(
  "list",
  nrow(mcar_sensitivity_grid)
)

for (s in seq_len(nrow(mcar_sensitivity_grid))) {

  row <- mcar_sensitivity_grid[s, , drop = FALSE]
  sensitivity_scenario_id <- row$sensitivity_scenario_id

  cat(
    "Combining MCAR scenario ",
    sensitivity_scenario_id,
    " / ",
    nrow(mcar_sensitivity_grid),
    ": ",
    row$scenario_label,
    "\n",
    sep = ""
  )

  scenario_chunk_dir <- file.path(
    mcar_chunk_dir,
    sprintf(
      "scenario_%02d_%s",
      sensitivity_scenario_id,
      row$scenario_key
    )
  )

  expected_files <- file.path(
    scenario_chunk_dir,
    sprintf(
      "mcar_scenario_%02d_chunk_%03d.rds",
      sensitivity_scenario_id,
      seq_len(chunks_per_mcar_scenario)
    )
  )

  missing_files <- expected_files[
    !file.exists(expected_files)
  ]

  if (length(missing_files) > 0L) {

    manifest_list[[s]] <- bind_cols(
      row,
      data.frame(
        chunks_expected = chunks_per_mcar_scenario,
        chunks_found =
          chunks_per_mcar_scenario - length(missing_files),
        repetitions_found = NA_integer_,
        result_rows = NA_integer_,
        complete = FALSE,
        problem = paste0(
          "Missing ",
          length(missing_files),
          " chunk file(s)"
        ),
        stringsAsFactors = FALSE
      )
    )

    next
  }

  chunk_objects <- lapply(
    expected_files,
    readRDS
  )

  results <- bind_rows(
    lapply(chunk_objects, function(x) x$results)
  )

  rng_seeds <- bind_rows(
    lapply(chunk_objects, function(x) x$rng_states)
  )

  dataset_qc <- bind_rows(
    lapply(chunk_objects, function(x) x$dataset_qc)
  )

  repetitions_found <- sort(unique(results$repetition))
  expected_repetitions <- seq_len(nsim_mcar)

  duplicate_keys <- results %>%
    count(
      sensitivity_scenario_id,
      repetition,
      method,
      name = "n_rows"
    ) %>%
    filter(n_rows != 1L)

  dataset_duplicate_keys <- dataset_qc %>%
    count(
      sensitivity_scenario_id,
      repetition,
      name = "n_rows"
    ) %>%
    filter(n_rows != 1L)

  rng_duplicate_keys <- rng_seeds %>%
    count(
      sensitivity_scenario_id,
      repetition,
      name = "n_rows"
    ) %>%
    filter(n_rows != 1L)

  expected_result_rows <-
    nsim_mcar * length(expected_methods)

  metadata_matches <- vapply(
    chunk_objects,
    function(x) {
      identical(
        x$metadata$sensitivity_scenario_id,
        sensitivity_scenario_id
      ) &&
        identical(
          x$metadata$primary_scenario_id,
          row$primary_scenario_id
        ) &&
        identical(
          x$metadata$nsim_mcar,
          nsim_mcar
        ) &&
        identical(
          x$metadata$reps_per_mcar_chunk,
          reps_per_mcar_chunk
        ) &&
        identical(x$metadata$selection, "mcar") &&
        identical(x$metadata$marker_error, "mvn") &&
        identical(x$metadata$nimp, nimp_primary) &&
        identical(x$metadata$mice_method, mice_method_primary) &&
        identical(x$metadata$mice_maxit, mice_maxit_primary) &&
        identical(x$metadata$jomo_nburn, jomo_nburn_primary) &&
        identical(x$metadata$jomo_nbetween, jomo_nbetween_primary)
    },
    logical(1)
  )

  complete <-
    setequal(repetitions_found, expected_repetitions) &&
    setequal(unique(dataset_qc$repetition), expected_repetitions) &&
    setequal(unique(rng_seeds$repetition), expected_repetitions) &&
    nrow(results) == expected_result_rows &&
    nrow(dataset_qc) == nsim_mcar &&
    nrow(rng_seeds) == nsim_mcar &&
    nrow(duplicate_keys) == 0L &&
    nrow(dataset_duplicate_keys) == 0L &&
    nrow(rng_duplicate_keys) == 0L &&
    setequal(unique(results$method), expected_methods) &&
    mcar_all_true(results$selection == "mcar") &&
    mcar_all_true(results$marker_error == "mvn") &&
    mcar_all_true(results$nimp == nimp_primary) &&
    mcar_all_true(results$mice_method == mice_method_primary) &&
    mcar_all_true(results$mice_maxit == mice_maxit_primary) &&
    mcar_all_true(results$jomo_nburn == jomo_nburn_primary) &&
    mcar_all_true(results$jomo_nbetween == jomo_nbetween_primary) &&
    mcar_all_true(dataset_qc$marker_blockwise) &&
    mcar_all_true(dataset_qc$phase2_matches_markers) &&
    mcar_all_true(dataset_qc$constant_true_probability) &&
    mcar_all_true(dataset_qc$sufficient_phase2_n) &&
    mcar_all_true(metadata_matches)

  problems <- character()

  if (!setequal(repetitions_found, expected_repetitions)) {
    problems <- c(problems, "repetition set incomplete")
  }

  if (nrow(results) != expected_result_rows) {
    problems <- c(
      problems,
      paste0(
        "expected ",
        expected_result_rows,
        " result rows, found ",
        nrow(results)
      )
    )
  }

  if (nrow(duplicate_keys) > 0L) {
    problems <- c(problems, "duplicate repetition-method keys")
  }

  if (
    nrow(dataset_duplicate_keys) > 0L ||
      nrow(rng_duplicate_keys) > 0L
  ) {
    problems <- c(problems, "duplicate dataset-QC or RNG keys")
  }

  if (!setequal(unique(results$method), expected_methods)) {
    problems <- c(problems, "unexpected or missing method labels")
  }

  if (!mcar_all_true(metadata_matches)) {
    problems <- c(problems, "chunk metadata mismatch")
  }

  if (
    !mcar_all_true(dataset_qc$marker_blockwise) ||
      !mcar_all_true(dataset_qc$phase2_matches_markers) ||
      !mcar_all_true(dataset_qc$constant_true_probability) ||
      !mcar_all_true(dataset_qc$sufficient_phase2_n)
  ) {
    problems <- c(problems, "dataset structural QC failure")
  }

  if (complete) {

    mcar_atomic_save_rds(
      results,
      file.path(
        mcar_scenario_dir,
        sprintf(
          "mcar_scenario_%02d_estimates.rds",
          sensitivity_scenario_id
        )
      ),
      compress = "xz"
    )

    mcar_atomic_save_rds(
      rng_seeds,
      file.path(
        mcar_rng_dir,
        sprintf(
          "mcar_scenario_%02d_rng_seeds.rds",
          sensitivity_scenario_id
        )
      ),
      compress = "xz"
    )

    mcar_atomic_save_rds(
      dataset_qc,
      file.path(
        mcar_dataset_qc_dir,
        sprintf(
          "mcar_scenario_%02d_dataset_qc.rds",
          sensitivity_scenario_id
        )
      ),
      compress = "xz"
    )
  }

  failed_rows <- results %>%
    filter(
      status != "ok" |
        !is.finite(estimate) |
        !is.finite(se) |
        se <= 0
    )

  if (nrow(failed_rows) > 0L) {

    failure_list[[failure_index]] <- failed_rows %>%
      select(
        sensitivity_scenario_id,
        primary_scenario_id,
        scenario_label,
        scenario_key,
        repetition,
        method,
        status,
        message,
        N,
        target_phase2_fraction,
        r2_a_marker,
        r2_y_marker,
        theta_true
      )

    failure_index <- failure_index + 1L
  }

  combined_result_list[[s]] <- results
  combined_rng_list[[s]] <- rng_seeds
  combined_dataset_qc_list[[s]] <- dataset_qc

  manifest_list[[s]] <- bind_cols(
    row,
    data.frame(
      chunks_expected = chunks_per_mcar_scenario,
      chunks_found = chunks_per_mcar_scenario,
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

mcar_atomic_write_csv(
  manifest,
  file.path(
    mcar_combined_dir,
    "phase2_mcar_production_completion_manifest.csv"
  )
)

failures <- if (length(failure_list) == 0L) {
  tibble::tibble(
    sensitivity_scenario_id = integer(),
    primary_scenario_id = integer(),
    scenario_label = character(),
    scenario_key = character(),
    repetition = integer(),
    method = character(),
    status = character(),
    message = character(),
    N = integer(),
    target_phase2_fraction = numeric(),
    r2_a_marker = numeric(),
    r2_y_marker = numeric(),
    theta_true = numeric()
  )
} else {
  bind_rows(failure_list)
}

mcar_atomic_write_csv(
  failures,
  file.path(
    mcar_combined_dir,
    "phase2_mcar_production_method_failures.csv"
  )
)

if (
  nrow(manifest) != nrow(mcar_sensitivity_grid) ||
    any(!manifest$complete)
) {

  cat(
    "\nMCAR production is incomplete.\n",
    "See:\n",
    file.path(
      mcar_combined_dir,
      "phase2_mcar_production_completion_manifest.csv"
    ),
    "\n",
    sep = ""
  )

  quit(
    save = "no",
    status = 2L
  )
}


# ------------------------------------------------------------------------------
# 3. Reproduce the validated canary from production repetitions 1-2
# ------------------------------------------------------------------------------

canary_results <- read_csv(
  canary_results_file,
  show_col_types = FALSE
)

comparison_primary_ids <- c(168L, 208L)

production_canary_rows <- bind_rows(combined_result_list) %>%
  filter(
    primary_scenario_id %in% comparison_primary_ids,
    repetition %in% 1:2
  )

canary_rows <- canary_results %>%
  filter(
    primary_scenario_id %in% comparison_primary_ids,
    repetition %in% 1:2
  )

excluded_comparison_columns <- c(
  "elapsed_seconds",
  "total_all_methods_elapsed_seconds",
  "task_id",
  "chunk_id",
  "mice_method",
  "message"
)

comparison_columns <- setdiff(
  intersect(
    names(canary_rows),
    names(production_canary_rows)
  ),
  excluded_comparison_columns
)

comparison_key <- c(
  "primary_scenario_id",
  "repetition",
  "method"
)

comparison_columns <- unique(
  c(
    comparison_key,
    comparison_columns
  )
)

canary_reproduction_list <- vector(
  "list",
  length(comparison_primary_ids)
)

for (i in seq_along(comparison_primary_ids)) {

  primary_id <- comparison_primary_ids[[i]]

  canary_subset <- canary_rows %>%
    filter(primary_scenario_id == primary_id) %>%
    select(all_of(comparison_columns)) %>%
    arrange(primary_scenario_id, repetition, method)

  production_subset <- production_canary_rows %>%
    filter(primary_scenario_id == primary_id) %>%
    select(all_of(comparison_columns)) %>%
    arrange(primary_scenario_id, repetition, method)

  rownames(canary_subset) <- NULL
  rownames(production_subset) <- NULL

  rows_complete <-
    nrow(canary_subset) == 12L &&
    nrow(production_subset) == 12L

  # Canonicalize canary/production key storage types before exact comparison.
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
    vapply(
      canary_subset,
      is.numeric,
      logical(1)
    )
  ]

  numeric_differences <- if (keys_complete) unlist(
    lapply(
      numeric_columns,
      function(column_name) {
        difference <- abs(
          production_subset[[column_name]] -
            canary_subset[[column_name]]
        )
        difference[is.finite(difference)]
      }
    ),
    use.names = FALSE
  ) else numeric()

  max_abs_numeric_difference <- if (
    length(numeric_differences) == 0L
  ) {
    NA_real_
  } else {
    max(numeric_differences)
  }

  canary_reproduction_list[[i]] <- data.frame(
    primary_scenario_id = primary_id,
    expected_rows = 12L,
    canary_rows = nrow(canary_subset),
    production_rows = nrow(production_subset),
    max_abs_numeric_difference = max_abs_numeric_difference,
    reproduced_within_tolerance =
      rows_complete && keys_complete && isTRUE(comparison),
    comparison_message = if (
      rows_complete && keys_complete && isTRUE(comparison)
    ) {
      ""
    } else {
      paste(comparison, collapse = "; ")
    },
    stringsAsFactors = FALSE
  )
}

canary_reproduction_qc <- bind_rows(
  canary_reproduction_list
)

mcar_atomic_write_csv(
  canary_reproduction_qc,
  file.path(
    mcar_combined_dir,
    "phase2_mcar_canary_reproduction_qc.csv"
  )
)


# ------------------------------------------------------------------------------
# 4. Verify paired Naive reproduction for all 2,000 primary repetitions
# ------------------------------------------------------------------------------

naive_reproduction_list <- vector(
  "list",
  nrow(mcar_sensitivity_grid)
)

for (i in seq_len(nrow(mcar_sensitivity_grid))) {

  row <- mcar_sensitivity_grid[i, , drop = FALSE]

  primary_result_file <- file.path(
    combined_dir,
    "scenarios",
    sprintf(
      "scenario_%03d_estimates.rds",
      row$primary_scenario_id
    )
  )

  if (!file.exists(primary_result_file)) {
    stop("Matched primary result file not found: ", primary_result_file)
  }

  primary_naive <- readRDS(primary_result_file) %>%
    filter(method == "Naive") %>%
    arrange(repetition)

  mcar_naive <- combined_result_list[[i]] %>%
    filter(method == "Naive") %>%
    arrange(repetition)

  comparison_columns <- intersect(
    c("repetition", "estimate", "se", "p_value", "conf_low", "conf_high"),
    intersect(names(primary_naive), names(mcar_naive))
  )

  comparison <- all.equal(
    mcar_naive[, comparison_columns, drop = FALSE],
    primary_naive[, comparison_columns, drop = FALSE],
    tolerance = 1e-12,
    check.attributes = FALSE
  )

  numeric_columns <- setdiff(comparison_columns, "repetition")

  rows_complete <-
    nrow(primary_naive) == nsim_mcar &&
    nrow(mcar_naive) == nsim_mcar

  keys_complete <-
    rows_complete &&
    !anyDuplicated(primary_naive$repetition) &&
    !anyDuplicated(mcar_naive$repetition) &&
    identical(primary_naive$repetition, mcar_naive$repetition)

  numeric_differences <- if (keys_complete) unlist(
    lapply(
      numeric_columns,
      function(column_name) {
        difference <- abs(
          mcar_naive[[column_name]] - primary_naive[[column_name]]
        )
        difference[is.finite(difference)]
      }
    ),
    use.names = FALSE
  ) else numeric()

  naive_reproduction_list[[i]] <- data.frame(
    sensitivity_scenario_id = row$sensitivity_scenario_id,
    primary_scenario_id = row$primary_scenario_id,
    expected_rows = nsim_mcar,
    primary_rows = nrow(primary_naive),
    mcar_rows = nrow(mcar_naive),
    max_abs_numeric_difference = if (
      length(numeric_differences) == 0L
    ) {
      NA_real_
    } else {
      max(numeric_differences)
    },
    reproduced_within_tolerance =
      rows_complete && keys_complete && isTRUE(comparison),
    comparison_message = if (
      rows_complete && keys_complete && isTRUE(comparison)
    ) {
      ""
    } else {
      paste(comparison, collapse = "; ")
    },
    stringsAsFactors = FALSE
  )
}

primary_naive_reproduction_qc <- bind_rows(naive_reproduction_list)

mcar_atomic_write_csv(
  primary_naive_reproduction_qc,
  file.path(
    mcar_combined_dir,
    "phase2_mcar_primary_naive_reproduction_qc.csv"
  )
)


# ------------------------------------------------------------------------------
# 5. Final production-content gate
# ------------------------------------------------------------------------------

production_qc <- data.frame(
  check = c(
    "All four MCAR production scenarios are complete",
    "Each scenario contains exactly 2,000 repetitions",
    "Each scenario contains exactly 12,000 method-result rows",
    "All chunk and dataset structural checks pass",
    "Production repetitions reproduce the completed canary",
    "Paired Naive results reproduce all matched primary repetitions"
  ),
  passed = c(
    nrow(manifest) == 4L && mcar_all_true(manifest$complete),
    mcar_all_true(manifest$repetitions_found == nsim_mcar),
    mcar_all_true(
      manifest$result_rows ==
        nsim_mcar * length(expected_methods)
    ),
    mcar_all_true(manifest$problem == ""),
    nrow(canary_reproduction_qc) == 2L &&
      mcar_all_true(canary_reproduction_qc$canary_rows == 12L) &&
      mcar_all_true(canary_reproduction_qc$production_rows == 12L) &&
      mcar_all_true(canary_reproduction_qc$reproduced_within_tolerance),
    nrow(primary_naive_reproduction_qc) == 4L &&
      mcar_all_true(primary_naive_reproduction_qc$primary_rows == nsim_mcar) &&
      mcar_all_true(primary_naive_reproduction_qc$mcar_rows == nsim_mcar) &&
      mcar_all_true(primary_naive_reproduction_qc$reproduced_within_tolerance)
  ),
  stringsAsFactors = FALSE
)

production_qc$overall_production_content_pass <- mcar_all_true(
  production_qc$passed
)

mcar_atomic_write_csv(
  production_qc,
  file.path(
    mcar_combined_dir,
    "phase2_mcar_production_qc.csv"
  )
)

if (!mcar_all_true(production_qc$passed)) {

  failed_checks <- production_qc$check[!production_qc$passed]

  stop(
    "Phase-2 MCAR production combination FAILED. Failed checks:\n- ",
    paste(failed_checks, collapse = "\n- "),
    "\nSee outputs in: ",
    mcar_combined_dir
  )
}

mcar_atomic_write_lines(
  c(
    "PHASE-2 MCAR PRODUCTION COMBINATION: PASS",
    paste("Completed:", Sys.time()),
    paste("Production scenarios:", nrow(manifest)),
    paste("Repetitions per scenario:", nsim_mcar),
    paste("Method-level failures recorded:", nrow(failures)),
    "Canary reproduction: PASS",
    "Primary-paired Naive reproduction: PASS"
  ),
  combine_pass_file
)

cat(
  "\nPhase-2 MCAR production combination PASSED.\n",
  "Complete scenarios: ",
  sum(manifest$complete),
  " / ",
  nrow(manifest),
  "\nMethod-level failures recorded: ",
  nrow(failures),
  "\nCombined directory: ",
  mcar_combined_dir,
  "\n",
  sep = ""
)
