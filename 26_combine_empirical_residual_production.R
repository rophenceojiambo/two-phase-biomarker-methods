################################################################################
# 26_combine_empirical_residual_production.R
#
# Validates and combines all empirical-residual production chunks.
#
# The script retains method estimates, RNG seeds, and dataset QC separately for
# each scenario. It also verifies that alternative-scenario repetitions 1 and 2
# reproduce the completed canary within numerical tolerance.
################################################################################

source("00_config.R")
source("22_empirical_residual_sensitivity_helpers.R")

library(dplyr)
library(readr)

read_positive_integer_env <- function(name, default) {
  value_text <- Sys.getenv(name, unset = as.character(default))
  value <- suppressWarnings(as.integer(value_text))
  
  if (
    length(value) != 1L ||
    is.na(value) ||
    value < 1L
  ) {
    stop(
      name,
      " must be a positive integer; received '",
      value_text,
      "'."
    )
  }
  
  value
}

all_true <- function(x) {
  isTRUE(all(x))
}

nsim_empirical <- read_positive_integer_env(
  "SIM_EMPIRICAL_NSIM",
  2000L
)

reps_per_empirical_chunk <- read_positive_integer_env(
  "SIM_EMPIRICAL_REPS_PER_CHUNK",
  10L
)

chunks_per_empirical_scenario <- as.integer(
  ceiling(
    nsim_empirical /
      reps_per_empirical_chunk
  )
)

expected_methods <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

# Set production input and combined-output paths.
empirical_results_dir <- file.path(
  results_dir,
  "empirical_residual_sensitivity"
)

empirical_canary_dir <- file.path(
  empirical_results_dir,
  "canary"
)

empirical_production_dir <- file.path(
  empirical_results_dir,
  "production"
)

empirical_cache_dir <- file.path(
  empirical_production_dir,
  "cache"
)

empirical_chunk_dir <- file.path(
  empirical_production_dir,
  "chunks"
)

empirical_combined_dir <- file.path(
  empirical_production_dir,
  "combined"
)

empirical_scenario_dir <- file.path(
  empirical_combined_dir,
  "scenarios"
)

empirical_rng_dir <- file.path(
  empirical_combined_dir,
  "rng_seeds"
)

empirical_dataset_qc_dir <- file.path(
  empirical_combined_dir,
  "dataset_qc"
)

invisible(
  lapply(
    c(
      empirical_combined_dir,
      empirical_scenario_dir,
      empirical_rng_dir,
      empirical_dataset_qc_dir
    ),
    dir.create,
    showWarnings = FALSE,
    recursive = TRUE
  )
)

combine_pass_file <- file.path(
  empirical_combined_dir,
  "EMPIRICAL_RESIDUAL_PRODUCTION_COMBINE_PASS.txt"
)

# Prevent a previous pass file from surviving a failed recombination.
if (file.exists(combine_pass_file)) {
  unlink(combine_pass_file)
}

cache_pass_file <- file.path(
  empirical_cache_dir,
  "EMPIRICAL_RESIDUAL_PRODUCTION_CACHE_PASS.txt"
)

cache_qc_file <- file.path(
  empirical_cache_dir,
  "empirical_residual_production_cache_qc.csv"
)

canary_pass_file <- file.path(
  empirical_canary_dir,
  "EMPIRICAL_RESIDUAL_CANARY_PASS.txt"
)

canary_qc_file <- file.path(
  empirical_canary_dir,
  "empirical_residual_canary_qc.csv"
)

canary_results_file <- file.path(
  empirical_canary_dir,
  "empirical_residual_canary_results.csv"
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
    "Required production gate files are missing:\n- ",
    paste(missing_gate_files, collapse = "\n- ")
  )
}

cache_qc <- read_csv(
  cache_qc_file,
  show_col_types = FALSE
)

canary_qc <- read_csv(
  canary_qc_file,
  show_col_types = FALSE
)

if (
  !all(c("passed", "overall_cache_pass") %in%
       names(cache_qc)) ||
  !all_true(cache_qc$passed) ||
  !all_true(cache_qc$overall_cache_pass)
) {
  stop("The production-cache QC file does not record a full pass.")
}

if (
  !all(c("passed", "overall_canary_pass") %in%
       names(canary_qc)) ||
  !all_true(canary_qc$passed) ||
  !all_true(canary_qc$overall_canary_pass)
) {
  stop("The canary QC file does not record a full pass.")
}

n_scenarios <- nrow(empirical_sensitivity_grid)

manifest_list <- vector("list", n_scenarios)
combined_result_list <- vector("list", n_scenarios)
combined_rng_list <- vector("list", n_scenarios)
combined_dataset_qc_list <- vector(
  "list",
  n_scenarios
)

failure_list <- list()
failure_index <- 0L

# Validate and combine each scenario.
for (scenario_index in seq_len(n_scenarios)) {
  scenario <- empirical_sensitivity_grid[
    scenario_index,
    ,
    drop = FALSE
  ]
  
  scenario_id <- scenario$sensitivity_scenario_id
  
  cat(
    "Combining empirical scenario ",
    scenario_id,
    " of ",
    n_scenarios,
    ": ",
    scenario$scenario_label,
    "\n",
    sep = ""
  )
  
  scenario_chunk_dir <- file.path(
    empirical_chunk_dir,
    sprintf(
      "scenario_%02d_%s",
      scenario_id,
      scenario$scenario_key
    )
  )
  
  expected_files <- file.path(
    scenario_chunk_dir,
    sprintf(
      "empirical_scenario_%02d_chunk_%03d.rds",
      scenario_id,
      seq_len(chunks_per_empirical_scenario)
    )
  )
  
  missing_files <- expected_files[
    !file.exists(expected_files)
  ]
  
  if (length(missing_files) > 0L) {
    manifest_list[[scenario_index]] <- bind_cols(
      scenario,
      data.frame(
        chunks_expected =
          chunks_per_empirical_scenario,
        chunks_found = (
          chunks_per_empirical_scenario -
            length(missing_files)
        ),
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
    function(path) {
      tryCatch(
        readRDS(path),
        error = function(e) {
          stop(
            "Could not read chunk:\n",
            path,
            "\n",
            conditionMessage(e)
          )
        }
      )
    }
  )
  
  valid_chunk_objects <- vapply(
    chunk_objects,
    function(object) {
      is.list(object) &&
        all(c(
          "metadata",
          "results",
          "rng_states",
          "dataset_qc"
        ) %in% names(object))
    },
    logical(1)
  )
  
  if (!all_true(valid_chunk_objects)) {
    manifest_list[[scenario_index]] <- bind_cols(
      scenario,
      data.frame(
        chunks_expected =
          chunks_per_empirical_scenario,
        chunks_found =
          sum(valid_chunk_objects),
        repetitions_found = NA_integer_,
        result_rows = NA_integer_,
        complete = FALSE,
        problem =
          "One or more chunk objects are malformed",
        stringsAsFactors = FALSE
      )
    )
    
    next
  }
  
  results <- bind_rows(
    lapply(chunk_objects, function(x) x$results)
  )
  
  rng_seeds <- bind_rows(
    lapply(chunk_objects, function(x) x$rng_states)
  )
  
  dataset_qc <- bind_rows(
    lapply(chunk_objects, function(x) x$dataset_qc)
  )
  
  required_result_columns <- c(
    "sensitivity_scenario_id",
    "primary_scenario_id",
    "scenario_label",
    "scenario_key",
    "repetition",
    "method",
    "status",
    "message",
    "estimate",
    "se",
    "marker_error",
    "nimp",
    "mice_method",
    "mice_maxit",
    "jomo_nburn",
    "jomo_nbetween"
  )
  
  required_rng_columns <- c(
    "sensitivity_scenario_id",
    "primary_scenario_id",
    "repetition",
    "data_seed",
    "method_seed"
  )
  
  required_dataset_columns <- c(
    "sensitivity_scenario_id",
    "primary_scenario_id",
    "repetition",
    "marker_blockwise",
    "phase2_matches_markers",
    "sufficient_phase2_n",
    "finite_probability_range"
  )
  
  missing_result_columns <- setdiff(
    required_result_columns,
    names(results)
  )
  
  missing_rng_columns <- setdiff(
    required_rng_columns,
    names(rng_seeds)
  )
  
  missing_dataset_columns <- setdiff(
    required_dataset_columns,
    names(dataset_qc)
  )
  
  if (
    length(missing_result_columns) > 0L ||
    length(missing_rng_columns) > 0L ||
    length(missing_dataset_columns) > 0L
  ) {
    missing_text <- c(
      if (length(missing_result_columns) > 0L) {
        paste0(
          "results: ",
          paste(
            missing_result_columns,
            collapse = ", "
          )
        )
      },
      if (length(missing_rng_columns) > 0L) {
        paste0(
          "RNG: ",
          paste(missing_rng_columns, collapse = ", ")
        )
      },
      if (length(missing_dataset_columns) > 0L) {
        paste0(
          "dataset QC: ",
          paste(
            missing_dataset_columns,
            collapse = ", "
          )
        )
      }
    )
    
    manifest_list[[scenario_index]] <- bind_cols(
      scenario,
      data.frame(
        chunks_expected =
          chunks_per_empirical_scenario,
        chunks_found =
          chunks_per_empirical_scenario,
        repetitions_found = NA_integer_,
        result_rows = nrow(results),
        complete = FALSE,
        problem = paste(
          missing_text,
          collapse = "; "
        ),
        stringsAsFactors = FALSE
      )
    )
    
    next
  }
  
  expected_repetitions <- seq_len(nsim_empirical)
  repetitions_found <- sort(unique(results$repetition))
  
  result_keys <- results %>%
    count(
      sensitivity_scenario_id,
      repetition,
      method,
      name = "n_rows"
    )
  
  dataset_keys <- dataset_qc %>%
    count(
      sensitivity_scenario_id,
      repetition,
      name = "n_rows"
    )
  
  rng_keys <- rng_seeds %>%
    count(
      sensitivity_scenario_id,
      repetition,
      name = "n_rows"
    )
  
  expected_result_rows <- (
    nsim_empirical * length(expected_methods)
  )
  
  metadata_matches <- vapply(
    seq_along(chunk_objects),
    function(chunk_index) {
      metadata <- chunk_objects[[chunk_index]]$metadata
      
      expected_start <- (
        (chunk_index - 1L) *
          reps_per_empirical_chunk
      ) + 1L
      
      expected_end <- min(
        chunk_index *
          reps_per_empirical_chunk,
        nsim_empirical
      )
      
      identical(
        metadata$sensitivity_scenario_id,
        scenario_id
      ) &&
        identical(
          metadata$primary_scenario_id,
          scenario$primary_scenario_id
        ) &&
        identical(
          metadata$chunk_id,
          as.integer(chunk_index)
        ) &&
        identical(
          metadata$rep_start,
          expected_start
        ) &&
        identical(
          metadata$rep_end,
          expected_end
        ) &&
        identical(
          metadata$nsim_empirical,
          nsim_empirical
        ) &&
        identical(
          metadata$reps_per_empirical_chunk,
          reps_per_empirical_chunk
        ) &&
        identical(metadata$marker_error, "empirical") &&
        identical(metadata$nimp, nimp_primary) &&
        identical(
          metadata$mice_method,
          mice_method_primary
        ) &&
        identical(
          metadata$mice_maxit,
          mice_maxit_primary
        ) &&
        identical(
          metadata$jomo_nburn,
          jomo_nburn_primary
        ) &&
        identical(
          metadata$jomo_nbetween,
          jomo_nbetween_primary
        )
    },
    logical(1)
  )
  
  checks <- c(
    result_repetitions_complete = setequal(
      repetitions_found,
      expected_repetitions
    ),
    dataset_repetitions_complete = setequal(
      dataset_qc$repetition,
      expected_repetitions
    ),
    rng_repetitions_complete = setequal(
      rng_seeds$repetition,
      expected_repetitions
    ),
    result_row_count = (
      nrow(results) == expected_result_rows
    ),
    dataset_row_count = (
      nrow(dataset_qc) == nsim_empirical
    ),
    rng_row_count = (
      nrow(rng_seeds) == nsim_empirical
    ),
    unique_result_keys = (
      nrow(result_keys) == expected_result_rows &&
        all_true(result_keys$n_rows == 1L)
    ),
    unique_dataset_keys = (
      nrow(dataset_keys) == nsim_empirical &&
        all_true(dataset_keys$n_rows == 1L)
    ),
    unique_rng_keys = (
      nrow(rng_keys) == nsim_empirical &&
        all_true(rng_keys$n_rows == 1L)
    ),
    expected_methods = setequal(
      results$method,
      expected_methods
    ),
    correct_result_setting = all_true(
      results$sensitivity_scenario_id ==
        scenario_id &
        results$primary_scenario_id ==
        scenario$primary_scenario_id
    ),
    correct_dataset_setting = all_true(
      dataset_qc$sensitivity_scenario_id ==
        scenario_id &
        dataset_qc$primary_scenario_id ==
        scenario$primary_scenario_id
    ),
    correct_rng_setting = all_true(
      rng_seeds$sensitivity_scenario_id ==
        scenario_id &
        rng_seeds$primary_scenario_id ==
        scenario$primary_scenario_id
    ),
    empirical_marker_error = all_true(
      results$marker_error == "empirical"
    ),
    production_mi_settings = all_true(
      results$nimp == nimp_primary &
        results$mice_method ==
        mice_method_primary &
        results$mice_maxit ==
        mice_maxit_primary &
        results$jomo_nburn ==
        jomo_nburn_primary &
        results$jomo_nbetween ==
        jomo_nbetween_primary
    ),
    marker_missingness_blockwise = all_true(
      dataset_qc$marker_blockwise
    ),
    phase2_matches_markers = all_true(
      dataset_qc$phase2_matches_markers
    ),
    sufficient_phase2_n = all_true(
      dataset_qc$sufficient_phase2_n
    ),
    finite_probability_range = all_true(
      dataset_qc$finite_probability_range
    ),
    chunk_metadata_match = all_true(
      metadata_matches
    )
  )
  
  complete <- all_true(checks)
  
  problem_labels <- c(
    result_repetitions_complete =
      "result repetition set incomplete",
    dataset_repetitions_complete =
      "dataset-QC repetition set incomplete",
    rng_repetitions_complete =
      "RNG repetition set incomplete",
    result_row_count =
      "unexpected result-row count",
    dataset_row_count =
      "unexpected dataset-QC row count",
    rng_row_count =
      "unexpected RNG row count",
    unique_result_keys =
      "duplicate or missing repetition-method keys",
    unique_dataset_keys =
      "duplicate or missing dataset-QC keys",
    unique_rng_keys =
      "duplicate or missing RNG keys",
    expected_methods =
      "unexpected or missing methods",
    correct_result_setting =
      "incorrect setting in results",
    correct_dataset_setting =
      "incorrect setting in dataset QC",
    correct_rng_setting =
      "incorrect setting in RNG data",
    empirical_marker_error =
      "non-empirical marker-error label",
    production_mi_settings =
      "production MI settings do not match",
    marker_missingness_blockwise =
      "marker missingness is not blockwise",
    phase2_matches_markers =
      "phase2 disagrees with marker observation",
    sufficient_phase2_n =
      "insufficient Phase-2 sample",
    finite_probability_range =
      "invalid true selection probability",
    chunk_metadata_match =
      "chunk metadata mismatch"
  )
  
  problems <- unname(
    problem_labels[names(checks)[!checks]]
  )
  
  # Save validated scenario-level components.
  if (complete) {
    saveRDS(
      results,
      file.path(
        empirical_scenario_dir,
        sprintf(
          "empirical_scenario_%02d_estimates.rds",
          scenario_id
        )
      ),
      compress = "xz"
    )
    
    saveRDS(
      rng_seeds,
      file.path(
        empirical_rng_dir,
        sprintf(
          "empirical_scenario_%02d_rng_seeds.rds",
          scenario_id
        )
      ),
      compress = "xz"
    )
    
    saveRDS(
      dataset_qc,
      file.path(
        empirical_dataset_qc_dir,
        sprintf(
          "empirical_scenario_%02d_dataset_qc.rds",
          scenario_id
        )
      ),
      compress = "xz"
    )
  }
  
  # Record estimator failures without treating them as structural failure.
  failed_rows <- results %>%
    filter(
      is.na(status) |
        status != "ok" |
        !is.finite(estimate) |
        !is.finite(se) |
        se <= 0
    )
  
  if (nrow(failed_rows) > 0L) {
    failure_index <- failure_index + 1L
    
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
  }
  
  combined_result_list[[scenario_index]] <- results
  combined_rng_list[[scenario_index]] <- rng_seeds
  combined_dataset_qc_list[[scenario_index]] <-
    dataset_qc
  
  manifest_list[[scenario_index]] <- bind_cols(
    scenario,
    data.frame(
      chunks_expected =
        chunks_per_empirical_scenario,
      chunks_found =
        chunks_per_empirical_scenario,
      repetitions_found =
        length(repetitions_found),
      result_rows = nrow(results),
      complete = complete,
      problem = paste(
        problems,
        collapse = "; "
      ),
      stringsAsFactors = FALSE
    )
  )
}

manifest <- bind_rows(manifest_list)

manifest_file <- file.path(
  empirical_combined_dir,
  "empirical_residual_production_completion_manifest.csv"
)

write_csv(manifest, manifest_file)

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

write_csv(
  failures,
  file.path(
    empirical_combined_dir,
    "empirical_residual_production_method_failures.csv"
  )
)

if (
  nrow(manifest) != n_scenarios ||
  !all_true(manifest$complete)
) {
  cat(
    "\nEmpirical production is incomplete.\n",
    "See: ",
    manifest_file,
    "\n",
    sep = ""
  )
  
  quit(save = "no", status = 2L)
}

# Verify that production repetitions 1 and 2 reproduce the canary.
canary_results <- read_csv(
  canary_results_file,
  show_col_types = FALSE
)

comparison_primary_ids <- c(168L, 208L)
comparison_repetitions <- 1:2

production_canary_rows <- bind_rows(
  combined_result_list
) %>%
  filter(
    primary_scenario_id %in%
      comparison_primary_ids,
    repetition %in% comparison_repetitions
  )

canary_rows <- canary_results %>%
  filter(
    primary_scenario_id %in%
      comparison_primary_ids,
    repetition %in% comparison_repetitions
  )

excluded_comparison_columns <- c(
  "elapsed_seconds",
  "total_all_methods_elapsed_seconds",
  "task_id",
  "chunk_id",
  "mice_method",
  "message"
)

comparison_key <- c(
  "primary_scenario_id",
  "repetition",
  "method"
)

comparison_columns <- setdiff(
  intersect(
    names(canary_rows),
    names(production_canary_rows)
  ),
  excluded_comparison_columns
)

comparison_columns <- unique(
  c(comparison_key, comparison_columns)
)

expected_comparison_rows <- (
  length(comparison_repetitions) *
    length(expected_methods)
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
  
  canary_keys <- canary_subset %>%
    count(
      primary_scenario_id,
      repetition,
      method,
      name = "n"
    )
  
  production_keys <- production_subset %>%
    count(
      primary_scenario_id,
      repetition,
      method,
      name = "n"
    )
  
  # Canonicalize comparison-key storage types before requiring exact equality.
  #
  # Canary results are read from CSV, so integer-valued identifiers are
  # imported as doubles. Production results are read from RDS and retain
  # integer storage. Convert both sides to the same explicit key types before
  # using identical(). This changes only storage representation, not values.
  canary_keys_compare <- canary_keys %>%
    transmute(
      primary_scenario_id =
        as.integer(primary_scenario_id),
      repetition =
        as.integer(repetition),
      method =
        as.character(method)
    ) %>%
    as.data.frame()

  production_keys_compare <- production_keys %>%
    transmute(
      primary_scenario_id =
        as.integer(primary_scenario_id),
      repetition =
        as.integer(repetition),
      method =
        as.character(method)
    ) %>%
    as.data.frame()

  keys_valid <- (
    nrow(canary_subset) ==
      expected_comparison_rows &&
      nrow(production_subset) ==
      expected_comparison_rows &&
      nrow(canary_keys) ==
      expected_comparison_rows &&
      nrow(production_keys) ==
      expected_comparison_rows &&
      all_true(canary_keys$n == 1L) &&
      all_true(production_keys$n == 1L) &&
      identical(
        canary_keys_compare,
        production_keys_compare
      )
  )
  
  rownames(canary_subset) <- NULL
  rownames(production_subset) <- NULL
  
  comparison <- if (keys_valid) {
    all.equal(
      production_subset,
      canary_subset,
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  } else {
    "Canary or production comparison keys are incomplete."
  }
  
  max_abs_numeric_difference <- NA_real_
  
  if (
    keys_valid &&
    nrow(canary_subset) ==
    nrow(production_subset)
  ) {
    numeric_columns <- comparison_columns[
      vapply(
        canary_subset,
        is.numeric,
        logical(1)
      ) &
        vapply(
          production_subset,
          is.numeric,
          logical(1)
        )
    ]
    
    numeric_differences <- unlist(
      lapply(
        numeric_columns,
        function(column_name) {
          difference <- abs(
            production_subset[[column_name]] -
              canary_subset[[column_name]]
          )
          
          difference[
            is.finite(difference)
          ]
        }
      ),
      use.names = FALSE
    )
    
    max_abs_numeric_difference <- if (
      length(numeric_differences) == 0L
    ) {
      0
    } else {
      max(numeric_differences)
    }
  }
  
  canary_reproduction_list[[i]] <- data.frame(
    primary_scenario_id = primary_id,
    expected_rows = expected_comparison_rows,
    canary_rows = nrow(canary_subset),
    production_rows = nrow(production_subset),
    max_abs_numeric_difference =
      max_abs_numeric_difference,
    reproduced_within_tolerance =
      isTRUE(comparison),
    comparison_message = if (
      isTRUE(comparison)
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

write_csv(
  canary_reproduction_qc,
  file.path(
    empirical_combined_dir,
    "empirical_residual_canary_reproduction_qc.csv"
  )
)

# Make the final production-content decision.
production_qc <- data.frame(
  check = c(
    "All empirical production scenarios are complete",
    "Each scenario has the expected repetitions",
    "Each scenario has the expected method-result rows",
    "All chunk and dataset structural checks pass",
    "Production repetitions reproduce the canary"
  ),
  passed = c(
    nrow(manifest) == n_scenarios &&
      all_true(manifest$complete),
    all_true(
      manifest$repetitions_found ==
        nsim_empirical
    ),
    all_true(
      manifest$result_rows ==
        nsim_empirical *
        length(expected_methods)
    ),
    all_true(manifest$problem == ""),
    nrow(canary_reproduction_qc) ==
      length(comparison_primary_ids) &&
      all_true(
        canary_reproduction_qc$canary_rows ==
          expected_comparison_rows
      ) &&
      all_true(
        canary_reproduction_qc$production_rows ==
          expected_comparison_rows
      ) &&
      all_true(
        canary_reproduction_qc$
          reproduced_within_tolerance
      )
  ),
  stringsAsFactors = FALSE
)

production_qc$overall_production_content_pass <-
  all_true(production_qc$passed)

write_csv(
  production_qc,
  file.path(
    empirical_combined_dir,
    "empirical_residual_production_qc.csv"
  )
)

if (!all_true(production_qc$passed)) {
  failed_checks <- production_qc$check[
    !production_qc$passed
  ]
  
  stop(
    "Empirical residual production combination failed. Failed checks:\n- ",
    paste(failed_checks, collapse = "\n- "),
    "\nSee outputs in: ",
    empirical_combined_dir
  )
}

writeLines(
  c(
    "EMPIRICAL RESIDUAL PRODUCTION COMBINATION: PASS",
    paste("Completed:", Sys.time()),
    paste("Production scenarios:", n_scenarios),
    paste(
      "Repetitions per scenario:",
      nsim_empirical
    ),
    paste(
      "Method-level failures recorded:",
      nrow(failures)
    ),
    "Canary reproduction: PASS"
  ),
  combine_pass_file
)

cat(
  "\nEmpirical residual production combination passed.\n",
  "Complete scenarios: ",
  sum(manifest$complete),
  " of ",
  nrow(manifest),
  "\nMethod-level failures: ",
  nrow(failures),
  "\nCombined directory: ",
  empirical_combined_dir,
  "\n",
  sep = ""
)
