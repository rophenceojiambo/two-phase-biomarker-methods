################################################################################
# 42_combine_ipw_weight_capping_results.R
#
# Validates and combines all IPW weight-capping production chunks. Production
# repetitions 1-2 in the two canary scenarios must reproduce the canary exactly.
################################################################################

source("00_config.R")
source("39_ipw_weight_capping_helpers.R")

library(dplyr)
library(readr)


# ------------------------------------------------------------------------------
# 1. Settings, paths, and gates
# ------------------------------------------------------------------------------

nsim_ipw_cap <- read_ipw_cap_positive_integer_env(
  "SIM_IPW_CAP_NSIM",
  2000L
)

reps_per_ipw_cap_chunk <- read_ipw_cap_positive_integer_env(
  "SIM_IPW_CAP_REPS_PER_CHUNK",
  100L
)

if (nsim_ipw_cap != 2000L || reps_per_ipw_cap_chunk != 100L) {
  stop(
    "The prespecified IPW-capping combination requires ",
    "SIM_IPW_CAP_NSIM=2000 and SIM_IPW_CAP_REPS_PER_CHUNK=100."
  )
}

chunks_per_ipw_cap_scenario <- as.integer(
  ceiling(nsim_ipw_cap / reps_per_ipw_cap_chunk)
)

ipw_cap_results_dir <- file.path(results_dir, "ipw_weight_capping_sensitivity")
ipw_cap_canary_dir <- file.path(ipw_cap_results_dir, "canary")
ipw_cap_production_dir <- file.path(ipw_cap_results_dir, "production")
ipw_cap_chunk_dir <- file.path(ipw_cap_production_dir, "chunks")
ipw_cap_combined_dir <- file.path(ipw_cap_production_dir, "combined")
ipw_cap_scenario_dir <- file.path(ipw_cap_combined_dir, "scenarios")
ipw_cap_rng_dir <- file.path(ipw_cap_combined_dir, "rng_seeds")
ipw_cap_dataset_qc_dir <- file.path(ipw_cap_combined_dir, "dataset_qc")

invisible(
  lapply(
    c(
      ipw_cap_combined_dir,
      ipw_cap_scenario_dir,
      ipw_cap_rng_dir,
      ipw_cap_dataset_qc_dir
    ),
    dir.create,
    showWarnings = FALSE,
    recursive = TRUE
  )
)

combine_pass_file <- file.path(
  ipw_cap_combined_dir,
  "IPW_WEIGHT_CAPPING_PRODUCTION_COMBINE_PASS.txt"
)

if (file.exists(combine_pass_file)) {
  unlink(combine_pass_file)
}

fcs_pmm_summary_pass_file <- file.path(
  results_dir,
  "fcs_pmm_sensitivity",
  "production",
  "summary",
  "FCS_PMM_PRODUCTION_SUMMARY_PASS.txt"
)

fcs_pmm_summary_qc_file <- file.path(
  results_dir,
  "fcs_pmm_sensitivity",
  "production",
  "summary",
  "fcs_pmm_summary_qc.csv"
)

canary_pass_file <- file.path(
  ipw_cap_canary_dir,
  "IPW_WEIGHT_CAPPING_CANARY_PASS.txt"
)

canary_results_file <- file.path(
  ipw_cap_canary_dir,
  "ipw_weight_capping_canary_results.csv"
)

canary_qc_file <- file.path(
  ipw_cap_canary_dir,
  "ipw_weight_capping_canary_qc.csv"
)

canary_specification_file <- file.path(
  ipw_cap_canary_dir,
  "ipw_weight_capping_specification.csv"
)

required_gate_files <- c(
  fcs_pmm_summary_pass_file,
  fcs_pmm_summary_qc_file,
  canary_pass_file,
  canary_results_file,
  canary_qc_file,
  canary_specification_file
)

missing_gate_files <- required_gate_files[!file.exists(required_gate_files)]

if (length(missing_gate_files) > 0L) {
  stop(
    "Required IPW-capping production gate files are missing:\n- ",
    paste(missing_gate_files, collapse = "\n- ")
  )
}

if (!any(grepl(
  "FCS-PMM PRODUCTION SUMMARY: PASS",
  readLines(fcs_pmm_summary_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The FCS-PMM production summary PASS marker is invalid.")
}

fcs_pmm_summary_qc <- read_csv(
  fcs_pmm_summary_qc_file,
  show_col_types = FALSE
)

if (
  nrow(fcs_pmm_summary_qc) == 0L ||
    !ipw_cap_all_true(fcs_pmm_summary_qc$passed) ||
    !ipw_cap_all_true(fcs_pmm_summary_qc$overall_summary_pass)
) {
  stop("The FCS-PMM summary QC file does not record a full pass.")
}

if (!any(grepl(
  "IPW WEIGHT CAPPING CANARY: PASS",
  readLines(canary_pass_file, warn = FALSE),
  fixed = TRUE
))) {
  stop("The IPW weight-capping canary PASS marker is invalid.")
}

canary_qc <- read_csv(canary_qc_file, show_col_types = FALSE)
canary_specification <- read_csv(
  canary_specification_file,
  show_col_types = FALSE
)

if (
  nrow(canary_qc) == 0L ||
    !ipw_cap_all_true(canary_qc$passed) ||
    !ipw_cap_all_true(canary_qc$overall_canary_pass)
) {
  stop("The IPW-capping canary QC file does not record a full pass.")
}

current_specification <- data.frame(
  lower_probability = ipw_cap_lower_probability,
  upper_probability = ipw_cap_upper_probability,
  quantile_type = ipw_cap_quantile_type,
  threshold_population = "estimated raw weights among Phase-2 observations",
  operation = "two-sided winsorization",
  stringsAsFactors = FALSE
)

if (!isTRUE(all.equal(
  current_specification,
  canary_specification,
  tolerance = 0,
  check.attributes = FALSE
))) {
  stop("The current IPW cap specification differs from the validated canary.")
}


# ------------------------------------------------------------------------------
# 2. Combine and validate each scenario
# ------------------------------------------------------------------------------

manifest_list <- vector("list", nrow(ipw_cap_sensitivity_grid))
combined_result_list <- vector("list", nrow(ipw_cap_sensitivity_grid))
failure_list <- list()
failure_index <- 1L

for (s in seq_len(nrow(ipw_cap_sensitivity_grid))) {

  row <- ipw_cap_sensitivity_grid[s, , drop = FALSE]
  sensitivity_scenario_id <- row$sensitivity_scenario_id

  cat(
    "Combining IPW-capping scenario ",
    sensitivity_scenario_id,
    " / ",
    nrow(ipw_cap_sensitivity_grid),
    ": ",
    row$scenario_label,
    "\n",
    sep = ""
  )

  scenario_chunk_dir <- file.path(
    ipw_cap_chunk_dir,
    sprintf(
      "scenario_%02d_%s",
      sensitivity_scenario_id,
      row$scenario_key
    )
  )

  expected_files <- file.path(
    scenario_chunk_dir,
    sprintf(
      "ipw_cap_scenario_%02d_chunk_%02d.rds",
      sensitivity_scenario_id,
      seq_len(chunks_per_ipw_cap_scenario)
    )
  )

  missing_files <- expected_files[!file.exists(expected_files)]

  if (length(missing_files) > 0L) {

    manifest_list[[s]] <- bind_cols(
      row,
      data.frame(
        chunks_expected = chunks_per_ipw_cap_scenario,
        chunks_found =
          chunks_per_ipw_cap_scenario - length(missing_files),
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
    function(path) tryCatch(readRDS(path), error = function(e) NULL)
  )

  readable_chunks <- !vapply(chunk_objects, is.null, logical(1))

  if (!ipw_cap_all_true(readable_chunks)) {
    manifest_list[[s]] <- bind_cols(
      row,
      data.frame(
        chunks_expected = chunks_per_ipw_cap_scenario,
        chunks_found = sum(readable_chunks),
        repetitions_found = NA_integer_,
        result_rows = NA_integer_,
        complete = FALSE,
        problem = "One or more chunk files are unreadable",
        stringsAsFactors = FALSE
      )
    )
    next
  }

  required_components <- c("metadata", "results", "rng_states", "dataset_qc")
  components_present <- vapply(
    chunk_objects,
    function(x) {
      all(required_components %in% names(x)) &&
        is.list(x$metadata) &&
        is.data.frame(x$results) &&
        is.data.frame(x$rng_states) &&
        is.data.frame(x$dataset_qc)
    },
    logical(1)
  )

  if (!ipw_cap_all_true(components_present)) {
    manifest_list[[s]] <- bind_cols(
      row,
      data.frame(
        chunks_expected = chunks_per_ipw_cap_scenario,
        chunks_found = chunks_per_ipw_cap_scenario,
        repetitions_found = NA_integer_,
        result_rows = NA_integer_,
        complete = FALSE,
        problem = "One or more chunks have missing or invalid components",
        stringsAsFactors = FALSE
      )
    )
    next
  }

  results <- bind_rows(lapply(chunk_objects, function(x) x$results))
  rng_seeds <- bind_rows(lapply(chunk_objects, function(x) x$rng_states))
  dataset_qc <- bind_rows(lapply(chunk_objects, function(x) x$dataset_qc))

  required_result_columns <- c(
    "sensitivity_scenario_id", "primary_scenario_id", "repetition", "method",
    "status", "estimate", "se", "selection", "marker_error",
    "cap_lower_probability", "cap_upper_probability", "cap_quantile_type",
    "n_capped_lower", "n_capped_upper", "capped_weight_min",
    "capped_weight_max", "cap_lower_threshold", "cap_upper_threshold"
  )
  required_dataset_columns <- c(
    "sensitivity_scenario_id", "repetition", "marker_blockwise",
    "phase2_matches_markers", "finite_valid_probabilities",
    "sufficient_phase2_n"
  )
  required_rng_columns <- c(
    "sensitivity_scenario_id", "repetition", "data_rng_state",
    "method_rng_state"
  )

  content_columns_present <-
    all(required_result_columns %in% names(results)) &&
    all(required_dataset_columns %in% names(dataset_qc)) &&
    all(required_rng_columns %in% names(rng_seeds))

  if (!content_columns_present) {
    manifest_list[[s]] <- bind_cols(
      row,
      data.frame(
        chunks_expected = chunks_per_ipw_cap_scenario,
        chunks_found = chunks_per_ipw_cap_scenario,
        repetitions_found = NA_integer_,
        result_rows = nrow(results),
        complete = FALSE,
        problem = "Combined chunks are missing required columns",
        stringsAsFactors = FALSE
      )
    )
    next
  }

  repetitions_found <- sort(unique(results$repetition))
  expected_repetitions <- seq_len(nsim_ipw_cap)

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
    function(chunk_index) {
      x <- chunk_objects[[chunk_index]]
      expected_rep_start <-
        (chunk_index - 1L) * reps_per_ipw_cap_chunk + 1L
      expected_rep_end <- min(
        chunk_index * reps_per_ipw_cap_chunk,
        nsim_ipw_cap
      )

      identical(
        x$metadata$sensitivity_scenario_id,
        sensitivity_scenario_id
      ) &&
        identical(x$metadata$primary_scenario_id, row$primary_scenario_id) &&
        identical(x$metadata$chunk_id, as.integer(chunk_index)) &&
        identical(x$metadata$rep_start, expected_rep_start) &&
        identical(x$metadata$rep_end, expected_rep_end) &&
        identical(x$metadata$nsim_ipw_cap, nsim_ipw_cap) &&
        identical(
          x$metadata$reps_per_ipw_cap_chunk,
          reps_per_ipw_cap_chunk
        ) &&
        identical(x$metadata$selection, "mar") &&
        identical(x$metadata$marker_error, "mvn") &&
        identical(
          x$metadata$cap_probabilities,
          c(ipw_cap_lower_probability, ipw_cap_upper_probability)
        ) &&
        identical(x$metadata$cap_quantile_type, ipw_cap_quantile_type) &&
        identical(x$metadata$cap_operation, "two-sided winsorization")
    },
    logical(1)
  )

  successful_results <- results %>%
    filter(status == "ok")

  successful_cap_qc <- nrow(successful_results) == 0L || ipw_cap_all_true(
    successful_results$cap_lower_probability ==
      ipw_cap_lower_probability &
      successful_results$cap_upper_probability ==
        ipw_cap_upper_probability &
      successful_results$cap_quantile_type == ipw_cap_quantile_type &
      successful_results$n_capped_lower > 0L &
      successful_results$n_capped_upper > 0L &
      successful_results$capped_weight_min >=
        successful_results$cap_lower_threshold - 1e-12 &
      successful_results$capped_weight_max <=
        successful_results$cap_upper_threshold + 1e-12
  )

  complete <-
    setequal(repetitions_found, expected_repetitions) &&
    setequal(unique(dataset_qc$repetition), expected_repetitions) &&
    setequal(unique(rng_seeds$repetition), expected_repetitions) &&
    nrow(results) == nsim_ipw_cap &&
    nrow(dataset_qc) == nsim_ipw_cap &&
    nrow(rng_seeds) == nsim_ipw_cap &&
    nrow(duplicate_keys) == 0L &&
    nrow(dataset_duplicate_keys) == 0L &&
    nrow(rng_duplicate_keys) == 0L &&
    identical(unique(results$method), ipw_capped_label) &&
    ipw_cap_all_true(results$selection == "mar") &&
    ipw_cap_all_true(results$marker_error == "mvn") &&
    ipw_cap_all_true(dataset_qc$marker_blockwise) &&
    ipw_cap_all_true(dataset_qc$phase2_matches_markers) &&
    ipw_cap_all_true(dataset_qc$finite_valid_probabilities) &&
    ipw_cap_all_true(dataset_qc$sufficient_phase2_n) &&
    ipw_cap_all_true(metadata_matches) &&
    successful_cap_qc

  problems <- character()

  if (!setequal(repetitions_found, expected_repetitions)) {
    problems <- c(problems, "repetition set incomplete")
  }

  if (!setequal(unique(dataset_qc$repetition), expected_repetitions)) {
    problems <- c(problems, "dataset-QC repetition set incomplete")
  }

  if (!setequal(unique(rng_seeds$repetition), expected_repetitions)) {
    problems <- c(problems, "RNG repetition set incomplete")
  }

  if (nrow(results) != nsim_ipw_cap) {
    problems <- c(
      problems,
      paste0("expected ", nsim_ipw_cap, " rows, found ", nrow(results))
    )
  }

  if (nrow(dataset_qc) != nsim_ipw_cap) {
    problems <- c(problems, "unexpected dataset-QC row count")
  }

  if (nrow(rng_seeds) != nsim_ipw_cap) {
    problems <- c(problems, "unexpected RNG row count")
  }

  if (
    nrow(duplicate_keys) > 0L ||
      nrow(dataset_duplicate_keys) > 0L ||
      nrow(rng_duplicate_keys) > 0L
  ) {
    problems <- c(problems, "duplicate result, dataset-QC, or RNG keys")
  }

  if (!identical(unique(results$method), ipw_capped_label)) {
    problems <- c(problems, "unexpected method label")
  }

  if (!ipw_cap_all_true(metadata_matches)) {
    problems <- c(problems, "chunk metadata mismatch")
  }

  if (!successful_cap_qc) {
    problems <- c(problems, "successful capped-weight diagnostic failure")
  }

  if (
    !ipw_cap_all_true(dataset_qc$marker_blockwise) ||
      !ipw_cap_all_true(dataset_qc$phase2_matches_markers) ||
      !ipw_cap_all_true(dataset_qc$finite_valid_probabilities) ||
      !ipw_cap_all_true(dataset_qc$sufficient_phase2_n)
  ) {
    problems <- c(problems, "dataset structural QC failure")
  }

  if (complete) {

    ipw_cap_atomic_save_rds(
      results,
      file.path(
        ipw_cap_scenario_dir,
        sprintf(
          "ipw_cap_scenario_%02d_estimates.rds",
          sensitivity_scenario_id
        )
      ),
      compress = "xz"
    )

    ipw_cap_atomic_save_rds(
      rng_seeds,
      file.path(
        ipw_cap_rng_dir,
        sprintf(
          "ipw_cap_scenario_%02d_rng_seeds.rds",
          sensitivity_scenario_id
        )
      ),
      compress = "xz"
    )

    ipw_cap_atomic_save_rds(
      dataset_qc,
      file.path(
        ipw_cap_dataset_qc_dir,
        sprintf(
          "ipw_cap_scenario_%02d_dataset_qc.rds",
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
      chunks_expected = chunks_per_ipw_cap_scenario,
      chunks_found = chunks_per_ipw_cap_scenario,
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

ipw_cap_atomic_write_csv(
  manifest,
  file.path(ipw_cap_combined_dir, "ipw_weight_capping_completion_manifest.csv")
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

ipw_cap_atomic_write_csv(
  failures,
  file.path(ipw_cap_combined_dir, "ipw_weight_capping_method_failures.csv")
)

if (
  nrow(manifest) != nrow(ipw_cap_sensitivity_grid) ||
    !ipw_cap_all_true(manifest$complete)
) {
  cat(
    "\nIPW-capping production is incomplete. See:\n",
    file.path(
      ipw_cap_combined_dir,
      "ipw_weight_capping_completion_manifest.csv"
    ),
    "\n",
    sep = ""
  )
  quit(save = "no", status = 2L)
}


# ------------------------------------------------------------------------------
# 3. Reproduce capped-IPW canary results from production repetitions 1-2
# ------------------------------------------------------------------------------

canary_capped <- read_csv(canary_results_file, show_col_types = FALSE) %>%
  filter(method == ipw_capped_label) %>%
  arrange(primary_scenario_id, repetition)

production_capped <- bind_rows(combined_result_list) %>%
  filter(
    primary_scenario_id %in% ipw_cap_canary_grid$primary_scenario_id,
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
  intersect(names(canary_capped), names(production_capped)),
  excluded_columns
)

comparison_key <- c("primary_scenario_id", "repetition", "method")
comparison_columns <- unique(c(comparison_key, comparison_columns))

required_canary_columns <- c(
  comparison_key,
  "estimate", "se", "p_value", "conf_low", "conf_high", "status",
  "cap_lower_threshold", "cap_upper_threshold", "n_capped_lower",
  "n_capped_upper", "raw_weight_min", "raw_weight_p99", "raw_weight_max",
  "raw_weight_cv", "raw_weight_ess", "capped_weight_min",
  "capped_weight_p99", "capped_weight_max", "capped_weight_cv",
  "capped_weight_ess"
)

if (!all(required_canary_columns %in% comparison_columns)) {
  stop("Canary and production results lack required comparison columns.")
}

canary_reproduction_list <- vector(
  "list",
  nrow(ipw_cap_canary_grid)
)

for (i in seq_len(nrow(ipw_cap_canary_grid))) {

  primary_id <- ipw_cap_canary_grid$primary_scenario_id[[i]]

  canary_subset <- canary_capped %>%
    filter(primary_scenario_id == primary_id) %>%
    select(all_of(comparison_columns)) %>%
    arrange(repetition)

  production_subset <- production_capped %>%
    filter(primary_scenario_id == primary_id) %>%
    select(all_of(comparison_columns)) %>%
    arrange(repetition)

  rownames(canary_subset) <- NULL
  rownames(production_subset) <- NULL

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

  canary_keys_complete <-
    nrow(canary_subset) == 2L &&
    nrow(production_subset) == 2L &&
    !anyDuplicated(canary_keys_compare) &&
    !anyDuplicated(production_keys_compare) &&
    identical(
      canary_keys_compare,
      production_keys_compare
    )

  comparison <- if (canary_keys_complete) {
    all.equal(
      production_subset,
      canary_subset,
      tolerance = 1e-12,
      check.attributes = FALSE
    )
  } else {
    "Canary and production keys are incomplete or unequal."
  }

  numeric_columns <- names(canary_subset)[
    vapply(canary_subset, is.numeric, logical(1))
  ]

  numeric_differences <- if (canary_keys_complete) {
    unlist(
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
    )
  } else {
    numeric()
  }

  canary_reproduction_list[[i]] <- data.frame(
    primary_scenario_id = primary_id,
    expected_rows = 2L,
    canary_rows = nrow(canary_subset),
    production_rows = nrow(production_subset),
    max_abs_numeric_difference = if (!canary_keys_complete) {
      NA_real_
    } else if (length(numeric_differences) == 0L) {
      0
    } else {
      max(numeric_differences)
    },
    reproduced_within_tolerance = isTRUE(comparison),
    comparison_message = if (isTRUE(comparison)) {
      ""
    } else {
      paste(comparison, collapse = "; ")
    },
    stringsAsFactors = FALSE
  )
}

canary_reproduction_qc <- bind_rows(canary_reproduction_list)

ipw_cap_atomic_write_csv(
  canary_reproduction_qc,
  file.path(ipw_cap_combined_dir, "ipw_cap_canary_reproduction_qc.csv")
)


# ------------------------------------------------------------------------------
# 4. Final combined-content QC
# ------------------------------------------------------------------------------

all_results <- bind_rows(combined_result_list)

scenario_counts <- all_results %>%
  count(sensitivity_scenario_id, name = "n_rows")

successful_all <- all_results %>%
  filter(status == "ok")

production_qc <- data.frame(
  check = c(
    "Four completion-manifest rows are present",
    "Every scenario is complete",
    "Exactly 8,000 capped-IPW result rows are present",
    "Each scenario contains 2,000 repetitions",
    "Only the capped-IPW method label is present",
    "All results retain the primary MAR/MVN DGM labels",
    "Every successful result uses 1st/99th type-7 capping",
    "Every successful result caps at least one weight in each tail",
    "Every successful result obeys its capping bounds",
    "Production repetitions 1-2 reproduce the canary exactly"
  ),
  passed = c(
    nrow(manifest) == 4L,
    ipw_cap_all_true(manifest$complete),
    nrow(all_results) == 4L * nsim_ipw_cap,
    nrow(scenario_counts) == 4L &
      ipw_cap_all_true(scenario_counts$n_rows == nsim_ipw_cap),
    identical(unique(all_results$method), ipw_capped_label),
    ipw_cap_all_true(all_results$selection == "mar") &
      ipw_cap_all_true(all_results$marker_error == "mvn"),
    nrow(successful_all) > 0L &
      ipw_cap_all_true(
        successful_all$cap_lower_probability == ipw_cap_lower_probability
      ) &
      ipw_cap_all_true(
        successful_all$cap_upper_probability == ipw_cap_upper_probability
      ) &
      ipw_cap_all_true(
        successful_all$cap_quantile_type == ipw_cap_quantile_type
      ),
    nrow(successful_all) > 0L &
      ipw_cap_all_true(successful_all$n_capped_lower > 0L) &
      ipw_cap_all_true(successful_all$n_capped_upper > 0L),
    nrow(successful_all) > 0L &
      ipw_cap_all_true(
        successful_all$capped_weight_min >=
          successful_all$cap_lower_threshold - 1e-12
      ) &
      ipw_cap_all_true(
        successful_all$capped_weight_max <=
          successful_all$cap_upper_threshold + 1e-12
      ),
    nrow(canary_reproduction_qc) == nrow(ipw_cap_canary_grid) &
      ipw_cap_all_true(canary_reproduction_qc$canary_rows == 2L) &
      ipw_cap_all_true(canary_reproduction_qc$production_rows == 2L) &
      ipw_cap_all_true(
        canary_reproduction_qc$reproduced_within_tolerance
      )
  ),
  stringsAsFactors = FALSE
)

production_qc$overall_production_content_pass <-
  ipw_cap_all_true(production_qc$passed)

ipw_cap_atomic_write_csv(
  production_qc,
  file.path(ipw_cap_combined_dir, "ipw_weight_capping_production_qc.csv")
)

ipw_cap_atomic_save_rds(
  list(
    scenarios = ipw_cap_sensitivity_grid,
    manifest = manifest,
    failures = failures,
    canary_reproduction_qc = canary_reproduction_qc,
    production_qc = production_qc
  ),
  file.path(
    ipw_cap_combined_dir,
    "ipw_weight_capping_combination_complete.rds"
  ),
  compress = "xz"
)

if (!ipw_cap_all_true(production_qc$passed)) {

  failed_checks <- production_qc$check[!production_qc$passed]

  stop(
    "IPW-capping production combination FAILED. Failed checks:\n- ",
    paste(failed_checks, collapse = "\n- ")
  )
}

ipw_cap_atomic_write_lines(
  c(
    "IPW WEIGHT CAPPING PRODUCTION COMBINATION: PASS",
    paste("Completed:", Sys.time()),
    paste("Scenarios:", nrow(ipw_cap_sensitivity_grid)),
    paste("Repetitions per scenario:", nsim_ipw_cap),
    paste("Total result rows:", nrow(all_results)),
    paste("Recorded failed rows:", nrow(failures)),
    "Canary reproduction: PASS"
  ),
  combine_pass_file
)

cat(
  "\nIPW weight-capping production combination PASSED.\n",
  "Combined directory: ",
  ipw_cap_combined_dir,
  "\n",
  sep = ""
)
