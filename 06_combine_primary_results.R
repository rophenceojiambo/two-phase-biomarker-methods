################################################################################
# 06_combine_primary_results.R
#
# Validates all chunk files and combines them into one compact RDS file per DGM.
#
# Raw repetition-level estimates are preserved. This is the "estimates dataset"
# recommended for reproducible analysis of a simulation study.
################################################################################

source("00_config.R")

library(dplyr)
library(readr)

scenario_dir <- file.path(
  combined_dir,
  "scenarios"
)

rng_dir <- file.path(
  combined_dir,
  "rng_states"
)

dir.create(
  scenario_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  rng_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

expected_methods <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

manifest <- vector(
  "list",
  nrow(primary_grid)
)

failure_list <- list()
failure_index <- 1L

for (
  s in seq_len(
    nrow(primary_grid)
  )
) {

  scenario_id <- primary_grid$scenario_id[s]

  cat(
    "Combining scenario ",
    scenario_id,
    " / ",
    nrow(primary_grid),
    "\n",
    sep = ""
  )

  scenario_chunk_dir <- file.path(
    chunk_dir,
    sprintf(
      "scenario_%03d",
      scenario_id
    )
  )

  expected_files <- file.path(
    scenario_chunk_dir,
    sprintf(
      "scenario_%03d_chunk_%02d.rds",
      scenario_id,
      seq_len(
        chunks_per_scenario
      )
    )
  )

  missing_files <- expected_files[
    !file.exists(
      expected_files
    )
  ]

  if (length(missing_files) > 0) {

    manifest[[s]] <- data.frame(
      scenario_id = scenario_id,
      chunks_expected =
        chunks_per_scenario,
      chunks_found =
        chunks_per_scenario -
        length(missing_files),
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

    next
  }

  chunk_objects <- lapply(
    expected_files,
    readRDS
  )

  results <- bind_rows(
    lapply(
      chunk_objects,
      function(x) x$results
    )
  )

  rng_states <- bind_rows(
    lapply(
      chunk_objects,
      function(x) x$rng_states
    )
  )

  # ---------------------------------------------------------------------------
  # Validation
  # ---------------------------------------------------------------------------

  repetitions_found <- sort(
    unique(
      results$repetition
    )
  )

  expected_reps <- seq_len(
    nsim_primary
  )

  duplicate_keys <- results %>%
    count(
      scenario_id,
      repetition,
      method,
      name = "n"
    ) %>%
    filter(
      n != 1L
    )

  missing_methods <- setdiff(
    expected_methods,
    unique(
      results$method
    )
  )

  unexpected_methods <- setdiff(
    unique(
      results$method
    ),
    expected_methods
  )

  expected_rows <- nsim_primary *
    length(
      expected_methods
    )

  complete <- (
    setequal(
      repetitions_found,
      expected_reps
    ) &&
      nrow(
        duplicate_keys
      ) == 0L &&
      length(
        missing_methods
      ) == 0L &&
      length(
        unexpected_methods
      ) == 0L &&
      nrow(
        results
      ) == expected_rows
  )

  problems <- character()

  if (
    !setequal(
      repetitions_found,
      expected_reps
    )
  ) {
    problems <- c(
      problems,
      "repetition set incomplete"
    )
  }

  if (
    nrow(
      duplicate_keys
    ) > 0L
  ) {
    problems <- c(
      problems,
      "duplicate repetition-method keys"
    )
  }

  if (
    length(
      missing_methods
    ) > 0L
  ) {
    problems <- c(
      problems,
      paste0(
        "missing methods: ",
        paste(
          missing_methods,
          collapse = ", "
        )
      )
    )
  }

  if (
    length(
      unexpected_methods
    ) > 0L
  ) {
    problems <- c(
      problems,
      paste0(
        "unexpected methods: ",
        paste(
          unexpected_methods,
          collapse = ", "
        )
      )
    )
  }

  if (
    nrow(
      results
    ) != expected_rows
  ) {
    problems <- c(
      problems,
      paste0(
        "expected ",
        expected_rows,
        " rows, found ",
        nrow(results)
      )
    )
  }

  # ---------------------------------------------------------------------------
  # Save validated scenario estimate dataset
  # ---------------------------------------------------------------------------

  scenario_file <- file.path(
    scenario_dir,
    sprintf(
      "scenario_%03d_estimates.rds",
      scenario_id
    )
  )

  rng_file <- file.path(
    rng_dir,
    sprintf(
      "scenario_%03d_rng_states.rds",
      scenario_id
    )
  )

  if (complete) {

    saveRDS(
      results,
      scenario_file,
      compress = "xz"
    )

    saveRDS(
      rng_states,
      rng_file,
      compress = "xz"
    )
  }

  failed_rows <- results %>%
    filter(
      status != "ok" |
        !is.finite(
          estimate
        ) |
        !is.finite(
          se
        )
    )

  if (
    nrow(
      failed_rows
    ) > 0L
  ) {

    failure_list[[failure_index]] <- failed_rows %>%
      select(
        scenario_id,
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

    failure_index <- failure_index +
      1L
  }

  manifest[[s]] <- data.frame(
    scenario_id = scenario_id,
    chunks_expected =
      chunks_per_scenario,
    chunks_found =
      chunks_per_scenario,
    repetitions_found =
      length(
        repetitions_found
      ),
    result_rows =
      nrow(
        results
      ),
    complete = complete,
    problem = if (
      length(
        problems
      ) == 0L
    ) {
      ""
    } else {
      paste(
        problems,
        collapse = "; "
      )
    },
    stringsAsFactors = FALSE
  )
}

manifest <- bind_rows(
  manifest
)

write_csv(
  manifest,
  file.path(
    combined_dir,
    "primary_completion_manifest.csv"
  )
)

failures <- if (
  length(
    failure_list
  ) == 0L
) {
  tibble::tibble(
    scenario_id = integer(),
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
  bind_rows(
    failure_list
  )
}

write_csv(
  failures,
  file.path(
    combined_dir,
    "primary_method_failures.csv"
  )
)

n_incomplete <- sum(
  !manifest$complete
)

cat(
  "\nCombination check complete.\n",
  "Complete scenarios: ",
  sum(
    manifest$complete
  ),
  " / ",
  nrow(
    manifest
  ),
  "\nIncomplete scenarios: ",
  n_incomplete,
  "\nMethod-level failures recorded: ",
  nrow(
    failures
  ),
  "\n",
  sep = ""
)

if (n_incomplete > 0L) {

  cat(
    "\nPrimary simulation is NOT complete.\n",
    "See:\n",
    file.path(
      combined_dir,
      "primary_completion_manifest.csv"
    ),
    "\n"
  )

  quit(
    save = "no",
    status = 2L
  )
}
