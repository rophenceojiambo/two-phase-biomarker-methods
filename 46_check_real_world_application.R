################################################################################
# 46_check_real_world_application.R
#
# Final structural and numerical QC for the completed MIDUS real-world
# application. This script does not rerun any statistical model.
################################################################################

library(dplyr)
library(readr)
library(tidyr)

project_dir <- normalizePath(".", winslash = "/", mustWork = FALSE)

real_world_dir <- file.path(
  project_dir,
  "results",
  "real_world_application"
)

results_file <- file.path(
  real_world_dir,
  "midus_real_world_method_results.csv"
)

qc_file <- file.path(
  real_world_dir,
  "midus_real_world_result_qc.csv"
)

sample_file <- file.path(
  real_world_dir,
  "midus_real_world_sample_summary.csv"
)

marker_file <- file.path(
  real_world_dir,
  "midus_real_world_marker_missingness.csv"
)

weight_file <- file.path(
  real_world_dir,
  "midus_real_world_weight_diagnostics.csv"
)

required_files <- c(
  results_file,
  qc_file,
  sample_file,
  marker_file,
  weight_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required real-world output file(s) not found:\n",
    paste0("  - ", missing_files, collapse = "\n")
  )
}

results <- read_csv(
  results_file,
  show_col_types = FALSE
)

result_qc <- read_csv(
  qc_file,
  show_col_types = FALSE
)

sample_summary <- read_csv(
  sample_file,
  show_col_types = FALSE
)

marker_missingness <- read_csv(
  marker_file,
  show_col_types = FALSE
)

weights <- read_csv(
  weight_file,
  show_col_types = FALSE
)

expected_methods <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

expected_outcomes <- c(
  "GrimAge2",
  "DunedinPACE"
)

expected_keys <- tidyr::expand_grid(
  outcome = expected_outcomes,
  method = expected_methods
)

observed_keys <- results %>%
  count(
    outcome,
    method,
    name = "n"
  )

key_join <- expected_keys %>%
  left_join(
    observed_keys,
    by = c("outcome", "method")
  )

sample_n1 <- sample_summary$N_phase1[[1]]
sample_n2 <- sample_summary$N_phase2[[1]]
sample_fraction <- sample_summary$phase2_fraction[[1]]
sample_not2 <- sample_summary$N_not_phase2[[1]]
sample_partial <- sample_summary$N_partial_marker_rows[[1]]

weight_checks <- weights %>%
  summarise(
    all_finite = all(
      is.finite(weight_min) &
        is.finite(weight_p99) &
        is.finite(weight_max) &
        is.finite(weight_cv) &
        is.finite(weight_ess)
    ),
    all_positive = all(
      weight_min > 0 &
        weight_p99 > 0 &
        weight_max > 0
    ),
    p99_not_above_max = all(
      weight_p99 <= weight_max
    ),
    ess_valid = all(
      weight_ess > 0 &
        weight_ess <= sample_n2
    ),
    max_weight_below_10 = all(
      weight_max < 10
    )
  )

marker_checks <- marker_missingness %>%
  summarise(
    eight_markers = n() == 8L,
    same_missing_count = n_distinct(n_missing) == 1L,
    missing_equals_nonphase2 = all(
      n_missing == sample_not2
    )
  )

qc_checks <- tibble::tibble(
  check = c(
    "Exactly two outcomes",
    "Exactly six methods",
    "Exactly 12 outcome-method rows",
    "One row per expected outcome-method combination",
    "All method fits have status ok",
    "All estimates are finite",
    "All standard errors are finite and positive",
    "All confidence limits are finite",
    "Saved result-QC file passes for all 12 fits",
    "Phase-1 N is positive",
    "Phase-2 N is positive and no larger than Phase-1 N",
    "Saved Phase-2 fraction equals N2/N1",
    "No partially observed RNA-marker blocks",
    "Exactly eight RNA markers are represented",
    "All eight markers have the same missing count",
    "RNA missing count equals number outside Phase 2",
    "All weight diagnostics are finite",
    "All weights are positive",
    "Weight 99th percentile does not exceed maximum",
    "Weight ESS is positive and no larger than Phase-2 N",
    "Maximum observed weight is below 10"
  ),
  passed = c(
    setequal(unique(results$outcome), expected_outcomes),
    setequal(unique(results$method), expected_methods),
    nrow(results) == 12L,
    nrow(key_join) == 12L &&
      all(!is.na(key_join$n)) &&
      all(key_join$n == 1L),
    all(results$status == "ok"),
    all(is.finite(results$estimate)),
    all(is.finite(results$se) & results$se > 0),
    all(
      is.finite(results$conf_low) &
        is.finite(results$conf_high)
    ),
    nrow(result_qc) == 12L &&
      all(result_qc$passed),
    is.finite(sample_n1) &&
      sample_n1 > 0,
    is.finite(sample_n2) &&
      sample_n2 > 0 &&
      sample_n2 <= sample_n1,
    isTRUE(
      all.equal(
        sample_fraction,
        sample_n2 / sample_n1,
        tolerance = 1e-12
      )
    ),
    sample_partial == 0L,
    marker_checks$eight_markers[[1]],
    marker_checks$same_missing_count[[1]],
    marker_checks$missing_equals_nonphase2[[1]],
    weight_checks$all_finite[[1]],
    weight_checks$all_positive[[1]],
    weight_checks$p99_not_above_max[[1]],
    weight_checks$ess_valid[[1]],
    weight_checks$max_weight_below_10[[1]]
  )
)

qc_checks <- qc_checks %>%
  mutate(
    overall_pass = all(passed)
  )

output_file <- file.path(
  real_world_dir,
  "midus_real_world_final_qc.csv"
)

write_csv(
  qc_checks,
  output_file
)

cat(
  "\n============================================================\n",
  "MIDUS REAL-WORLD APPLICATION — FINAL QC\n",
  "============================================================\n\n",
  sep = ""
)

print(
  tibble::as_tibble(qc_checks),
  n = Inf
)

cat(
  "\nPhase-1 N: ", sample_n1,
  "\nPhase-2 N: ", sample_n2,
  "\nPhase-2 fraction: ",
  round(sample_fraction, 4),
  "\nMaximum observed IPW/AIPW weight: ",
  round(max(weights$weight_max), 4),
  "\nMinimum observed weight ESS fraction: ",
  round(min(weights$weight_ess / sample_n2), 4),
  "\n\n",
  sep = ""
)

if (!all(qc_checks$passed)) {
  failed_checks <- qc_checks$check[!qc_checks$passed]

  stop(
    "FINAL REAL-WORLD QC FAILED:\n- ",
    paste(failed_checks, collapse = "\n- ")
  )
}

pass_file <- file.path(
  real_world_dir,
  "MIDUS_REAL_WORLD_APPLICATION_PASS.txt"
)

writeLines(
  c(
    "MIDUS REAL-WORLD APPLICATION: PASS",
    paste("Completed:", Sys.time()),
    paste("Phase-1 N:", sample_n1),
    paste("Phase-2 N:", sample_n2),
    paste(
      "Phase-2 fraction:",
      round(sample_fraction, 6)
    ),
    paste(
      "Maximum observed weight:",
      round(max(weights$weight_max), 6)
    ),
    paste(
      "Minimum weight ESS fraction:",
      round(min(weights$weight_ess / sample_n2), 6)
    )
  ),
  pass_file
)

cat(
  "FINAL REAL-WORLD QC: PASS\n",
  "QC file: ", output_file, "\n",
  "PASS marker: ", pass_file, "\n",
  sep = ""
)
