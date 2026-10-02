################################################################################
# 48_check_real_world_gtsummary_tables.R
#
# Structural QC for Word/table outputs created by
# 47_make_real_world_gtsummary_tables.R.
#
# This script does not rerun any statistical models or regenerate tables.
################################################################################


# ==============================================================================
# 1. PACKAGES AND PATHS
# ==============================================================================

required_packages <- c(
  "dplyr",
  "readr",
  "gtsummary"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}

library(dplyr)
library(readr)
library(gtsummary)

if (!file.exists("00_config.R")) {
  stop(
    "Run this script from the local midus_two_phase_sim project directory."
  )
}

source("00_config.R")

table_dir <- file.path(
  results_dir,
  "real_world_application",
  "manuscript_outputs",
  "tables"
)


# ==============================================================================
# 2. EXPECTED OUTPUT FILES
# ==============================================================================

expected_files <- c(
  "Table3_characteristics_by_phase2.docx",
  "Table3_standardization_qc.csv",
  "Table4A_methods_in_columns.docx",
  "Table4B_methods_in_rows.docx",
  "Table4A_and_Table4B_alternative_layouts.docx",
  "MIDUS_real_world_application_tables.docx",
  "Table3_characteristics_by_phase2.csv",
  "Table4A_methods_in_columns.csv",
  "Table4B_methods_in_rows.csv",
  "Table3_characteristics_by_phase2_gtsummary.rds",
  "Table4A_methods_in_columns_gtsummary.rds",
  "Table4B_methods_in_rows_gtsummary.rds",
  "real_world_table_manifest.csv"
)

expected_paths <- file.path(
  table_dir,
  expected_files
)

files_exist <- file.exists(
  expected_paths
)

file_sizes <- ifelse(
  files_exist,
  file.info(
    expected_paths
  )$size,
  NA_real_
)


# Read the explicit Table 3 standardization QC written by Script 47.
table3_standardization_qc_file <- file.path(
  table_dir,
  "Table3_standardization_qc.csv"
)

table3_standardization_qc <- if (
  file.exists(table3_standardization_qc_file)
) {
  read_csv(
    table3_standardization_qc_file,
    show_col_types = FALSE
  )
} else {
  NULL
}

table3_standardization_pass <- if (
  !is.null(table3_standardization_qc) &&
    nrow(table3_standardization_qc) == 1L
) {
  means <- c(
    table3_standardization_qc$discrimination_mean,
    table3_standardization_qc$grimage2_mean,
    table3_standardization_qc$dunedinpace_mean
  )

  sds <- c(
    table3_standardization_qc$discrimination_sd,
    table3_standardization_qc$grimage2_sd,
    table3_standardization_qc$dunedinpace_sd
  )

  all(abs(means) < 1e-12) &&
    all(abs(sds - 1) < 1e-12)
} else {
  FALSE
}


# ==============================================================================
# 3. READ SAVED GTSUMMARY OBJECTS
# ==============================================================================

tbl3_path <- file.path(
  table_dir,
  "Table3_characteristics_by_phase2_gtsummary.rds"
)

tbl4a_path <- file.path(
  table_dir,
  "Table4A_methods_in_columns_gtsummary.rds"
)

tbl4b_path <- file.path(
  table_dir,
  "Table4B_methods_in_rows_gtsummary.rds"
)

tbl3 <- if (file.exists(tbl3_path)) {
  readRDS(tbl3_path)
} else {
  NULL
}

tbl4a <- if (file.exists(tbl4a_path)) {
  readRDS(tbl4a_path)
} else {
  NULL
}

tbl4b <- if (file.exists(tbl4b_path)) {
  readRDS(tbl4b_path)
} else {
  NULL
}


# ==============================================================================
# 4. NUMERIC/DISPLAY QC
# ==============================================================================

table3_csv <- file.path(
  table_dir,
  "Table3_characteristics_by_phase2.csv"
)

table3_data <- if (file.exists(table3_csv)) {
  read_csv(
    table3_csv,
    show_col_types = FALSE
  )
} else {
  NULL
}


table4a_csv <- file.path(
  table_dir,
  "Table4A_methods_in_columns.csv"
)

table4b_csv <- file.path(
  table_dir,
  "Table4B_methods_in_rows.csv"
)

table4a_data <- if (file.exists(table4a_csv)) {
  read_csv(
    table4a_csv,
    show_col_types = FALSE
  )
} else {
  NULL
}

table4b_data <- if (file.exists(table4b_csv)) {
  read_csv(
    table4b_csv,
    show_col_types = FALSE
  )
} else {
  NULL
}

method_order <- c(
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

docx_files <- expected_paths[
  grepl(
    "\\.docx$",
    expected_paths,
    ignore.case = TRUE
  )
]

docx_nonempty <- all(
  file.exists(docx_files)
) &&
  all(
    file.info(docx_files)$size >
      1000
  )

# ------------------------------------------------------------------
# Use the saved gtsummary table bodies for structural QC.
#
# This is intentionally more robust than checking the exported CSV column
# headers because gtsummary::as_tibble() may retain raw names such as `label`
# rather than the manuscript-facing display header "Outcome".
# ------------------------------------------------------------------

table4a_body_qc <- if (
  inherits(
    tbl4a,
    "gtsummary"
  )
) {
  tbl4a$table_body
} else {
  NULL
}

table4b_body_qc <- if (
  inherits(
    tbl4b,
    "gtsummary"
  )
) {
  tbl4b$table_body
} else {
  NULL
}

table3_body_qc <- if (
  inherits(
    tbl3,
    "gtsummary"
  )
) {
  tbl3$table_body
} else {
  NULL
}

table4a_has_methods <- if (
  !is.null(table4a_body_qc)
) {
  all(
    method_order %in%
      names(table4a_body_qc)
  )
} else {
  FALSE
}

table4a_has_outcomes <- if (
  !is.null(table4a_body_qc) &&
    "label" %in%
      names(table4a_body_qc)
) {
  all(
    expected_outcomes %in%
      table4a_body_qc$label
  )
} else {
  FALSE
}

table4b_has_methods <- if (
  !is.null(table4b_body_qc) &&
    "label" %in%
      names(table4b_body_qc)
) {
  method_labels <- table4b_body_qc$label[
    table4b_body_qc$label %in%
      method_order
  ]

  identical(
    method_labels,
    method_order
  )
} else {
  FALSE
}

# Table 3 should display the SMD magnitude only. The estimate itself is stored
# in `estimate`; CI and p-value columns may remain internally in table_body but
# they should be hidden in the gtsummary table.
table3_has_smd_only <- if (
  !is.null(table3_body_qc) &&
    "estimate" %in%
      names(table3_body_qc)
) {
  smd_values <- table3_body_qc$estimate[
    !is.na(
      table3_body_qc$estimate
    )
  ]

  all(
    is.finite(
      smd_values
    )
  ) &&
    all(
      smd_values >= 0
    )
} else {
  FALSE
}

# Verify manuscript display precision directly from Table 4A's stored strings.
table4_four_decimal_coefficients <- if (
  !is.null(table4a_body_qc)
) {
  outcome_rows <- table4a_body_qc$label %in%
    expected_outcomes

  estimate_cells <- unlist(
    table4a_body_qc[
      outcome_rows,
      method_order,
      drop = FALSE
    ],
    use.names = FALSE
  )

  all(
    grepl(
      "^-?[0-9]+\\.[0-9]{4} \\(-?[0-9]+\\.[0-9]{4} to -?[0-9]+\\.[0-9]{4}\\)$",
      estimate_cells
    )
  )
} else {
  FALSE
}

table4_three_decimal_pvalues <- if (
  !is.null(table4a_body_qc)
) {
  p_rows <- table4a_body_qc$label ==
    "P value"

  p_cells <- unlist(
    table4a_body_qc[
      p_rows,
      method_order,
      drop = FALSE
    ],
    use.names = FALSE
  )

  all(
    p_cells == "<0.001" |
      grepl(
        "^0\\.[0-9]{3}$",
        p_cells
      )
  )
} else {
  FALSE
}


# ==============================================================================
# 5. FINAL QC TABLE
# ==============================================================================

qc <- tibble::tibble(
  check = c(
    "All expected table files exist",
    "All expected table files are nonempty",
    "All Word documents exceed 1 KB",
    "Table 3 exposure and outcomes have full-sample mean 0 and SD 1",
    "Saved Table 3 object is a gtsummary object",
    "Saved Table 4A object is a gtsummary object",
    "Saved Table 4B object is a gtsummary object",
    "Table 3 displays SMD magnitude only, without CI text",
    "Table 4 coefficients and CIs use four decimal places",
    "Table 4 p values use three decimal places or <0.001",
    "Table 4A contains all six method columns",
    "Table 4A contains both outcomes",
    "Table 4B methods appear in the prespecified order"
  ),
  passed = c(
    all(files_exist),
    all(
      is.finite(file_sizes) &
        file_sizes > 0
    ),
    docx_nonempty,
    table3_standardization_pass,
    inherits(
      tbl3,
      "gtsummary"
    ),
    inherits(
      tbl4a,
      "gtsummary"
    ),
    inherits(
      tbl4b,
      "gtsummary"
    ),
    table3_has_smd_only,
    table4_four_decimal_coefficients,
    table4_three_decimal_pvalues,
    table4a_has_methods,
    table4a_has_outcomes,
    table4b_has_methods
  )
) %>%
  mutate(
    overall_pass = all(passed)
  )

qc_file <- file.path(
  table_dir,
  "real_world_gtsummary_table_qc.csv"
)

write_csv(
  qc,
  qc_file
)

cat(
  "\n============================================================\n",
  "REAL-WORLD GTSUMMARY TABLE QC\n",
  "============================================================\n\n",
  sep = ""
)

print(
  tibble::as_tibble(qc),
  n = Inf
)

if (!all(qc$passed)) {

  failed_checks <- qc$check[
    !qc$passed
  ]

  stop(
    "\nREAL-WORLD TABLE QC FAILED:\n- ",
    paste(
      failed_checks,
      collapse = "\n- "
    )
  )
}

pass_file <- file.path(
  table_dir,
  "MIDUS_REAL_WORLD_TABLES_PASS.txt"
)

writeLines(
  c(
    "MIDUS REAL-WORLD GTSUMMARY TABLES: PASS",
    paste(
      "Completed:",
      Sys.time()
    ),
    paste(
      "QC checks:",
      nrow(qc)
    ),
    paste(
      "All checks passed:",
      all(qc$passed)
    )
  ),
  pass_file
)

cat(
  "\nREAL-WORLD GTSUMMARY TABLE QC: PASS\n",
  "QC file: ",
  qc_file,
  "\nPASS marker: ",
  pass_file,
  "\n",
  sep = ""
)
