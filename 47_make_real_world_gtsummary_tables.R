################################################################################
# 47_make_real_world_gtsummary_tables.R
#
# Creates manuscript-ready Word tables for the MIDUS real-world application.
#
# Outputs:
#   1. Table 3: Characteristics by Phase-2 RNA-marker availability
#   2. Table 4: Lower-triangle Pearson correlation matrix in Phase 2
#   3. Table 4A: Alternative method-comparison layout with methods in columns
#   4. Table 4B: Alternative method-comparison layout with methods in rows
#   5. One Word document containing Table 4A and Table 4B together
#   6. One combined Word document containing all manuscript/application tables
#
# The script DOES NOT rerun any statistical model, multiple imputation,
# inverse probability weighting, or AIPW analysis.
#
# Formatting:
#   - gtsummary JAMA theme
#   - Arial in Word exports
#   - Regression coefficients and 95% CIs: 4 decimal places
#   - P values: 3 decimal places, with <0.001
#   - Table 3 reports absolute standardized mean differences (SMDs)
#   - Table 4 reports Pearson correlations to 2 decimals with significance stars
#   - Table 4 displays the strict lower triangle only
################################################################################


# ==============================================================================
# 1. REQUIRED PACKAGES
# ==============================================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "gtsummary",
  "flextable",
  "officer",
  "Hmisc",
  "smd"
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
    "Install the following package(s) before running Script 47:\n  ",
    paste(missing_packages, collapse = ", "),
    "\n\nRun:\ninstall.packages(c(",
    paste0('"', missing_packages, '"', collapse = ", "),
    "))"
  )
}

library(dplyr)
library(tidyr)
library(readr)
library(gtsummary)
library(flextable)
library(officer)


# ==============================================================================
# 2. PROJECT AND INPUT PATHS
# ==============================================================================

if (!file.exists("00_config.R")) {
  stop(
    "Run this script from the midus_two_phase_sim project directory."
  )
}

source("00_config.R")

# Real-world analysis results produced by Script 44.
real_world_dir <- file.path(
  results_dir,
  "real_world_application"
)

results_file <- file.path(
  real_world_dir,
  "midus_real_world_method_results.csv"
)

# Final MIDUS analytic dataset. The file is stored in the project data folder.
analytic_file <- file.path(
  "data",
  "MIDUS_discrimination_analysis.rds"
)

required_files <- c(
  analytic_file,
  results_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {
  stop(
    "Required input file(s) not found:\n",
    paste0(
      "  - ",
      missing_files,
      collapse = "\n"
    )
  )
}


# ==============================================================================
# 3. OUTPUT DIRECTORIES
# ==============================================================================

manuscript_dir <- file.path(
  real_world_dir,
  "manuscript_outputs"
)

table_dir <- file.path(
  manuscript_dir,
  "tables"
)

dir.create(
  table_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ==============================================================================
# 4. SET GTSUMMARY JAMA THEME
# ==============================================================================

gtsummary::reset_gtsummary_theme()

gtsummary::theme_gtsummary_journal(
  journal = "jama"
)

# Slightly tighter spacing is useful for the six-method comparison table.
gtsummary::theme_gtsummary_compact()


# ==============================================================================
# 5. READ AND VALIDATE INPUT DATA
# ==============================================================================

analytic_data <- readRDS(
  analytic_file
)

method_results <- read_csv(
  results_file,
  show_col_types = FALSE
)

raw_marker_names <- c(
  "log2_cd19",
  "log2_cd3d",
  "log2_cd3e",
  "log2_cd4",
  "log2_cd8a",
  "log2_cd14",
  "log2_fcgr3a",
  "log2_ncam1"
)

required_analytic_columns <- c(
  "age",
  "sex",
  "race_eth",
  "discrimination",
  "grimage2",
  "dunedinpace",
  "discrimination_STD",
  "grimage2_STD",
  "dunedinpace_STD",
  raw_marker_names,
  "phase2"
)

missing_analytic_columns <- setdiff(
  required_analytic_columns,
  names(analytic_data)
)

if (length(missing_analytic_columns) > 0L) {
  stop(
    "Analytic RDS is missing required variable(s): ",
    paste(
      missing_analytic_columns,
      collapse = ", "
    )
  )
}

required_result_columns <- c(
  "outcome",
  "method",
  "estimate",
  "conf_low",
  "conf_high",
  "p_value",
  "N_phase1",
  "N_phase2",
  "phase2_fraction",
  "status"
)

missing_result_columns <- setdiff(
  required_result_columns,
  names(method_results)
)

if (length(missing_result_columns) > 0L) {
  stop(
    "Real-world method-results file is missing required variable(s): ",
    paste(
      missing_result_columns,
      collapse = ", "
    )
  )
}

method_order <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

outcome_order <- c(
  "GrimAge2",
  "DunedinPACE"
)

if (
  !setequal(
    unique(method_results$method),
    method_order
  )
) {
  stop(
    "Method-results file does not contain exactly the six expected methods."
  )
}

if (
  !setequal(
    unique(method_results$outcome),
    outcome_order
  )
) {
  stop(
    "Method-results file does not contain exactly GrimAge2 and DunedinPACE."
  )
}

if (
  nrow(method_results) !=
    length(method_order) *
      length(outcome_order)
) {
  stop(
    "Expected 12 outcome-method rows in the real-world results."
  )
}

if (
  !all(method_results$status == "ok")
) {
  stop(
    "At least one real-world method fit did not have status = 'ok'."
  )
}


# ==============================================================================
# 6. FORMATTING HELPERS
# ==============================================================================

format_coef_ci <- function(
    estimate,
    conf_low,
    conf_high
) {
  sprintf(
    "%.4f (%.4f to %.4f)",
    estimate,
    conf_low,
    conf_high
  )
}

format_p_value <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

application_n_phase1 <- unique(
  method_results$N_phase1
)

application_n_phase2 <- unique(
  method_results$N_phase2
)

if (
  length(application_n_phase1) != 1L ||
    length(application_n_phase2) != 1L
) {
  stop(
    "Expected one common Phase-1 N and one common Phase-2 N across ",
    "the real-world method results."
  )
}

application_note <- paste0(
  "Note: Values are standardized regression coefficients (95% CIs), ",
  "representing the SD difference in the epigenetic-aging outcome associated ",
  "with a 1-SD higher everyday discrimination score. The Phase-1 analytic ",
  "sample included N = ", application_n_phase1,
  " participants; N = ", application_n_phase2,
  " had complete Phase-2 RNA-marker measurements. All approaches adjust for ",
  "age, sex, and race/ethnicity. The Naive analysis does not include the RNA ",
  "markers; CCA, FCS-MI, JM-MI, IPW, and AIPW incorporate the eight Phase-2 ",
  "RNA leukocyte-marker transcripts according to their respective procedures. ",
  "CCA indicates complete-case analysis; FCS-MI, fully conditional ",
  "specification multiple imputation; JM-MI, joint-model multiple imputation; ",
  "IPW, inverse probability weighting; AIPW, augmented inverse probability ",
  "weighting; CI, confidence interval."
)

phase2_note <- paste0(
  "Note: Phase 1 only denotes participants in the Phase-1 analytic sample ",
  "without complete measurements of all eight RNA leukocyte-marker transcripts; ",
  "Phase 2 denotes participants with complete measurements of all eight markers. ",
  "Everyday discrimination, GrimAge2, and DunedinPACE were standardized to ",
  "mean 0 and SD 1 in the full Phase-1 analytic sample. RNA-marker values are ",
  "shown on their original log2-normalized transcript-abundance scale for Phase 2. ",
  "Because the RNA markers are unobserved by design in the Phase 1 only group, ",
  "no Phase-1 summary or SMD is shown for those marker rows. SMD indicates the ",
  "absolute standardized mean difference."
)


# ==============================================================================
# 7. TABLE 3: CHARACTERISTICS BY PHASE-2 RNA-MARKER AVAILABILITY
# ==============================================================================

# The final analytic dataset already contains the standardized exposure and
# outcomes used in the real-world analysis. Use those variables directly so
# Table 3 reflects the exact analysis-ready data rather than recomputing them.

table3_data <- analytic_data %>%
  transmute(
    phase_group = factor(
      if_else(
        phase2 == 1L,
        "Phase 2",
        "Phase 1 only"
      ),
      levels = c(
        "Phase 1 only",
        "Phase 2"
      )
    ),
    age = age,
    sex = droplevels(factor(sex)),
    race_eth = droplevels(factor(race_eth)),
    discrimination_STD = discrimination_STD,
    grimage2_STD = grimage2_STD,
    dunedinpace_STD = dunedinpace_STD,
    log2_cd19 = log2_cd19,
    log2_cd3d = log2_cd3d,
    log2_cd3e = log2_cd3e,
    log2_cd4 = log2_cd4,
    log2_cd8a = log2_cd8a,
    log2_cd14 = log2_cd14,
    log2_fcgr3a = log2_fcgr3a,
    log2_ncam1 = log2_ncam1
  )

if (anyNA(table3_data$phase_group)) {
  stop(
    "phase2 contains missing or non-binary values in the analytic dataset."
  )
}

# QC the standardized variables already stored in the analytic dataset.
table3_standardization_qc <- table3_data %>%
  summarise(
    discrimination_mean = mean(discrimination_STD),
    discrimination_sd = sd(discrimination_STD),
    grimage2_mean = mean(grimage2_STD),
    grimage2_sd = sd(grimage2_STD),
    dunedinpace_mean = mean(dunedinpace_STD),
    dunedinpace_sd = sd(dunedinpace_STD)
  )

standardization_pass <- all(
  abs(
    c(
      table3_standardization_qc$discrimination_mean,
      table3_standardization_qc$grimage2_mean,
      table3_standardization_qc$dunedinpace_mean
    )
  ) < 1e-10
) &&
  all(
    abs(
      c(
        table3_standardization_qc$discrimination_sd,
        table3_standardization_qc$grimage2_sd,
        table3_standardization_qc$dunedinpace_sd
      ) - 1
    ) < 1e-10
  )

if (!standardization_pass) {
  stop(
    "Stored standardized variables failed QC: expected full-sample mean 0 ",
    "and SD 1 for discrimination_STD, grimage2_STD, and dunedinpace_STD."
  )
}

cat(
  "Stored standardized-variable QC: PASS\n",
  "  Everyday discrimination: mean = ",
  sprintf("%.6f", table3_standardization_qc$discrimination_mean),
  ", SD = ",
  sprintf("%.6f", table3_standardization_qc$discrimination_sd),
  "\n",
  "  GrimAge2: mean = ",
  sprintf("%.6f", table3_standardization_qc$grimage2_mean),
  ", SD = ",
  sprintf("%.6f", table3_standardization_qc$grimage2_sd),
  "\n",
  "  DunedinPACE: mean = ",
  sprintf("%.6f", table3_standardization_qc$dunedinpace_mean),
  ", SD = ",
  sprintf("%.6f", table3_standardization_qc$dunedinpace_sd),
  "\n\n",
  sep = ""
)

n_phase1_only <- sum(
  table3_data$phase_group ==
    "Phase 1 only"
)

n_phase2 <- sum(
  table3_data$phase_group ==
    "Phase 2"
)

n_total <- nrow(
  table3_data
)

tbl3 <- table3_data %>%
  tbl_summary(
    by = phase_group,
    include = c(
      age,
      sex,
      race_eth,
      discrimination_STD,
      grimage2_STD,
      dunedinpace_STD,
      all_of(raw_marker_names)
    ),
    label = list(
      age ~ "Age, years",
      sex ~ "Sex",
      race_eth ~ "Race/ethnicity",
      discrimination_STD ~ "Everyday discrimination, standardized",
      grimage2_STD ~ "GrimAge2, standardized",
      dunedinpace_STD ~ "DunedinPACE, standardized",
      log2_cd19 ~ "CD19",
      log2_cd3d ~ "CD3D",
      log2_cd3e ~ "CD3E",
      log2_cd4 ~ "CD4",
      log2_cd8a ~ "CD8A",
      log2_cd14 ~ "CD14",
      log2_fcgr3a ~ "FCGR3A",
      log2_ncam1 ~ "NCAM1"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p})"
    ),
    type = list(
      sex ~ "categorical",
      race_eth ~ "categorical"
    ),
    digits = list(
      age ~ 1,
      discrimination_STD ~ 2,
      grimage2_STD ~ 2,
      dunedinpace_STD ~ 2,
      all_of(raw_marker_names) ~ 2,
      all_categorical() ~ c(0, 1)
    ),
    missing = "no"
  ) %>%

  # SMDs are meaningful only for variables observed in both Phase-1 groups.
  # The RNA markers are observed only in Phase 2, so they are excluded here.
  add_difference(
    include = c(
      age,
      sex,
      race_eth,
      discrimination_STD,
      grimage2_STD,
      dunedinpace_STD
    ),
    test = everything() ~ "smd",
    estimate_fun =
      everything() ~
        gtsummary::label_style_number(
          digits = 2
        )
  ) %>%

  # Report only the magnitude of the SMD. No CI or p value is needed for this
  # descriptive balance table.
  gtsummary::remove_column_merge(
    columns = estimate
  ) %>%

  modify_table_body(
    ~ .x %>%
      mutate(
        estimate = if_else(
          is.na(estimate),
          estimate,
          abs(estimate)
        ),
        # RNA values are unavailable in Phase 1 only by design.
        # Leave that column blank rather than printing NA (NA).
        stat_1 = if_else(
          variable %in% raw_marker_names,
          "",
          stat_1
        )
      )
  ) %>%

  modify_column_hide(
    columns = dplyr::any_of(
      c(
        "conf.low",
        "conf.high",
        "p.value"
      )
    )
  ) %>%

  modify_header(
    label ~ "**Characteristic**",
    all_stat_cols() ~ "**{level}, N = {n}**",
    estimate ~ "**SMD**"
  ) %>%

  gtsummary::remove_footnote_header(
    columns = everything()
  ) %>%

  gtsummary::remove_abbreviation() %>%

  gtsummary::modify_source_note(
    phase2_note
  ) %>%

  modify_caption(
    paste0(
      "Table 3. Characteristics of the MIDUS analytic sample ",
      "by Phase-2 RNA-marker availability"
    )
  ) %>%

  bold_labels()


# ==============================================================================
# 8. TABLE 4: PHASE-2 PEARSON CORRELATION MATRIX
# ==============================================================================

# Correlations among Phase-2 participants for:
#   GrimAge2, DunedinPACE, age, daily discrimination, and the eight
#   original log2 RNA leukocyte-marker transcript-abundance variables.
#
# GrimAge2, DunedinPACE, and daily discrimination are shown as z-scores to
# match the scale used in the real-world regression analyses. Age remains in
# years for interpretability.

table4_corr_data <- analytic_data %>%
  filter(
    phase2 == 1L
  ) %>%
  transmute(
    `GrimAge2 z-score` = grimage2_STD,
    `DunedinPACE z-score` = dunedinpace_STD,
    Age = age,
    `Daily discrimination z-score` = discrimination_STD,
    CD19 = log2_cd19,
    CD3D = log2_cd3d,
    CD3E = log2_cd3e,
    CD4 = log2_cd4,
    CD8A = log2_cd8a,
    CD14 = log2_cd14,
    FCGR3A = log2_fcgr3a,
    NCAM1 = log2_ncam1
  )

if (
  nrow(table4_corr_data) !=
    as.integer(application_n_phase2)
) {
  stop(
    "Correlation-table Phase-2 N does not match the real-world analysis N."
  )
}

if (anyNA(table4_corr_data)) {
  missing_by_variable <- vapply(
    table4_corr_data,
    function(x) sum(is.na(x)),
    integer(1)
  )

  stop(
    "Unexpected missing values remain among the 12 variables in the Phase-2 ",
    "correlation table:\n",
    paste(
      names(missing_by_variable),
      missing_by_variable,
      sep = " = ",
      collapse = "\n"
    )
  )
}

if (
  !all(
    vapply(
      table4_corr_data,
      is.numeric,
      logical(1)
    )
  )
) {
  stop(
    "All variables in the Phase-2 correlation table must be numeric."
  )
}

cor_results <- Hmisc::rcorr(
  as.matrix(
    table4_corr_data
  ),
  type = "pearson"
)

r_mat <- cor_results$r
p_mat <- cor_results$P

# Significance stars:
# * p < 0.05; ** p < 0.01; *** p < 0.001
stars_mat <- matrix(
  "",
  nrow = nrow(p_mat),
  ncol = ncol(p_mat),
  dimnames = dimnames(p_mat)
)

stars_mat[
  !is.na(p_mat) &
    p_mat < 0.05
] <- "*"

stars_mat[
  !is.na(p_mat) &
    p_mat < 0.01
] <- "**"

stars_mat[
  !is.na(p_mat) &
    p_mat < 0.001
] <- "***"

# Display only the strict lower triangle. The diagonal and upper triangle are
# left blank because they contain redundant information.
corr_text_mat <- matrix(
  "",
  nrow = nrow(r_mat),
  ncol = ncol(r_mat),
  dimnames = dimnames(r_mat)
)

for (i in seq_len(nrow(r_mat))) {
  for (j in seq_len(ncol(r_mat))) {

    if (i > j) {

      # Avoid displaying negative zero after rounding.
      r_value <- r_mat[i, j]

      if (
        is.finite(r_value) &&
          abs(r_value) < 0.005
      ) {
        r_value <- 0
      }

      corr_text_mat[i, j] <- paste0(
        sprintf(
          "%.2f",
          r_value
        ),
        stars_mat[i, j]
      )
    }
  }
}

cor_var_labels <- colnames(
  table4_corr_data
)

n_cor_vars <- length(
  cor_var_labels
)

# The final (12th) numeric column would contain only a blank diagonal entry and
# blank upper-triangle entries, so only columns 1--11 are displayed.
corr_table_final <- as.data.frame(
  corr_text_mat[
    ,
    seq_len(n_cor_vars - 1L),
    drop = FALSE
  ],
  check.names = FALSE
) %>%
  mutate(
    Variable = paste0(
      seq_len(n_cor_vars),
      ". ",
      cor_var_labels
    ),
    .before = 1
  )

names(corr_table_final) <- c(
  "Variable",
  as.character(
    seq_len(n_cor_vars - 1L)
  )
)

correlation_note <- paste0(
  "Note: Pearson correlations among Phase-2 participants (N = ",
  nrow(table4_corr_data),
  "). GrimAge2, DunedinPACE, and daily discrimination are standardized ",
  "z-scores; age is reported in years. CD19, CD3D, CD3E, CD4, CD8A, CD14, ",
  "FCGR3A, and NCAM1 are log2-transformed normalized transcript-abundance ",
  "values. *p < 0.05; **p < 0.01; ***p < 0.001."
)

cor_header_labels <- c(
  Variable = "",
  stats::setNames(
    as.character(
      seq_len(n_cor_vars - 1L)
    ),
    as.character(
      seq_len(n_cor_vars - 1L)
    )
  )
)

ft4corr <- flextable::flextable(
  corr_table_final
) %>%

  flextable::set_header_labels(
    values = cor_header_labels
  ) %>%

  flextable::theme_booktabs() %>%

  # Bold variable labels and the numbered column headers.
  flextable::bold(
    j = "Variable",
    part = "body"
  ) %>%

  flextable::bold(
    part = "header"
  ) %>%

  flextable::align(
    j = "Variable",
    align = "left",
    part = "all"
  ) %>%

  flextable::align(
    j = 2:ncol(corr_table_final),
    align = "center",
    part = "all"
  ) %>%

  flextable::font(
    fontname = "Arial",
    part = "all"
  ) %>%

  flextable::fontsize(
    size = 8,
    part = "all"
  ) %>%

  flextable::padding(
    padding.top = 2,
    padding.bottom = 2,
    padding.left = 3,
    padding.right = 3,
    part = "all"
  ) %>%

  flextable::width(
    j = "Variable",
    width = 2.55
  ) %>%

  flextable::width(
    j = 2:ncol(corr_table_final),
    width = 0.57
  ) %>%

  flextable::set_table_properties(
    layout = "fixed",
    width = 1
  ) %>%

  flextable::add_footer_lines(
    values = correlation_note
  ) %>%

  flextable::align(
    align = "left",
    part = "footer"
  ) %>%

  flextable::fontsize(
    size = 8,
    part = "footer"
  ) %>%

  flextable::set_caption(
    caption = paste0(
      "Table 4. Pearson correlations among Phase-2 characteristics and ",
      "RNA leukocyte-marker transcripts in the MIDUS Refresher analytic sample"
    )
  )


# ==============================================================================
# 9. PREPARE MODEL-COMPARISON VALUES
# ==============================================================================

model_display <- method_results %>%
  mutate(
    method = factor(
      method,
      levels = method_order
    ),
    outcome = factor(
      outcome,
      levels = outcome_order
    ),
    estimate_ci = format_coef_ci(
      estimate,
      conf_low,
      conf_high
    ),
    p_display = format_p_value(
      p_value
    )
  ) %>%
  arrange(
    outcome,
    method
  )


# ==============================================================================
# 10. TABLE 4A: METHODS IN COLUMNS
# ==============================================================================

build_table4a_rows <- function(
    data,
    outcome_name
) {

  dd <- data %>%
    filter(
      as.character(outcome) ==
        outcome_name
    ) %>%
    arrange(method)

  estimate_values <- setNames(
    dd$estimate_ci,
    as.character(dd$method)
  )

  p_values <- setNames(
    dd$p_display,
    as.character(dd$method)
  )

  estimate_row <- tibble::tibble(
    variable = outcome_name,
    row_type = "label",
    label = outcome_name
  )

  p_row <- tibble::tibble(
    variable = outcome_name,
    row_type = "level",
    label = "P value"
  )

  for (method_name in method_order) {
    estimate_row[[method_name]] <-
      estimate_values[[method_name]]

    p_row[[method_name]] <-
      p_values[[method_name]]
  }

  bind_rows(
    estimate_row,
    p_row
  )
}

table4a_body <- bind_rows(
  lapply(
    outcome_order,
    function(outcome_name) {
      build_table4a_rows(
        model_display,
        outcome_name
      )
    }
  )
)

tbl4a <- gtsummary::as_gtsummary(
  table4a_body
) %>%

  modify_column_hide(
    columns = c(
      variable,
      row_type
    )
  ) %>%

  modify_header(
    label ~ "**Outcome**",
    Naive ~ "**Naive**",
    CCA ~ "**CCA**",
    `FCS-MI` ~ "**FCS-MI**",
    `JM-MI` ~ "**JM-MI**",
    IPW ~ "**IPW**",
    AIPW ~ "**AIPW**"
  ) %>%

  modify_column_alignment(
    columns = all_of(
      method_order
    ),
    align = "center"
  ) %>%

  modify_italic(
    columns = label,
    rows = row_type == "level"
  ) %>%

  bold_labels() %>%

  gtsummary::remove_footnote_header(
    columns = everything()
  ) %>%

  gtsummary::remove_abbreviation() %>%

  gtsummary::modify_source_note(
    application_note
  ) %>%

  modify_caption(
    paste0(
      "Table 4A. Comparison of statistical methods for estimating ",
      "associations between everyday discrimination and epigenetic aging ",
      "in MIDUS"
    )
  )


# ==============================================================================
# 11. TABLE 4B: METHODS IN ROWS
# ==============================================================================

table4b_body <- model_display %>%
  select(
    outcome,
    method,
    estimate_ci,
    p_display
  ) %>%

  pivot_wider(
    names_from = outcome,
    values_from = c(
      estimate_ci,
      p_display
    )
  ) %>%

  mutate(
    method = factor(
      method,
      levels = method_order
    )
  ) %>%

  arrange(method) %>%

  transmute(
    variable = "method",
    row_type = "level",
    label = as.character(method),
    grim_estimate =
      estimate_ci_GrimAge2,
    grim_p =
      p_display_GrimAge2,
    pace_estimate =
      estimate_ci_DunedinPACE,
    pace_p =
      p_display_DunedinPACE
  )

tbl4b <- gtsummary::as_gtsummary(
  table4b_body
) %>%

  modify_column_hide(
    columns = c(
      variable,
      row_type
    )
  ) %>%

  modify_header(
    label ~ "**Method**",
    grim_estimate ~ "**Coefficient (95% CI)**",
    grim_p ~ "**P value**",
    pace_estimate ~ "**Coefficient (95% CI)**",
    pace_p ~ "**P value**"
  ) %>%

  modify_spanning_header(
    c(
      grim_estimate,
      grim_p
    ) ~ "**GrimAge2**"
  ) %>%

  modify_spanning_header(
    c(
      pace_estimate,
      pace_p
    ) ~ "**DunedinPACE**"
  ) %>%

  modify_column_alignment(
    columns = c(
      grim_estimate,
      grim_p,
      pace_estimate,
      pace_p
    ),
    align = "center"
  ) %>%

  gtsummary::remove_footnote_header(
    columns = everything()
  ) %>%

  gtsummary::remove_abbreviation() %>%

  gtsummary::modify_source_note(
    application_note
  ) %>%

  modify_caption(
    paste0(
      "Table 4B. Comparison of statistical methods for estimating ",
      "associations between everyday discrimination and epigenetic aging ",
      "in MIDUS"
    )
  )


# ==============================================================================
# 12. SAVE TABLE OBJECTS AND MACHINE-READABLE DATA
# ==============================================================================

# Save Table 3 standardization QC explicitly.
write_csv(
  table3_standardization_qc,
  file.path(
    table_dir,
    "Table3_standardization_qc.csv"
  )
)

# Save the formatted Table 4 body plus the underlying correlation and p-value
# matrices for reproducibility and manuscript cross-checking.
write_csv(
  corr_table_final,
  file.path(
    table_dir,
    "Table4_phase2_correlations_display.csv"
  )
)

write_csv(
  as.data.frame(
    r_mat,
    check.names = FALSE
  ) %>%
    mutate(
      Variable = rownames(r_mat),
      .before = 1
    ),
  file.path(
    table_dir,
    "Table4_phase2_correlations_r.csv"
  )
)

write_csv(
  as.data.frame(
    p_mat,
    check.names = FALSE
  ) %>%
    mutate(
      Variable = rownames(p_mat),
      .before = 1
    ),
  file.path(
    table_dir,
    "Table4_phase2_correlations_p.csv"
  )
)

saveRDS(
  tbl3,
  file.path(
    table_dir,
    "Table3_characteristics_by_phase2_gtsummary.rds"
  )
)

saveRDS(
  tbl4a,
  file.path(
    table_dir,
    "Table4A_methods_in_columns_gtsummary.rds"
  )
)

saveRDS(
  tbl4b,
  file.path(
    table_dir,
    "Table4B_methods_in_rows_gtsummary.rds"
  )
)

write_csv(
  gtsummary::as_tibble(
    tbl3,
    col_labels = FALSE
  ),
  file.path(
    table_dir,
    "Table3_characteristics_by_phase2.csv"
  )
)

write_csv(
  gtsummary::as_tibble(
    tbl4a,
    col_labels = FALSE
  ),
  file.path(
    table_dir,
    "Table4A_methods_in_columns.csv"
  )
)

write_csv(
  gtsummary::as_tibble(
    tbl4b,
    col_labels = FALSE
  ),
  file.path(
    table_dir,
    "Table4B_methods_in_rows.csv"
  )
)


# ==============================================================================
# 13. CONVERT GTSUMMARY TABLES TO FLEXTABLE FOR WORD
# ==============================================================================

style_word_table <- function(
    tbl,
    font_size = 9.5
) {

  ft <- gtsummary::as_flex_table(
    tbl
  )

  ft <- ft %>%
    flextable::font(
      fontname = "Arial",
      part = "all"
    ) %>%
    flextable::fontsize(
      size = font_size,
      part = "all"
    ) %>%
    flextable::padding(
      padding.top = 3,
      padding.bottom = 3,
      padding.left = 4,
      padding.right = 4,
      part = "all"
    ) %>%
    flextable::set_table_properties(
      layout = "autofit",
      width = 1
    ) %>%
    flextable::autofit()

  ft
}

ft3 <- style_word_table(
  tbl3,
  font_size = 9.5
)

ft4a <- style_word_table(
  tbl4a,
  font_size = 8.5
)

ft4b <- style_word_table(
  tbl4b,
  font_size = 9.0
)


# ==============================================================================
# 14. WORD PAGE SETTINGS
# ==============================================================================

portrait_section <- officer::prop_section(
  page_size = officer::page_size(
    orient = "portrait",
    width = 8.5,
    height = 11
  ),
  page_margins = officer::page_mar(
    top = 0.6,
    bottom = 0.6,
    left = 0.65,
    right = 0.65
  ),
  type = "continuous"
)

landscape_section <- officer::prop_section(
  page_size = officer::page_size(
    orient = "landscape",
    width = 8.5,
    height = 11
  ),
  page_margins = officer::page_mar(
    top = 0.55,
    bottom = 0.55,
    left = 0.55,
    right = 0.55
  ),
  type = "continuous"
)


# ==============================================================================
# 15. EXPORT INDIVIDUAL WORD DOCUMENTS
# ==============================================================================

table3_docx <- file.path(
  table_dir,
  "Table3_characteristics_by_phase2.docx"
)

table4_corr_docx <- file.path(
  table_dir,
  "Table4_phase2_correlations.docx"
)

table4a_docx <- file.path(
  table_dir,
  "Table4A_methods_in_columns.docx"
)

table4b_docx <- file.path(
  table_dir,
  "Table4B_methods_in_rows.docx"
)

flextable::save_as_docx(
  ft3,
  path = table3_docx,
  pr_section = portrait_section,
  align = "center"
)

flextable::save_as_docx(
  ft4corr,
  path = table4_corr_docx,
  pr_section = landscape_section,
  align = "center"
)

flextable::save_as_docx(
  ft4a,
  path = table4a_docx,
  pr_section = landscape_section,
  align = "center"
)

flextable::save_as_docx(
  ft4b,
  path = table4b_docx,
  pr_section = portrait_section,
  align = "center"
)


# ==============================================================================
# 16. EXPORT TABLE 4A AND TABLE 4B TOGETHER
#
# This document is intended for side-by-side review of the two alternative
# manuscript layouts. Landscape orientation is used because Table 4A requires
# the additional horizontal space.
# ==============================================================================

table4_comparison_docx <- file.path(
  table_dir,
  "Table4A_and_Table4B_alternative_layouts.docx"
)

flextable::save_as_docx(
  values = list(
    "Table 4A. Methods in columns" =
      ft4a,
    "Table 4B. Methods in rows" =
      ft4b
  ),
  path = table4_comparison_docx,
  pr_section = landscape_section,
  align = "center"
)


# ==============================================================================
# 17. EXPORT ONE COMBINED WORD DOCUMENT
#
# All tables are placed in a landscape document so Table 4A fits comfortably.
# Individual portrait versions of Tables 3 and 4B are saved above.
# ==============================================================================

combined_docx <- file.path(
  table_dir,
  "MIDUS_real_world_application_tables.docx"
)

flextable::save_as_docx(
  values = list(
    "Table 3. Characteristics by Phase-2 RNA-marker availability" =
      ft3,
    "Table 4. Phase-2 Pearson correlation matrix" =
      ft4corr,
    "Table 4A. Method comparison — methods in columns" =
      ft4a,
    "Table 4B. Method comparison — methods in rows" =
      ft4b
  ),
  path = combined_docx,
  pr_section = landscape_section,
  align = "center"
)


# ==============================================================================
# 18. TABLE MANIFEST
# ==============================================================================

table_manifest <- tibble::tibble(
  table = c(
    "Table 3",
    "Table 4",
    "Table 4A",
    "Table 4B"
  ),
  purpose = c(
    "Characteristics by Phase-2 RNA-marker availability with absolute SMD",
    "Lower-triangle Pearson correlations among 12 Phase-2 variables",
    "Alternative real-world method comparison with methods in columns",
    "Alternative real-world method comparison with methods in rows"
  ),
  word_file = c(
    basename(table3_docx),
    basename(table4_corr_docx),
    basename(table4a_docx),
    basename(table4b_docx)
  ),
  orientation = c(
    "portrait",
    "landscape",
    "landscape",
    "portrait"
  )
)

write_csv(
  table_manifest,
  file.path(
    table_dir,
    "real_world_table_manifest.csv"
  )
)


# ==============================================================================
# 19. FINAL CONSOLE OUTPUT
# ==============================================================================

cat(
  "\n============================================================\n",
  "REAL-WORLD GTSUMMARY TABLES COMPLETE\n",
  "============================================================\n",
  "Analytic sample N: ", n_total, "\n",
  "Phase 1 only N: ", n_phase1_only, "\n",
  "Phase 2 N: ", n_phase2, "\n\n",
  "Word tables:\n",
  "  ", table3_docx, "\n",
  "  ", table4_corr_docx, "\n",
  "  ", table4a_docx, "\n",
  "  ", table4b_docx, "\n\n",
  "Table 4A + 4B comparison document:\n",
  "  ", table4_comparison_docx, "\n\n",
  "Combined all-table Word document:\n",
  "  ", combined_docx, "\n\n",
  "Table directory:\n",
  "  ", table_dir, "\n",
  "============================================================\n",
  sep = ""
)

# Restore package defaults for the remainder of the R session.
gtsummary::reset_gtsummary_theme()
