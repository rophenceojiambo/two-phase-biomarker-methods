################################################################################
# MIDUS REAL-WORLD APPLICATION: ANALYTIC SAMPLE
#
# Exposure:
#   Everyday discrimination
#
# Outcomes:
#   GrimAge2
#   DunedinPACE
#
# Phase-1 covariates:
#   Age
#   Sex
#   Race/ethnicity
#
# BMI is extracted from Project 4 for source-data provenance but is not
# included in the final analytic dataset or analysis models.
#
# Phase-2 covariates:
#   Eight RNA-based leukocyte-marker transcripts
#
# Plate variables are kept separate:
#   dnam_plate = BRA6DMPLATE
#   rna_plate  = RA6RPLATE
################################################################################


# ==============================================================================
# 1. PACKAGES AND DATA DIRECTORY
# ==============================================================================

library(haven)
library(dplyr)
library(tidyr)
library(readr)

# Directory containing the authorized MIDUS source files.
# For reproducibility, set MIDUS_SOURCE_DIR to the local restricted-data location.
# The default data/source directory is ignored by Git.
data_dir <- Sys.getenv(
  "MIDUS_SOURCE_DIR",
  unset = file.path("data", "source")
)

# Analysis-ready files are written to data/ by default so that scripts 44-48
# can find data/MIDUS_discrimination_analysis.rds directly.
# Set MIDUS_ANALYTIC_DIR to override this location.
output_dir <- Sys.getenv(
  "MIDUS_ANALYTIC_DIR",
  unset = "data"
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# 2. READ SOURCE DATA
# ==============================================================================

# ------------------------------------------------------------------------------
# DNA methylation age scores
#
# Purpose:
#   Provides the two epigenetic aging outcomes:
#     - GrimAge2
#     - DunedinPACE
#
#   Also provides the DNA methylation assay plate:
#     - BRA6DMPLATE
#
#   Participants with available methylation-age scores define the base
#   Phase-1 cohort.
# ------------------------------------------------------------------------------

dna <- read_sav(file.path(data_dir, "M2MR1_GEN_DNAmAge_N2118_20230822.sav"),
                user_na = FALSE)


# ------------------------------------------------------------------------------
# Project 4 biomarker data
#
# Purpose:
#   Provides BMI measured during the biomarker assessment:
#     - RA4PBMI
#
#   Project 1 age is used rather than age at the biomarker visit.
# ------------------------------------------------------------------------------

p4 <- read_sav(file.path(data_dir, "MR1_BIO_AGGREGATE_N863_20240405.sav"),
               user_na = FALSE)


# ------------------------------------------------------------------------------
# MIDUS Refresher Main Project 1 survey
#
# Purpose:
#   Provides Phase-1 variables for the Main Refresher sample:
#     - Project 1 age
#     - Sex
#     - Race/ethnicity
#     - Everyday discrimination
#
#   Everyday discrimination:
#     - RA1SDAYDI
# ------------------------------------------------------------------------------

p1 <- read_sav(file.path(data_dir, "MR1_P1_SURVEY_N3577_20250407.sav"),
               user_na = FALSE)


# ------------------------------------------------------------------------------
# Milwaukee Refresher 1 Project 1 survey
#
# Purpose:
#   Provides corresponding Phase-1 variables for the Milwaukee Refresher:
#     - Project 1 age
#     - Sex
#     - Race/ethnicity
#     - Everyday discrimination
#
#   Everyday discrimination:
#     - RAACDAYDI
#
#   These variables are harmonized with the Main Refresher variables.
# ------------------------------------------------------------------------------

mke <- read_sav(file.path(data_dir, "36722-0001-Data-REST.sav"),
                user_na = FALSE)


# ------------------------------------------------------------------------------
# MIDUS Refresher RNA / gene-expression data
#
# Purpose:
#   Provides the eight partially observed leukocyte-marker transcripts:
#     - CD19
#     - CD3D
#     - CD3E
#     - CD4
#     - CD8A
#     - CD14
#     - FCGR3A
#     - NCAM1
#
#   Also provides RNA assay plate:
#     - RA6RPLATE
#
#   Availability of all eight markers defines Phase 2.
# ------------------------------------------------------------------------------

rna <- read_sav(file.path(data_dir, "MR1_P6_RNA_SCORE_N863_20210923.sav"),
                user_na = FALSE)


# ------------------------------------------------------------------------------
# RNA file read again with SPSS user-defined missing codes retained
#
# Purpose:
#   Used only to inspect the original coding of RA6RPLATE, including:
#     - valid assay plate numbers
#     - 98 = MISSING
#     - 99 = INAPP
#
#   This object is diagnostic only.
# ------------------------------------------------------------------------------

rna_codes <- read_sav(file.path(data_dir, "MR1_P6_RNA_SCORE_N863_20210923.sav"),
                      user_na = TRUE)


# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

# Remove SPSS value labels while retaining numeric values.
num <- function(x) as.numeric(zap_labels(x))


# Everyday discrimination has a valid constructed range of 9-36.
clean_discrimination <- function(x) {
  x <- num(x)
  x[x < 9 | x > 36] <- NA_real_
  x
}


# Safeguard for RNA special missing codes.
clean_rna <- function(x) {
  x <- num(x)
  x[x %in% c(98, 99)] <- NA_real_
  x
}


# ==============================================================================
# 4. DNA METHYLATION OUTCOMES
# ==============================================================================

# Keep MIDUS Refresher participants with available methylation-age scores.
# This creates the base Phase-1 cohort.

dnam <- dna %>%
  transmute(
    id = num(M2MRID),
    sample = num(SAMPLMAJ),
    mr_case = num(M2MRCASE),
    dnam_available = num(BRA6DMAVAIL),
    grimage2 = num(BRA6DMAGEGRIMAGE2),
    dunedinpace = num(BRA6DMAGEDUNEDINPACE),
    dnam_plate = num(BRA6DMPLATE)
  ) %>%
  filter(mr_case == 2, sample %in% c(20, 21), dnam_available == 1) %>%
  select(id, sample, grimage2, dunedinpace, dnam_plate)


# ==============================================================================
# 5. PROJECT 4 BMI
# ==============================================================================

# BMI is taken from the biomarker examination.
# Age is intentionally not taken from Project 4.

p4_vars <- p4 %>%
  transmute(
    id = num(MRID),
    sample = num(SAMPLMAJ),
    bmi = num(RA4PBMI)
  ) %>%
  filter(sample %in% c(20, 21))


# ==============================================================================
# 6. MAIN REFRESHER PROJECT 1 VARIABLES
# ==============================================================================

# Main Refresher Phase-1 variables:
#   age            = Project 1 age
#   sex            = Project 1 sex
#   race_eth       = harmonized race/ethnicity
#   discrimination = everyday discrimination

p1_main <- p1 %>%
  transmute(
    id = num(MRID),
    sample = num(SAMPLMAJ),
    age = num(RA1PRAGE),
    sex_code = num(RA1PRSEX),
    hispanic = num(RA1PF1),
    race = num(RA1PF7A),
    discrimination = clean_discrimination(RA1SDAYDI)
  ) %>%
  filter(sample == 20) %>%
  mutate(
    sex = factor(sex_code, levels = c(1, 2), labels = c("Male", "Female")),
    race_eth = case_when(
      hispanic %in% 2:7 ~ "Other",
      hispanic == 1 & race == 1 ~ "non-Hispanic White",
      hispanic == 1 & race == 2 ~ "non-Hispanic Black",
      hispanic == 1 & race %in% 3:6 ~ "Other",
      TRUE ~ NA_character_
    )
  ) %>%
  select(id, sample, age, sex, race_eth, discrimination)


# ==============================================================================
# 7. MILWAUKEE REFRESHER PROJECT 1 VARIABLES
# ==============================================================================

# Milwaukee Refresher uses corresponding RAAC variables.
#
# The Milwaukee Refresher sample was designed as an African-American sample.
# Non-Hispanic participants are therefore coded as non-Hispanic Black;
# Hispanic participants are included in the "Other" category.

p1_mke <- mke %>%
  transmute(
    id = num(MRID),
    sample = num(SAMPLMAJ),
    age = num(RAACRAGE),
    sex_code = num(RAACRSEX),
    hispanic = num(RAACF1),
    discrimination = clean_discrimination(RAACDAYDI)
  ) %>%
  filter(sample == 21) %>%
  mutate(
    sex = factor(sex_code, levels = c(1, 2), labels = c("Male", "Female")),
    race_eth = case_when(
      hispanic == 1 ~ "non-Hispanic Black",
      hispanic %in% 2:7 ~ "Other",
      TRUE ~ NA_character_
    )
  ) %>%
  select(id, sample, age, sex, race_eth, discrimination)


# Combine Main and Milwaukee Project 1 variables.
survey <- bind_rows(p1_main, p1_mke) %>%
  mutate(
    race_eth = factor(
      race_eth,
      levels = c("non-Hispanic White", "non-Hispanic Black", "Other")
    )
  )


# ==============================================================================
# 8. RNA LEUKOCYTE-MARKER VARIABLES
# ==============================================================================

# The RNA variables are already log2-transformed transcript abundances.

rna_vars <- rna %>%
  transmute(
    id = num(MRID),
    sample = num(SAMPLMAJ),
    rna_plate = num(RA6RPLATE),
    log2_cd19 = clean_rna(RA6RCD19),
    log2_cd3d = clean_rna(RA6RCD3D),
    log2_cd3e = clean_rna(RA6RCD3E),
    log2_cd4 = clean_rna(RA6RCD4),
    log2_cd8a = clean_rna(RA6RCD8A),
    log2_cd14 = clean_rna(RA6RCD14),
    log2_fcgr3a = clean_rna(RA6RFCGR3A),
    log2_ncam1 = clean_rna(RA6RNCAM1)
  ) %>%
  filter(sample %in% c(20, 21))


marker_vars <- c(
  "log2_cd19", "log2_cd3d", "log2_cd3e", "log2_cd4",
  "log2_cd8a", "log2_cd14", "log2_fcgr3a", "log2_ncam1"
)


# ==============================================================================
# 9. MERGE ANALYTIC DATA
# ==============================================================================

# Start from participants with available methylation-age scores.
# Add BMI, Project 1 variables, and the partially observed RNA markers.

df_analytic <- dnam %>%
  left_join(p4_vars, by = c("id", "sample")) %>%
  left_join(survey, by = c("id", "sample")) %>%
  left_join(rna_vars, by = c("id", "sample")) %>%
  mutate(
    race_eth = factor(
      race_eth,
      levels = c("non-Hispanic White", "non-Hispanic Black", "Other")
    ),
    discrimination_STD = as.numeric(scale(discrimination)),
    n_rna_observed = rowSums(!is.na(across(all_of(marker_vars)))),
    phase2 = as.integer(n_rna_observed == 8)
  ) %>%
  select(id, grimage2, dunedinpace, dnam_plate, rna_plate, 
         age, sex, race_eth, discrimination,
         discrimination_STD, all_of(marker_vars), phase2)


# ==============================================================================
# 10. PHASE-1 VARIABLE MISSINGNESS
# ==============================================================================

# RNA markers are intentionally excluded from this table because their
# missingness is the two-phase problem of interest.

phase1_missingness <- df_analytic %>%
  summarise(
    across(
      c(grimage2, dunedinpace, discrimination, age, sex, race_eth),
      ~ sum(is.na(.x))
    )
  ) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(percent_missing = round(100 * n_missing / nrow(df_analytic), 2))

print(phase1_missingness)


# ==============================================================================
# 11. CREATE FINAL DISCRIMINATION ANALYSIS SAMPLE
# ==============================================================================

# The methodological problem concerns missing RNA covariates.
#
# Therefore, the exposure, outcomes, and Phase-1 covariates are required
# to be observed. RNA markers are intentionally NOT included in this filter.
#
# Continuous variables are standardized after defining the final analytic
# sample so that all scaling constants correspond to the same target sample.

df_discrimination <- df_analytic %>%
  filter(
    !is.na(discrimination),
    !is.na(grimage2),
    !is.na(dunedinpace),
    !is.na(age),
    !is.na(sex),
    !is.na(race_eth)
  )


# Helper function for z-standardization.
#
# na.rm = TRUE is important for the RNA markers because they are observed
# only in the Phase-2 subset. Missing RNA values remain NA.

standardize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}


# Standardize all continuous variables used to calibrate the simulation.
df_discrimination <- df_discrimination %>%
  mutate(
    # Primary exposure
    discrimination_STD = standardize(discrimination),
    
    # Epigenetic aging outcomes
    grimage2_STD = standardize(grimage2),
    dunedinpace_STD = standardize(dunedinpace),
    
    # Continuous Phase-1 covariates
    age_STD = standardize(age),
    
    # Eight partially observed RNA markers
    across(
      all_of(marker_vars),
      standardize,
      .names = "{.col}_STD"
    )
  )


# ==============================================================================
# 12. SAVE DATASETS
# ==============================================================================

# Base two-phase analytic dataset.
saveRDS(df_analytic, file.path(output_dir, "MIDUS_two_phase_analytic.rds"))
write_csv(df_analytic, file.path(output_dir, "MIDUS_two_phase_analytic.csv"))


# Complete Phase-1 dataset for the discrimination application.
saveRDS(df_discrimination, file.path(output_dir, "MIDUS_discrimination_analysis.rds"))
write_csv(df_discrimination, file.path(output_dir, "MIDUS_discrimination_analysis.csv"))


# ==============================================================================
# 13. FINAL DATASET STRUCTURE
# ==============================================================================
glimpse(df_analytic)
glimpse(df_discrimination)
