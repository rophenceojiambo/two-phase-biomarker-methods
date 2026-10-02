################################################################################
# 08_make_manuscript_outputs.R
#
# Produces ADEMP tables and manuscript/presentation-ready figures from the
# rsimsum summary outputs. No simulation is rerun.
################################################################################

source("00_config.R")

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

performance_tidy <- read_csv(
  file.path(
    summary_dir,
    "primary_rsimsum_tidy.csv"
  ),
  show_col_types = FALSE
)

performance_wide <- read_csv(
  file.path(
    summary_dir,
    "primary_performance_wide.csv"
  ),
  show_col_types = FALSE
)

type1 <- read_csv(
  file.path(
    summary_dir,
    "primary_type1_error.csv"
  ),
  show_col_types = FALSE
)

failures <- read_csv(
  file.path(
    summary_dir,
    "primary_failure_rates.csv"
  ),
  show_col_types = FALSE
)

runtime <- read_csv(
  file.path(
    summary_dir,
    "primary_runtime.csv"
  ),
  show_col_types = FALSE
)

phase2_diag <- read_csv(
  file.path(
    summary_dir,
    "primary_phase2_diagnostics.csv"
  ),
  show_col_types = FALSE
)

weights <- read_csv(
  file.path(
    summary_dir,
    "primary_weight_diagnostics.csv"
  ),
  show_col_types = FALSE
)

method_order <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

# ------------------------------------------------------------------------------
# 1. ADEMP design table
# ------------------------------------------------------------------------------

ademp_table <- tibble::tribble(
  ~component, ~item, ~specification,
  "Aims", "Primary question",
  "Compare six approaches for estimating the adjusted exposure coefficient when eight correlated RNA-based leukocyte-marker transcripts are observed only in a nested Phase-2 subsample.",
  "Data-generating mechanism", "Phase-1 sample size N",
  "500, 800, 1500",
  "Data-generating mechanism", "Expected Phase-2 fraction",
  "0.25, 0.50, 0.65",
  "Data-generating mechanism", "Exposure-marker partial R2",
  "0.01, 0.05, 0.10",
  "Data-generating mechanism", "Outcome-marker partial R2",
  "0, 0.02, 0.05, 0.10",
  "Data-generating mechanism", "True exposure coefficient theta",
  "0, 0.15, 0.30",
  "Data-generating mechanism", "Primary marker residual distribution",
  "Multivariate normal calibrated to MIDUS",
  "Data-generating mechanism", "Phase-2 availability",
  "MAR given A, Y, age, sex, race/ethnicity; intercept calibrated to target fraction",
  "Estimand", "Target",
  "Exposure coefficient theta from the full-data linear regression adjusted for age, sex, race/ethnicity, and the eight transcript markers",
  "Methods", "Compared methods",
  "Naive Phase-1, CCA, FCS-MI, JM-MI, IPW, AIPW",
  "Performance", "Primary measures",
  "Bias, empirical SE, model-based SE, SE ratio, 95% coverage, bias-eliminated coverage, Type I error when theta=0, failure rate",
  "Performance", "Additional diagnostics",
  "MSE, relative precision, runtime, Phase-2 sample size, selection probabilities, IPW/AIPW weight diagnostics",
  "Monte Carlo", "Repetitions",
  "2,000 per DGM; 324 primary DGMs; 648,000 simulated datasets"
)

write_csv(
  ademp_table,
  file.path(
    table_dir,
    "table_ADEMP_design.csv"
  )
)

# ------------------------------------------------------------------------------
# 2. Main scenario tables
# ------------------------------------------------------------------------------

scenario_labels <- primary_grid %>%
  mutate(
    scenario_label = case_when(
      N == 800L &
        abs(
          phase2_fraction -
            0.65
        ) < 1e-12 &
        abs(
          r2_a_marker -
            0.05
        ) < 1e-12 &
        abs(
          r2_y_marker -
            0.05
        ) < 1e-12 &
        abs(
          theta -
            0.15
        ) < 1e-12 ~
        "MIDUS-like",
      N == 500L &
        abs(
          phase2_fraction -
            0.25
        ) < 1e-12 &
        abs(
          r2_a_marker -
            0.10
        ) < 1e-12 &
        abs(
          r2_y_marker -
            0.10
        ) < 1e-12 &
        abs(
          theta -
            0.15
        ) < 1e-12 ~
        "Stress",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(
      scenario_label
    )
  ) %>%
  select(
    scenario_id,
    scenario_label
  )

main_table <- performance_wide %>%
  inner_join(
    scenario_labels,
    by = "scenario_id"
  ) %>%
  left_join(
    failures %>%
      select(
        scenario_id,
        method,
        failure_rate,
        mcse_failure
      ),
    by = c(
      "scenario_id",
      "method"
    )
  ) %>%
  left_join(
    runtime %>%
      select(
        scenario_id,
        method,
        mean_runtime_seconds
      ),
    by = c(
      "scenario_id",
      "method"
    )
  ) %>%
  transmute(
    scenario = scenario_label,
    method = factor(
      method,
      levels = method_order
    ),
    mean_estimate =
      est__thetamean,
    bias =
      est__bias,
    mcse_bias =
      mcse__bias,
    empirical_se =
      est__empse,
    mcse_empirical_se =
      mcse__empse,
    model_se =
      est__modelse,
    mcse_model_se =
      mcse__modelse,
    se_ratio =
      se_ratio,
    mcse_se_ratio =
      mcse_se_ratio,
    coverage =
      est__cover,
    mcse_coverage =
      mcse__cover,
    bias_eliminated_coverage =
      est__becover,
    mcse_bias_eliminated_coverage =
      mcse__becover,
    mse =
      est__mse,
    failure_rate =
      failure_rate,
    mean_runtime_seconds =
      mean_runtime_seconds
  ) %>%
  arrange(
    scenario,
    method
  )

write_csv(
  main_table,
  file.path(
    table_dir,
    "table_main_MIDUS_like_and_stress.csv"
  )
)

# ------------------------------------------------------------------------------
# 3. Null / Type I error table for MIDUS-like and stress-style design points
# ------------------------------------------------------------------------------

null_reference_ids <- primary_grid %>%
  filter(
    theta == 0,
    (
      N == 800L &
        phase2_fraction == 0.65 &
        r2_a_marker == 0.05 &
        r2_y_marker == 0.05
    ) |
      (
        N == 500L &
          phase2_fraction == 0.25 &
          r2_a_marker == 0.10 &
          r2_y_marker == 0.10
      )
  ) %>%
  mutate(
    scenario_label = ifelse(
      N == 800L,
      "MIDUS-like null",
      "Stress null"
    )
  ) %>%
  select(
    scenario_id,
    scenario_label
  )

type1_table <- type1 %>%
  inner_join(
    null_reference_ids,
    by = "scenario_id"
  ) %>%
  mutate(
    method = factor(
      method,
      levels = method_order
    )
  ) %>%
  arrange(
    scenario_label,
    method
  )

write_csv(
  type1_table,
  file.path(
    table_dir,
    "table_type1_reference_scenarios.csv"
  )
)

# ------------------------------------------------------------------------------
# 4. Calibration / Phase-2 diagnostic tables
# ------------------------------------------------------------------------------

write_csv(
  phase2_diag,
  file.path(
    table_dir,
    "table_phase2_diagnostics_all_scenarios.csv"
  )
)

write_csv(
  weights,
  file.path(
    table_dir,
    "table_weight_diagnostics_all_scenarios.csv"
  )
)

# ------------------------------------------------------------------------------
# 5. Deterioration flags
# ------------------------------------------------------------------------------

deterioration <- performance_wide %>%
  left_join(
    type1 %>%
      select(
        scenario_id,
        method,
        null_rejection_rate,
        mcse_null_rejection
      ),
    by = c(
      "scenario_id",
      "method"
    )
  ) %>%
  left_join(
    failures %>%
      select(
        scenario_id,
        method,
        failure_rate
      ),
    by = c(
      "scenario_id",
      "method"
    )
  ) %>%
  mutate(
    flag_bias =
      abs(
        est__bias
      ) >
        deterioration_thresholds$abs_bias,
    flag_se_ratio =
      se_ratio <
        deterioration_thresholds$se_ratio_low |
        se_ratio >
          deterioration_thresholds$se_ratio_high,
    flag_coverage =
      est__cover <
        deterioration_thresholds$coverage_low |
        est__cover >
          deterioration_thresholds$coverage_high,
    flag_type1 =
      ifelse(
        theta_true == 0,
        null_rejection_rate >
          deterioration_thresholds$type1_high,
        FALSE
      ),
    flag_failure =
      failure_rate >
        deterioration_thresholds$failure_rate,
    n_flags =
      rowSums(
        cbind(
          flag_bias,
          flag_se_ratio,
          flag_coverage,
          flag_type1,
          flag_failure
        ),
        na.rm = TRUE
      )
  )

write_csv(
  deterioration,
  file.path(
    table_dir,
    "table_deterioration_flags.csv"
  )
)

deterioration_by_method <- deterioration %>%
  group_by(
    method
  ) %>%
  summarise(
    scenarios = n(),
    pct_flag_bias =
      100 *
        mean(
          flag_bias,
          na.rm = TRUE
        ),
    pct_flag_se_ratio =
      100 *
        mean(
          flag_se_ratio,
          na.rm = TRUE
        ),
    pct_flag_coverage =
      100 *
        mean(
          flag_coverage,
          na.rm = TRUE
        ),
    pct_flag_type1 =
      100 *
        mean(
          flag_type1[
            theta_true == 0
          ],
          na.rm = TRUE
        ),
    pct_flag_failure =
      100 *
        mean(
          flag_failure,
          na.rm = TRUE
        ),
    .groups = "drop"
  )

write_csv(
  deterioration_by_method,
  file.path(
    table_dir,
    "table_deterioration_summary_by_method.csv"
  )
)

# ------------------------------------------------------------------------------
# 6. Plot helper
# ------------------------------------------------------------------------------

save_plot_both <- function(
  p,
  stem,
  width = 12,
  height = 8
) {

  ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        stem,
        ".pdf"
      )
    ),
    plot = p,
    width = width,
    height = height
  )

  ggsave(
    filename = file.path(
      figure_dir,
      paste0(
        stem,
        ".png"
      )
    ),
    plot = p,
    width = width,
    height = height,
    dpi = 300
  )
}

# ------------------------------------------------------------------------------
# 7. Bias figures
#
# One figure per theta x exposure-marker R2. This keeps the full factorial
# readable while allowing N, Phase-2 fraction, and outcome-marker R2 to remain
# visible in each figure.
# ------------------------------------------------------------------------------

bias_data <- performance_tidy %>%
  filter(
    stat == "bias"
  ) %>%
  mutate(
    method = factor(
      method,
      levels = method_order
    )
  )

for (
  theta_value in theta_levels
) {

  for (
    rA_value in r2_a_marker_levels
  ) {

    dd <- bias_data %>%
      filter(
        abs(
          theta_true -
            theta_value
        ) < 1e-12,
        abs(
          r2_a_marker -
            rA_value
        ) < 1e-12
      )

    p <- ggplot(
      dd,
      aes(
        x = r2_y_marker,
        y = est,
        group = method,
        linetype = method,
        shape = method
      )
    ) +
      geom_hline(
        yintercept = 0,
        linewidth = 0.4
      ) +
      geom_errorbar(
        aes(
          ymin =
            est -
              1.96 *
                mcse,
          ymax =
            est +
              1.96 *
                mcse
        ),
        width = 0.005,
        position = position_dodge(
          width = 0.006
        )
      ) +
      geom_line() +
      geom_point(
        size = 2
      ) +
      facet_grid(
        N ~ phase2_fraction,
        labeller = label_both
      ) +
      labs(
        x = "Outcome-marker partial R-squared",
        y = "Bias",
        linetype = "Method",
        shape = "Method",
        title = paste0(
          "Bias: theta = ",
          theta_value,
          ", exposure-marker partial R2 = ",
          rA_value
        )
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom"
      )

    save_plot_both(
      p,
      sprintf(
        "bias_theta_%03d_rA_%03d",
        round(
          100 *
            theta_value
        ),
        round(
          100 *
            rA_value
        )
      )
    )
  }
}

# ------------------------------------------------------------------------------
# 8. SE-ratio figures
# ------------------------------------------------------------------------------

se_ratio_data <- performance_tidy %>%
  filter(
    stat == "relerror"
  ) %>%
  mutate(
    se_ratio =
      1 +
        est /
          100,
    mcse_se_ratio =
      mcse /
        100,
    method = factor(
      method,
      levels = method_order
    )
  )

for (
  theta_value in theta_levels
) {

  for (
    rA_value in r2_a_marker_levels
  ) {

    dd <- se_ratio_data %>%
      filter(
        abs(
          theta_true -
            theta_value
        ) < 1e-12,
        abs(
          r2_a_marker -
            rA_value
        ) < 1e-12
      )

    p <- ggplot(
      dd,
      aes(
        x = r2_y_marker,
        y = se_ratio,
        group = method,
        linetype = method,
        shape = method
      )
    ) +
      geom_hline(
        yintercept = 1,
        linewidth = 0.4
      ) +
      geom_errorbar(
        aes(
          ymin =
            se_ratio -
              1.96 *
                mcse_se_ratio,
          ymax =
            se_ratio +
              1.96 *
                mcse_se_ratio
        ),
        width = 0.005,
        position = position_dodge(
          width = 0.006
        )
      ) +
      geom_line() +
      geom_point(
        size = 2
      ) +
      facet_grid(
        N ~ phase2_fraction,
        labeller = label_both
      ) +
      labs(
        x = "Outcome-marker partial R-squared",
        y = "Model SE / empirical SE",
        linetype = "Method",
        shape = "Method",
        title = paste0(
          "Standard-error calibration: theta = ",
          theta_value,
          ", exposure-marker partial R2 = ",
          rA_value
        )
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom"
      )

    save_plot_both(
      p,
      sprintf(
        "se_ratio_theta_%03d_rA_%03d",
        round(
          100 *
            theta_value
        ),
        round(
          100 *
            rA_value
        )
      )
    )
  }
}

# ------------------------------------------------------------------------------
# 9. Coverage + bias-eliminated coverage
# ------------------------------------------------------------------------------

coverage_data <- performance_tidy %>%
  filter(
    stat %in% c(
      "cover",
      "becover"
    )
  ) %>%
  mutate(
    coverage_type = recode(
      stat,
      cover = "Coverage",
      becover =
        "Bias-eliminated coverage"
    ),
    method = factor(
      method,
      levels = method_order
    )
  )

for (
  theta_value in theta_levels
) {

  for (
    rA_value in r2_a_marker_levels
  ) {

    dd <- coverage_data %>%
      filter(
        abs(
          theta_true -
            theta_value
        ) < 1e-12,
        abs(
          r2_a_marker -
            rA_value
        ) < 1e-12
      )

    p <- ggplot(
      dd,
      aes(
        x = r2_y_marker,
        y = est,
        group = interaction(
          method,
          coverage_type
        ),
        linetype = coverage_type,
        shape = method
      )
    ) +
      geom_hline(
        yintercept = nominal_level,
        linewidth = 0.4
      ) +
      geom_errorbar(
        aes(
          ymin =
            est -
              1.96 *
                mcse,
          ymax =
            est +
              1.96 *
                mcse
        ),
        width = 0.005
      ) +
      geom_line() +
      geom_point(
        size = 2
      ) +
      facet_grid(
        N ~ phase2_fraction,
        labeller = label_both
      ) +
      labs(
        x = "Outcome-marker partial R-squared",
        y = "Probability",
        linetype = "Coverage measure",
        shape = "Method",
        title = paste0(
          "95% interval performance: theta = ",
          theta_value,
          ", exposure-marker partial R2 = ",
          rA_value
        )
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom"
      )

    save_plot_both(
      p,
      sprintf(
        "coverage_theta_%03d_rA_%03d",
        round(
          100 *
            theta_value
        ),
        round(
          100 *
            rA_value
        )
      )
    )
  }
}

# ------------------------------------------------------------------------------
# 10. Type I error figures (theta = 0 only)
# ------------------------------------------------------------------------------

type1 <- type1 %>%
  mutate(
    method = factor(
      method,
      levels = method_order
    )
  )

for (
  rA_value in r2_a_marker_levels
) {

  dd <- type1 %>%
    filter(
      abs(
        r2_a_marker -
          rA_value
      ) < 1e-12
    )

  p <- ggplot(
    dd,
    aes(
      x = r2_y_marker,
      y = null_rejection_rate,
      group = method,
      linetype = method,
      shape = method
    )
  ) +
    geom_hline(
      yintercept = alpha_level,
      linewidth = 0.4
    ) +
    geom_errorbar(
      aes(
        ymin =
          null_rejection_rate -
            1.96 *
              mcse_null_rejection,
        ymax =
          null_rejection_rate +
            1.96 *
              mcse_null_rejection
      ),
      width = 0.005
    ) +
    geom_line() +
    geom_point(
      size = 2
    ) +
    facet_grid(
      N ~ phase2_fraction,
      labeller = label_both
    ) +
    labs(
      x = "Outcome-marker partial R-squared",
      y = "Null rejection rate",
      linetype = "Method",
      shape = "Method",
      title = paste0(
        "Null rejection rate: theta = 0, exposure-marker partial R2 = ",
        rA_value
      )
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom"
    )

  save_plot_both(
    p,
    sprintf(
      "type1_rA_%03d",
      round(
        100 *
          rA_value
      )
    )
  )
}

# ------------------------------------------------------------------------------
# 11. Empirical vs model-based SE, all scenarios
# ------------------------------------------------------------------------------

se_scatter <- performance_wide %>%
  mutate(
    method = factor(
      method,
      levels = method_order
    )
  )

p_se <- ggplot(
  se_scatter,
  aes(
    x = est__empse,
    y = est__modelse,
    shape = method
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linewidth = 0.4
  ) +
  geom_point(
    alpha = 0.6
  ) +
  facet_wrap(
    ~ method,
    scales = "free"
  ) +
  labs(
    x = "Empirical SE",
    y = "Model-based SE",
    title =
      "Empirical versus model-based standard errors across all primary DGMs"
  ) +
  theme_bw()

save_plot_both(
  p_se,
  "empirical_vs_model_se_all_scenarios",
  width = 12,
  height = 8
)

# ------------------------------------------------------------------------------
# 12. Failure-rate and runtime summaries
# ------------------------------------------------------------------------------

runtime_overview <- runtime %>%
  group_by(
    method,
    N
  ) %>%
  summarise(
    median_seconds =
      median(
        mean_runtime_seconds,
        na.rm = TRUE
      ),
    p95_seconds =
      unname(
        quantile(
          mean_runtime_seconds,
          0.95,
          na.rm = TRUE
        )
      ),
    .groups = "drop"
  )

write_csv(
  runtime_overview,
  file.path(
    table_dir,
    "table_runtime_overview.csv"
  )
)

p_runtime <- ggplot(
  runtime_overview,
  aes(
    x = factor(N),
    y = median_seconds,
    group = method,
    linetype = method,
    shape = method
  )
) +
  geom_line() +
  geom_point(
    size = 2
  ) +
  scale_y_log10() +
  labs(
    x = "Phase-1 sample size",
    y = "Median runtime per dataset (seconds, log scale)",
    linetype = "Method",
    shape = "Method",
    title = "Computational time by method and sample size"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )

save_plot_both(
  p_runtime,
  "runtime_by_method_and_N",
  width = 10,
  height = 7
)

failure_overview <- failures %>%
  group_by(
    method
  ) %>%
  summarise(
    mean_failure_rate =
      mean(
        failure_rate,
        na.rm = TRUE
      ),
    max_failure_rate =
      max(
        failure_rate,
        na.rm = TRUE
      ),
    scenarios_with_any_failure =
      sum(
        failure_rate > 0,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

write_csv(
  failure_overview,
  file.path(
    table_dir,
    "table_failure_overview.csv"
  )
)

# ------------------------------------------------------------------------------
# 13. Deterioration overview plot
# ------------------------------------------------------------------------------

deterioration_plot_data <- deterioration_by_method %>%
  pivot_longer(
    cols = starts_with(
      "pct_flag_"
    ),
    names_to = "criterion",
    values_to = "percent_scenarios"
  ) %>%
  mutate(
    criterion = recode(
      criterion,
      pct_flag_bias = "Bias",
      pct_flag_se_ratio = "SE ratio",
      pct_flag_coverage = "Coverage",
      pct_flag_type1 = "Type I error",
      pct_flag_failure = "Failure"
    ),
    method = factor(
      method,
      levels = method_order
    )
  )

p_det <- ggplot(
  deterioration_plot_data,
  aes(
    x = criterion,
    y = percent_scenarios,
    group = method,
    linetype = method,
    shape = method
  )
) +
  geom_line() +
  geom_point(
    size = 2
  ) +
  labs(
    x = NULL,
    y = "Primary DGMs crossing prespecified deterioration threshold (%)",
    linetype = "Method",
    shape = "Method",
    title = "Conditions in which method performance deteriorated"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

save_plot_both(
  p_det,
  "deterioration_summary",
  width = 10,
  height = 7
)

cat(
  "\nManuscript outputs complete.\n",
  "Tables: ",
  table_dir,
  "\nFigures: ",
  figure_dir,
  "\n",
  sep = ""
)
