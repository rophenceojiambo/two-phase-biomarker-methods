################################################################################
# 45_make_real_world_manuscript_outputs.R
#
# Creates manuscriptMtables and figures for the MIDUS real-world
# application using the already-completed six-method analysis from Script 44.
#
# This script uses the saved results from Script 44 and does not rerun multiple
# imputation, IPW, or AIPW.
################################################################################

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------------------------

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

weight_file <- file.path(
  real_world_dir,
  "midus_real_world_weight_diagnostics.csv"
)

sample_file <- file.path(
  real_world_dir,
  "midus_real_world_sample_summary.csv"
)

required_files <- c(
  results_file,
  weight_file,
  sample_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    "Required real-world output file(s) not found:\n",
    paste0("  - ", missing_files, collapse = "\n"),
    "\nRun 44_run_midus_real_world_application.R first."
  )
}

manuscript_dir <- file.path(
  real_world_dir,
  "manuscript_outputs"
)

figure_dir <- file.path(
  manuscript_dir,
  "figures"
)

table_dir <- file.path(
  manuscript_dir,
  "tables"
)

diagnostic_dir <- file.path(
  manuscript_dir,
  "diagnostics"
)

for (d in c(manuscript_dir, figure_dir, table_dir, diagnostic_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------------------------------------------------------
# 2. Read completed analysis outputs
# ------------------------------------------------------------------------------

results <- read_csv(
  results_file,
  show_col_types = FALSE
)

weights <- read_csv(
  weight_file,
  show_col_types = FALSE
)

sample_summary <- read_csv(
  sample_file,
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

outcome_order <- c(
  "GrimAge2",
  "DunedinPACE"
)

results <- results %>%
  mutate(
    method = factor(method, levels = method_order),
    outcome = factor(outcome, levels = outcome_order)
  )

weights <- weights %>%
  mutate(
    method = factor(method, levels = method_order),
    outcome = factor(outcome, levels = outcome_order)
  )

# ------------------------------------------------------------------------------
# 3. Method aesthetics
#
# These deliberately match the simulation figure aesthetics.
# ------------------------------------------------------------------------------

method_colors <- c(
  "Naive" = "#1B9E77",
  "CCA" = "#D95F02",
  "FCS-MI" = "#7570B3",
  "JM-MI" = "#E7298A",
  "IPW" = "#66A61E",
  "AIPW" = "#E6AB02"
)

method_shapes <- c(
  "Naive"  = 21,
  "CCA"    = 24,
  "FCS-MI" = 22,
  "JM-MI"  = 23,
  "IPW"    = 25,
  "AIPW"   = 10
)

method_linetypes <- c(
  "Naive"  = "solid",
  "CCA"    = "dashed",
  "FCS-MI" = "dotted",
  "JM-MI"  = "dotdash",
  "IPW"    = "longdash",
  "AIPW"   = "twodash"
)

figure_font <- if (.Platform$OS.type == "windows") {
  "Arial"
} else {
  "sans"
}

theme_real_world <- theme_bw(
  base_size = 12,
  base_family = figure_font
) +
  theme(
    text = element_text(
      family = figure_font,
      colour = "black"
    ),

    plot.title = element_blank(),

    strip.background = element_rect(
      fill = "grey94",
      color = "grey35",
      linewidth = 0.55
    ),

    strip.text.x = element_text(
      size = 10.4,
      colour = "black",
      lineheight = 1.05,
      margin = margin(
        t = 5,
        b = 5
      )
    ),

    strip.text.y = element_text(
      size = 10.8,
      colour = "black",
      margin = margin(
        l = 5,
        r = 5
      )
    ),

    panel.border = element_rect(
      color = "grey35",
      linewidth = 0.55
    ),

    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    axis.title = element_text(
      size = 12.2,
      colour = "black"
    ),

    axis.text = element_text(
      size = 10.5,
      colour = "black"
    ),

    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_blank(),

    legend.text = element_text(
      size = 10.2,
      colour = "black"
    ),

    legend.key.width = grid::unit(
      1.45,
      "lines"
    ),

    legend.spacing.x = grid::unit(
      0.20,
      "cm"
    ),

    legend.margin = margin(
      t = 2,
      b = 0
    ),

    plot.margin = margin(
      7,
      8,
      4,
      8
    )
  )

save_figure <- function(
    plot,
    stem,
    directory,
    width = 9.0,
    height = 5.6
) {

  ggsave(
    filename = file.path(
      directory,
      paste0(stem, ".png")
    ),
    plot = plot,
    width = width,
    height = height,
    dpi = 320,
    bg = "white"
  )

  ggsave(
    filename = file.path(
      directory,
      paste0(stem, ".pdf")
    ),
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    bg = "white"
  )
}

# ------------------------------------------------------------------------------
# 4. Manuscript numeric table
# ------------------------------------------------------------------------------

format_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(
      p < 0.001,
      "<0.001",
      sprintf("%.3f", p)
    )
  )
}

manuscript_table <- results %>%
  transmute(
    Outcome = as.character(outcome),
    Method = as.character(method),
    Estimate = round(estimate, 3),
    `Standard error` = round(se, 3),
    `95% CI` = sprintf(
      "%.3f, %.3f",
      conf_low,
      conf_high
    ),
    `P-value` = format_p(p_value),
    `Phase-1 N` = N_phase1,
    `Phase-2 N` = N_phase2,
    `Phase-2 fraction` = round(phase2_fraction, 3)
  )

write_csv(
  manuscript_table,
  file.path(
    table_dir,
    "Table_real_world_application.csv"
  )
)

# A second table keeps unrounded values.
write_csv(
  results %>%
    select(
      outcome,
      method,
      estimate,
      se,
      conf_low,
      conf_high,
      p_value,
      N_phase1,
      N_phase2,
      phase2_fraction
    ),
  file.path(
    table_dir,
    "Table_real_world_application_unrounded.csv"
  )
)

# ------------------------------------------------------------------------------
# 5A. Primary coefficient / confidence-interval figure
#
# This unannotated side-by-side version follows the simulation figure style.
# The annotated version below is used as Figure 6 in the manuscript.
# ------------------------------------------------------------------------------

plot_data_simple <- results %>%
  mutate(
    method_plot = factor(
      method,
      levels = rev(method_order)
    )
  )

base_plot_data_simple <- plot_data_simple %>%
  filter(
    method != "AIPW"
  )

aipw_plot_data_simple <- plot_data_simple %>%
  filter(
    method == "AIPW"
  )

coef_plot_simple <- ggplot(
  plot_data_simple,
  aes(
    x = estimate,
    y = method_plot,
    color = method
  )
) +

  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.65,
    color = "grey20"
  ) +

  geom_segment(
    aes(
      x = conf_low,
      xend = conf_high,
      y = method_plot,
      yend = method_plot
    ),
    linewidth = 0.85
  ) +

  geom_point(
    data = base_plot_data_simple,
    aes(
      shape = method,
      fill = method
    ),
    size = 3.4,
    stroke = 1.15
  ) +

  geom_point(
    data = aipw_plot_data_simple,
    shape = 21,
    color = "white",
    fill = "white",
    size = 4.5,
    stroke = 0,
    show.legend = FALSE
  ) +

  geom_point(
    data = aipw_plot_data_simple,
    aes(
      shape = method
    ),
    size = 3.6,
    stroke = 1.20,
    show.legend = TRUE
  ) +

  facet_wrap(
    ~ outcome,
    nrow = 1
  ) +

  scale_y_discrete(
    limits = rev(
      method_order
    )
  ) +

  scale_color_manual(
    values = method_colors,
    breaks = method_order,
    limits = method_order,
    drop = FALSE
  ) +

  scale_fill_manual(
    values = setNames(
      rep(
        "white",
        length(method_order)
      ),
      method_order
    ),
    breaks = method_order,
    limits = method_order,
    drop = FALSE
  ) +

  scale_shape_manual(
    values = method_shapes,
    breaks = method_order,
    limits = method_order,
    drop = FALSE
  ) +

  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      label.position = "right",
      override.aes = list(
        shape = unname(
          method_shapes
        ),
        fill = rep(
          "white",
          length(method_order)
        ),
        linewidth = 1.0,
        size = rep(
          3.8,
          length(method_order)
        )
      )
    ),
    fill = "none",
    shape = "none"
  ) +

  labs(
    x = paste0(
      "Standardized association with daily discrimination ",
      "(estimate and 95% CI)"
    ),
    y = NULL
  ) +

  theme_real_world

save_figure(
  plot = coef_plot_simple,
  stem = "Figure_real_world_method_comparison",
  directory = figure_dir,
  width = 9.6,
  height = 5.8
)


# ------------------------------------------------------------------------------
# 5B. Annotated side-by-side forest figure
#
# This version adds the adjusted coefficient (95% CI) and P value within each
# outcome panel and is used as Figure 6 in the manuscript.
#
# Output files:
#   Figure_real_world_method_comparison.*                  = unannotated figure
#   Figure_real_world_method_comparison_with_annotations.* = annotated Figure 6
# ------------------------------------------------------------------------------

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop(
    "Package 'patchwork' is required for the annotated manuscript figure.\n",
    "Install it with: install.packages('patchwork')"
  )
}

format_plot_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(
      p < 0.001,
      "<0.001",
      sprintf("%.3f", p)
    )
  )
}

method_y_map <- setNames(
  rev(
    seq_along(
      method_order
    )
  ),
  method_order
)

plot_data_annotated <- results %>%
  mutate(
    method = factor(
      as.character(method),
      levels = method_order
    ),
    outcome = factor(
      as.character(outcome),
      levels = outcome_order
    ),
    method_y = unname(
      method_y_map[
        as.character(method)
      ]
    ),
    estimate_ci_display = sprintf(
      "%.3f (%.3f to %.3f)",
      estimate,
      conf_low,
      conf_high
    ),
    p_display = format_plot_p(
      p_value
    )
  )

build_annotated_outcome_panel <- function(
    outcome_name,
    show_method_labels = TRUE
) {

  dd <- plot_data_annotated %>%
    filter(
      as.character(outcome) ==
        outcome_name
    )

  dd_base <- dd %>%
    filter(
      method != "AIPW"
    )

  dd_aipw <- dd %>%
    filter(
      method == "AIPW"
    )

  # Use a COMMON graphical x-axis through 0.20 for both outcomes.
  max_ci <- max(
    dd$conf_high,
    na.rm = TRUE
  )

  graph_plot_end <- max(
    0.235,
    max_ci + 0.003
  )

  graph_breaks <- seq(
    0,
    0.20,
    by = 0.05
  )

  estimate_x <- 0.265
  p_x <- 0.375
  panel_max <- 0.430

  ggplot(
    dd,
    aes(
      x = estimate,
      y = method_y,
      color = method
    )
  ) +

    geom_blank(
      data = tibble::tibble(
        x = panel_max,
        y = 1
      ),
      aes(
        x = x,
        y = y
      ),
      inherit.aes = FALSE
    ) +

    geom_vline(
      xintercept = 0,
      linetype = "dotted",
      linewidth = 0.65,
      color = "grey20"
    ) +

    geom_segment(
      aes(
        x = conf_low,
        xend = conf_high,
        y = method_y,
        yend = method_y
      ),
      linewidth = 0.85
    ) +

    geom_point(
      data = dd_base,
      aes(
        shape = method,
        fill = method
      ),
      size = 3.4,
      stroke = 1.15
    ) +

    geom_point(
      data = dd_aipw,
      shape = 21,
      color = "white",
      fill = "white",
      size = 4.5,
      stroke = 0,
      show.legend = FALSE
    ) +

    geom_point(
      data = dd_aipw,
      aes(
        shape = method
      ),
      size = 3.6,
      stroke = 1.20,
      show.legend = TRUE
    ) +

    geom_text(
      data = dd,
      aes(
        x = estimate_x,
        y = method_y,
        label = estimate_ci_display
      ),
      inherit.aes = FALSE,
      hjust = 0,
      color = "black",
      family = figure_font,
      size = 3.00
    ) +

    geom_text(
      data = dd,
      aes(
        x = p_x,
        y = method_y,
        label = p_display
      ),
      inherit.aes = FALSE,
      hjust = 0,
      color = "black",
      family = figure_font,
      size = 3.00
    ) +

    annotate(
      "text",
      x = estimate_x,
      y = 6.28,
      label = "Adjusted beta (95% CI)",
      hjust = 0,
      color = "black",
      family = figure_font,
      fontface = "bold",
      size = 3.05
    ) +

    annotate(
      "text",
      x = p_x,
      y = 6.28,
      label = "p-value",
      hjust = 0,
      color = "black",
      family = figure_font,
      fontface = "bold",
      size = 3.05
    ) +

    facet_wrap(
      ~ outcome,
      nrow = 1
    ) +

    scale_y_continuous(
      limits = c(
        0.45,
        6.55
      ),
      breaks = rev(
        seq_along(
          method_order
        )
      ),
      labels = if (
        show_method_labels
      ) {
        method_order
      } else {
        rep(
          "",
          length(method_order)
        )
      },
      expand = expansion(
        mult = c(
          0.01,
          0.02
        )
      )
    ) +

    scale_x_continuous(
      limits = c(
        -0.01,
        panel_max
      ),
      breaks = graph_breaks,
      labels = function(x) {
        sprintf(
          "%.2f",
          x
        )
      },
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +

    scale_color_manual(
      values = method_colors,
      breaks = method_order,
      limits = method_order,
      drop = FALSE
    ) +

    scale_fill_manual(
      values = setNames(
        rep(
          "white",
          length(method_order)
        ),
        method_order
      ),
      breaks = method_order,
      limits = method_order,
      drop = FALSE
    ) +

    scale_shape_manual(
      values = method_shapes,
      breaks = method_order,
      limits = method_order,
      drop = FALSE
    ) +

    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        label.position = "right",
        override.aes = list(
          shape = unname(
            method_shapes
          ),
          fill = rep(
            "white",
            length(method_order)
          ),
          linewidth = 1.0,
          size = rep(
            3.8,
            length(method_order)
          )
        )
      ),
      fill = "none",
      shape = "none"
    ) +

    labs(
      x = NULL,
      y = NULL
    ) +

    theme_real_world +

    theme(
      axis.text.x = element_text(
        size = 9.6
      ),

      axis.text.y = element_text(
        size = 10.5,
        margin = margin(
          r = 3
        )
      ),

      axis.ticks.y = element_line(
        color = "black"
      ),

      plot.margin = margin(
        2,
        2,
        1,
        if (
          show_method_labels
        ) {
          3
        } else {
          1
        }
      )
    )
}

grim_panel_annotated <- build_annotated_outcome_panel(
  outcome_name = "GrimAge2",
  show_method_labels = TRUE
)

pace_panel_annotated <- build_annotated_outcome_panel(
  outcome_name = "DunedinPACE",
  show_method_labels = FALSE
)

x_title_plot_annotated <- ggplot() +

  annotate(
    "text",
    x = 0.5,
    y = 0.5,
    label = paste0(
      "Standardized association with daily discrimination ",
      "(estimate and 95% CI)"
    ),
    family = figure_font,
    size = 4.10,
    color = "black"
  ) +

  coord_cartesian(
    xlim = c(
      0,
      1
    ),
    ylim = c(
      0,
      1
    )
  ) +

  theme_void(
    base_family = figure_font
  ) +

  theme(
    plot.margin = margin(
      0,
      0,
      0,
      0
    )
  )

annotated_panel_row <- patchwork::wrap_plots(
  grim_panel_annotated,
  pace_panel_annotated,
  nrow = 1,
  widths = c(
    1,
    1
  ),
  guides = "collect"
) &
  theme(
    legend.position = "bottom"
  )

coef_plot_annotated <- patchwork::wrap_plots(
  annotated_panel_row,
  x_title_plot_annotated,
  ncol = 1,
  heights = c(
    1,
    0.045
  ),
  guides = "collect"
) &
  theme(
    legend.position = "bottom"
  )

save_figure(
  plot = coef_plot_annotated,
  stem = "Figure_real_world_method_comparison_with_annotations",
  directory = figure_dir,
  width = 12.8,
  height = 5.9
)

# ------------------------------------------------------------------------------
# 6. Method differences relative to AIPW
# ------------------------------------------------------------------------------

aipw_reference <- results %>%
  filter(method == "AIPW") %>%
  select(
    outcome,
    aipw_estimate = estimate,
    aipw_se = se
  )

method_differences <- results %>%
  left_join(
    aipw_reference,
    by = "outcome"
  ) %>%
  mutate(
    estimate_difference_from_AIPW =
      estimate - aipw_estimate,

    absolute_difference_from_AIPW =
      abs(estimate - aipw_estimate),

    percent_difference_from_AIPW =
      100 * (estimate - aipw_estimate) /
      abs(aipw_estimate),

    se_difference_from_AIPW =
      se - aipw_se
  ) %>%
  select(
    outcome,
    method,
    estimate,
    aipw_estimate,
    estimate_difference_from_AIPW,
    absolute_difference_from_AIPW,
    percent_difference_from_AIPW,
    se,
    aipw_se,
    se_difference_from_AIPW
  ) %>%
  arrange(
    outcome,
    method
  )

write_csv(
  method_differences,
  file.path(
    table_dir,
    "real_world_method_differences_from_AIPW.csv"
  )
)

# ------------------------------------------------------------------------------
# 7. Weight diagnostics
# ------------------------------------------------------------------------------

weight_summary <- weights %>%
  mutate(
    ess_fraction_of_phase2 =
      weight_ess / sample_summary$N_phase2[[1]]
  ) %>%
  select(
    outcome,
    method,
    weight_min,
    weight_p99,
    weight_max,
    weight_cv,
    weight_ess,
    ess_fraction_of_phase2
  )

write_csv(
  weight_summary,
  file.path(
    diagnostic_dir,
    "real_world_weight_diagnostics_summary.csv"
  )
)

# Diagnostic figure: maximum and 99th-percentile weights.

weight_plot_data <- weights %>%
  select(
    outcome,
    method,
    weight_p99,
    weight_max
  ) %>%
  pivot_longer(
    cols = c(weight_p99, weight_max),
    names_to = "weight_statistic",
    values_to = "weight"
  ) %>%
  mutate(
    statistic_index = case_when(
      weight_statistic == "weight_p99" ~ 1,
      weight_statistic == "weight_max" ~ 2,
      TRUE ~ NA_real_
    ),
    method_offset = case_when(
      as.character(method) == "IPW" ~ -0.08,
      as.character(method) == "AIPW" ~ 0.08,
      TRUE ~ 0
    ),
    x_plot = statistic_index + method_offset
  )

weight_base <- weight_plot_data %>%
  filter(method != "AIPW")

weight_aipw <- weight_plot_data %>%
  filter(method == "AIPW")

weight_plot <- ggplot(
  weight_plot_data,
  aes(
    x = x_plot,
    y = weight,
    color = method
  )
) +

  geom_point(
    data = weight_base,
    aes(
      shape = method,
      fill = method
    ),
    size = 3.4,
    stroke = 1.15
  ) +

  # White mask behind native AIPW shape 10.
  geom_point(
    data = weight_aipw,
    shape = 21,
    color = "white",
    fill = "white",
    size = 4.5,
    stroke = 0,
    show.legend = FALSE
  ) +

  geom_point(
    data = weight_aipw,
    aes(
      shape = method
    ),
    size = 3.6,
    stroke = 1.20,
    show.legend = TRUE
  ) +

  facet_wrap(
    ~ outcome,
    nrow = 1
  ) +

  scale_x_continuous(
    breaks = c(1, 2),
    labels = c(
      "99th percentile",
      "Maximum"
    ),
    limits = c(0.65, 2.35)
  ) +

  scale_color_manual(
    values = method_colors,
    breaks = c("IPW", "AIPW")
  ) +

  scale_fill_manual(
    values = c(
      "IPW" = "white",
      "AIPW" = "white"
    ),
    breaks = c("IPW", "AIPW")
  ) +

  scale_shape_manual(
    values = method_shapes,
    breaks = c("IPW", "AIPW")
  ) +

  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      label.position = "right",
      override.aes = list(
        shape = unname(
          method_shapes[c("IPW", "AIPW")]
        ),
        fill = c("white", "white"),
        size = c(3.8, 3.8)
      )
    ),
    fill = "none",
    shape = "none"
  ) +

  labs(
    x = NULL,
    y = "Inverse-probability weight"
  ) +

  theme_real_world

save_figure(
  plot = weight_plot,
  stem = "Diagnostic_real_world_weights",
  directory = diagnostic_dir,
  width = 8.0,
  height = 4.8
)

# ------------------------------------------------------------------------------
# 8. Compact key-findings file
# ------------------------------------------------------------------------------

key_findings <- results %>%
  filter(
    method %in% c(
      "FCS-MI",
      "JM-MI",
      "AIPW"
    )
  ) %>%
  select(
    outcome,
    method,
    estimate,
    se,
    conf_low,
    conf_high,
    p_value
  )

write_csv(
  key_findings,
  file.path(
    table_dir,
    "real_world_FCS_JM_AIPW_comparison.csv"
  )
)

# ------------------------------------------------------------------------------
# 9. Output manifest / README
# ------------------------------------------------------------------------------

readme_lines <- c(
  "MIDUS REAL-WORLD APPLICATION MANUSCRIPT OUTPUTS",
  "",
  "This directory contains manuscript-ready outputs created from the completed",
  "real-world analysis. Script 45 does not rerun MI, IPW, or AIPW.",
  "",
  "tables/Table_real_world_application.csv",
  "  Rounded manuscript-ready comparison of the six methods.",
  "",
  "tables/Table_real_world_application_unrounded.csv",
  "  Machine-readable unrounded results.",
  "",
  "tables/real_world_method_differences_from_AIPW.csv",
  "  Method-specific estimate and SE differences relative to AIPW.",
  "",
  "tables/real_world_FCS_JM_AIPW_comparison.csv",
  "  Focused comparison of FCS-MI, JM-MI, and AIPW.",
  "",
  "figures/Figure_real_world_method_comparison.png/.pdf",
  "  Primary coefficient/95% CI figure.",
  "",
  "diagnostics/real_world_weight_diagnostics_summary.csv",
  "  IPW/AIPW weight diagnostics including ESS fraction.",
  "",
  "diagnostics/Diagnostic_real_world_weights.png/.pdf",
  "  Weight-diagnostic figure."
)

writeLines(
  readme_lines,
  file.path(
    manuscript_dir,
    "README_REAL_WORLD_OUTPUTS.txt"
  )
)

cat(
  "\nREAL-WORLD MANUSCRIPT OUTPUTS COMPLETE\n",
  "Output directory: ", manuscript_dir, "\n",
  "Primary table: ",
  file.path(
    table_dir,
    "Table_real_world_application.csv"
  ),
  "\nPrimary figure: ",
  file.path(
    figure_dir,
    "Figure_real_world_method_comparison.png"
  ),
  "\n",
  sep = ""
)
