################################################################################
# 09_make_simulation_figures.R
#
# MIDUS two-phase simulation study
#
# --------------------------------
# This script:
#   1. Uses the validated primary simulation summaries already produced.
#   2. Generates a complete full-factorial figure library.
#   3. Generates weak / MIDUS-like / strong manuscript-candidate figures.
#   4. Uses style visual grammar:
#        - one performance measure per figure
#        - Phase-1 sample size on the x-axis
#        - Phase-2 sampling condition in columns
#        - true effect size in rows when applicable
#        - consistent method color + open symbol + line type
#        - one-row method legend
#        - target and acceptable-boundary reference lines
#   5. Uses fixed y-axis scales within each performance measure so figures are
#      directly comparable across R_A^2 x R_Y^2 settings.
#   6. Uses plotmath labels plus cairo_pdf so theta and phi_2 render correctly
#      in PDF files.
#
# IMPORTANT
# ---------
# Run this script from the ROOT of the transferred project:
#
#   midus_two_phase_sim/
#
# Example:
#   source("09_make_simulation_figures.R")
#
# Nothing in this script reruns the simulations or changes the validated results.
################################################################################


# ==============================================================================
# 0. PACKAGES
# ==============================================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "ggplot2",
  "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install these packages first: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(scales)


# ==============================================================================
# 1. INPUT FILES
# ==============================================================================

performance_file <- "results/summary/primary_performance_wide.csv"
failure_file     <- "results/summary/primary_failure_rates.csv"
runtime_file     <- "results/tables/table_runtime_overview.csv"

required_files <- c(
  performance_file,
  failure_file,
  runtime_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0L) {
  stop(
    "Required file(s) not found:\n",
    paste0("  - ", missing_files, collapse = "\n"),
    "\n\nRun this script from the root of midus_two_phase_sim."
  )
}


# ==============================================================================
# 2. OUTPUT DIRECTORIES
# ==============================================================================

figure_root <- "results/figures"

full_root <- file.path(
  figure_root,
  "01_full_factorial"
)

candidate_root <- file.path(
  figure_root,
  "02_manuscript_candidates"
)

final_root <- file.path(
  figure_root,
  "03_final_manuscript"
)

diagnostic_root <- file.path(
  figure_root,
  "04_diagnostics"
)

data_root <- "results/figure_data"

for (d in c(
  figure_root,
  full_root,
  candidate_root,
  final_root,
  diagnostic_root,
  data_root
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ==============================================================================
# 3. READ DATA
# ==============================================================================

perf <- read_csv(
  performance_file,
  show_col_types = FALSE
)

failures <- read_csv(
  failure_file,
  show_col_types = FALSE
)

runtime <- read_csv(
  runtime_file,
  show_col_types = FALSE
)


# ==============================================================================
# 4. EXPECTED METHODS AND ORDER
# ==============================================================================

method_order <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

perf <- perf %>%
  mutate(
    method = factor(
      method,
      levels = method_order
    )
  )

failures <- failures %>%
  mutate(
    method = factor(
      method,
      levels = method_order
    )
  )

runtime <- runtime %>%
  mutate(
    method = factor(
      method,
      levels = method_order
    )
  )


# ==============================================================================
# 5. METHOD COLORS, OPEN SYMBOLS, AND LINE TYPES
#
# Colors:
#   RColorBrewer Dark2 sequence
#
# Symbols:
#   The first five methods use white-filled open plotting symbols so the
#   connecting line does not pass visibly through the symbol.
#
#   Naive  = open circle
#   CCA    = open triangle up
#   FCS-MI = open square
#   JM-MI  = open diamond
#   IPW    = open triangle down
#
#   AIPW uses native plotting symbol 10 (open circle-plus), drawn over a white
#   masking circle so the connecting line is visually blocked in the center.
#
# Line types:
#   Use standard base ggplot line types only.
# ==============================================================================

method_colors <- c(
  "Naive"  = "#1B9E77",
  "CCA"    = "#D95F02",
  "FCS-MI" = "#7570B3",
  "JM-MI"  = "#E7298A",
  "IPW"    = "#66A61E",
  "AIPW"   = "#E6AB02"
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
# 6. DESIGN LEVELS
# ==============================================================================

N_levels <- sort(
  unique(perf$N)
)

phase2_levels <- sort(
  unique(perf$phase2_fraction)
)

theta_levels <- sort(
  unique(perf$theta_true)
)

rA_levels <- sort(
  unique(perf$r2_a_marker)
)

rY_levels <- sort(
  unique(perf$r2_y_marker)
)


# ==============================================================================
# 7. PLOT-SAFE DESIGN LABELS
#
# N is treated as a DESIGN FACTOR rather than a continuous numerical x-axis.
# This makes 500, 800, and 1,500 equally spaced, as in a factorial simulation
# figure.
#
# phi_2 is formally used as the target Phase-2 sampling fraction in the figure
# strips:
#
#   Low Phase-2 sampling
#   phi_2 = 0.25
#
#   Moderate Phase-2 sampling
#   phi_2 = 0.50
#
#   High Phase-2 sampling
#   phi_2 = 0.65
#
# The strip expressions are parsed by plotmath so phi[2] renders in PNG and PDF.
# ==============================================================================

perf <- perf %>%
  mutate(
    N_plot = factor(
      N,
      levels = c(500, 800, 1500),
      labels = c("500", "800", "1,500")
    ),
    
    phase2_plot = factor(
      phase2_fraction,
      levels = c(0.25, 0.50, 0.65)
    ),
    
    theta_plot = factor(
      theta_true,
      levels = theta_levels
    )
  )


phase2_labeller <- as_labeller(
  c(
    "0.25" = 'atop("Low Phase-2 sampling", phi[2] == "0.25")',
    "0.5"  = 'atop("Moderate Phase-2 sampling", phi[2] == "0.50")',
    "0.65" = 'atop("High Phase-2 sampling", phi[2] == "0.65")'
  ),
  label_parsed
)


theta_labeller <- as_labeller(
  setNames(
    paste0(
      'theta == "',
      sprintf("%.2f", theta_levels),
      '"'
    ),
    as.character(theta_levels)
  ),
  label_parsed
)


# ==============================================================================
# 8. FIGURE THEME
#
# Important:
#   - white background
#   - gray strip headers
#   - no panel grid
#   - one-row method legend
#   - strip font sized so "Moderate Phase-2 sampling" remains on ONE line
# ==============================================================================

figure_font <- if (.Platform$OS.type == "windows") "Arial" else "sans"

theme_simulation <- theme_bw(
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


# ==============================================================================
# 9. HELPERS
# ==============================================================================

num_code <- function(
    x,
    digits = 2
) {
  
  gsub(
    "\\.",
    "",
    formatC(
      x,
      format = "f",
      digits = digits
    )
  )
}


save_both <- function(
    plot,
    filename,
    directory,
    width = 9.6,
    height = 5.8
) {
  
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  ggsave(
    filename = file.path(
      directory,
      paste0(filename, ".png")
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
      paste0(filename, ".pdf")
    ),
    plot = plot,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    bg = "white"
  )
}


add_method_scales <- function(
    p
) {
  
  p +
    scale_color_manual(
      values = method_colors,
      breaks = method_order,
      limits = method_order,
      drop = FALSE
    ) +
    
    scale_fill_manual(
      values = setNames(
        rep("white", length(method_order)),
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
    
    scale_linetype_manual(
      values = method_linetypes,
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
          shape = unname(method_shapes),
          fill = rep("white", length(method_order)),
          linewidth = 1.0,
          size = rep(3.8, length(method_order))
        )
      ),
      fill = "none",
      shape = "none",
      linetype = "none"
    )
}

add_method_geoms <- function(
    p,
    data
) {
  
  base_data <- data %>%
    filter(method != "AIPW")
  
  aipw_data <- data %>%
    filter(method == "AIPW")
  
  p <- p +
    geom_line(
      linewidth = 1.00
    ) +
    
    # White-filled symbols 21-25 mask the connecting line in their interiors.
    geom_point(
      data = base_data,
      aes(
        shape = method,
        fill = method
      ),
      size = 3.4,
      stroke = 1.15
    )
  
  if (nrow(aipw_data) > 0L) {
    
    # AIPW uses native plotting symbol 10 (open circle-plus).
    # First draw a white mask so the connecting line does not show through
    # the center, then draw shape 10 on top. The top layer participates in
    # the legend so AIPW has a visible point symbol in the one-row key.
    p <- p +
      geom_point(
        data = aipw_data,
        shape = 21,
        color = "white",
        fill = "white",
        size = 4.5,
        stroke = 0,
        show.legend = FALSE
      ) +
      
      geom_point(
        data = aipw_data,
        aes(
          shape = method
        ),
        size = 3.6,
        stroke = 1.20,
        show.legend = TRUE
      )
  }
  
  p
}


# ==============================================================================
# 10. FIXED AXIS SCALES
#
# Same performance measure = same scale across every R_A^2 x R_Y^2 figure.
#
# These limits cover the observed primary production results.
# ==============================================================================

apply_metric_axis <- function(
    p,
    metric
) {
  
  if (metric == "bias") {
    
    p <- p +
      scale_y_continuous(
        limits = c(
          -0.01,
          0.025
        ),
        breaks = c(
          -0.01,
          0,
          0.01,
          0.02
        ),
        labels = label_number(
          accuracy = 0.01
        )
      )
  }
  
  
  if (metric == "se_ratio") {
    
    p <- p +
      scale_y_continuous(
        limits = c(
          0.85,
          1.15
        ),
        breaks = c(
          0.90,
          1.00,
          1.10
        ),
        labels = label_number(
          accuracy = 0.01
        )
      )
  }
  
  
  if (metric %in% c(
    "coverage",
    "bias_eliminated_coverage"
  )) {
    
    p <- p +
      scale_y_continuous(
        limits = c(
          0.84,
          0.98
        ),
        breaks = c(
          0.85,
          0.90,
          0.93,
          0.95,
          0.97
        ),
        labels = percent_format(
          accuracy = 1
        )
      )
  }
  
  
  if (metric == "empirical_se") {
    
    p <- p +
      scale_y_continuous(
        limits = c(
          0.02,
          0.105
        ),
        breaks = c(
          0.02,
          0.04,
          0.06,
          0.08,
          0.10
        )
      )
  }
  
  
  if (metric == "model_se") {
    
    p <- p +
      scale_y_continuous(
        limits = c(
          0.02,
          0.105
        ),
        breaks = c(
          0.02,
          0.04,
          0.06,
          0.08,
          0.10
        )
      )
  }
  
  
  if (metric == "mse") {
    
    p <- p +
      scale_y_continuous(
        limits = c(
          0,
          0.011
        ),
        breaks = c(
          0,
          0.0025,
          0.0050,
          0.0075,
          0.0100
        )
      )
  }
  
  
  if (metric == "relative_error") {
    
    p <- p +
      scale_y_continuous(
        limits = c(
          -15,
          15
        ),
        breaks = c(
          -15,
          -10,
          0,
          10,
          15
        ),
        labels = function(x) {
          paste0(
            x,
            "%"
          )
        }
      )
  }
  
  
  # Relative precision has a much wider observed range, including negative
  # values and values >400%. It remains available but is intentionally not
  # forced into an artificially narrow fixed scale.
  if (metric == "relative_precision") {
    
    p <- p +
      scale_y_continuous(
        limits = c(
          -25,
          450
        ),
        breaks = c(
          0,
          100,
          200,
          300,
          400
        ),
        labels = function(x) {
          paste0(
            x,
            "%"
          )
        }
      )
  }
  
  
  p
}


# ==============================================================================
# 11. REFERENCE LINES
#
# No shaded bands.
#
# Bias:
#   dotted 0
#
# SE ratio:
#   dotted 1.00
#   solid gray 0.90 and 1.10
#
# Coverage:
#   dotted 0.95
#   solid gray 0.93 and 0.97
#
# Relative SE error:
#   dotted 0
#   solid gray -10 and +10
# ==============================================================================

add_reference_lines <- function(
    p,
    metric
) {
  
  if (metric == "bias") {
    
    p <- p +
      geom_hline(
        yintercept = 0,
        linetype = "dotted",
        linewidth = 0.65,
        color = "grey20"
      )
  }
  
  
  if (metric == "se_ratio") {
    
    p <- p +
      geom_hline(
        yintercept = c(
          0.90,
          1.10
        ),
        linetype = "solid",
        linewidth = 0.55,
        color = "grey55"
      ) +
      
      geom_hline(
        yintercept = 1.00,
        linetype = "dotted",
        linewidth = 0.65,
        color = "grey20"
      )
  }
  
  
  if (metric %in% c(
    "coverage",
    "bias_eliminated_coverage"
  )) {
    
    p <- p +
      geom_hline(
        yintercept = c(
          0.93,
          0.97
        ),
        linetype = "solid",
        linewidth = 0.55,
        color = "grey55"
      ) +
      
      geom_hline(
        yintercept = 0.95,
        linetype = "dotted",
        linewidth = 0.65,
        color = "grey20"
      )
  }
  
  
  if (metric == "relative_error") {
    
    p <- p +
      geom_hline(
        yintercept = c(
          -10,
          10
        ),
        linetype = "solid",
        linewidth = 0.55,
        color = "grey55"
      ) +
      
      geom_hline(
        yintercept = 0,
        linetype = "dotted",
        linewidth = 0.65,
        color = "grey20"
      )
  }
  
  
  if (metric == "relative_precision") {
    
    p <- p +
      geom_hline(
        yintercept = 0,
        linetype = "dotted",
        linewidth = 0.65,
        color = "grey20"
      )
  }
  
  
  p
}


# ==============================================================================
# 12. PERFORMANCE-MEASURE SPECIFICATIONS
#
# These are the prespecified performance measures already generated by the
# existing rsimsum/post-processing pipeline.
# ==============================================================================

metric_specs <- list(
  
  mean_estimate = list(
    variable = "est__thetamean",
    y_label = expression(
      paste(
        "Mean estimate of ",
        theta
      )
    )
  ),
  
  bias = list(
    variable = "est__bias",
    y_label = "Empirical bias"
  ),
  
  empirical_se = list(
    variable = "est__empse",
    y_label = "Empirical SE (EmpSE)"
  ),
  
  model_se = list(
    variable = "est__modelse",
    y_label = "Average model-based SE (ModSE)"
  ),
  
  se_ratio = list(
    variable = "se_ratio",
    y_label = "SE ratio (ModSE / EmpSE)"
  ),
  
  mse = list(
    variable = "est__mse",
    y_label = "Mean squared error"
  ),
  
  relative_precision = list(
    variable = "est__relprec",
    y_label = "Relative precision"
  ),
  
  relative_error = list(
    variable = "est__relerror",
    y_label = "Relative error in model-based SE"
  ),
  
  coverage = list(
    variable = "est__cover",
    y_label = "Coverage rate"
  ),
  
  bias_eliminated_coverage = list(
    variable = "est__becover",
    y_label = "Bias-eliminated coverage rate"
  )
)


# ==============================================================================
# 13. GENERIC FACTORIAL FIGURE
#
# One figure = one R_A^2 x R_Y^2 combination.
#
# x-axis  = N (equally spaced design levels)
# columns = Phase-2 sampling condition
# rows    = theta
# lines   = methods
# ==============================================================================

make_factorial_metric_plot <- function(
    data,
    rA,
    rY,
    metric_name,
    theta_values = theta_levels,
    output_directory,
    file_prefix = NULL,
    save_plot_data = TRUE
) {
  
  spec <- metric_specs[[metric_name]]
  
  if (is.null(spec)) {
    stop(
      "Unknown metric: ",
      metric_name
    )
  }
  
  
  dat <- data %>%
    filter(
      near(
        r2_a_marker,
        rA
      ),
      near(
        r2_y_marker,
        rY
      ),
      theta_true %in% theta_values
    ) %>%
    mutate(
      plot_value = .data[[spec$variable]],
      
      theta_plot = factor(
        theta_true,
        levels = theta_values
      )
    )
  
  
  if (nrow(dat) == 0L) {
    return(
      invisible(NULL)
    )
  }
  
  
  # Build thphi labels dynamically for the subset used in this plot.
  theta_labeller_use <- as_labeller(
    setNames(
      paste0(
        'theta == "',
        sprintf("%.2f", theta_values),
        '"'
      ),
      as.character(theta_values)
    ),
    label_parsed
  )
  
  
  p <- ggplot(
    dat,
    aes(
      x = N_plot,
      y = plot_value,
      color = method,
      linetype = method,
      group = method
    )
  )
  
  
  # Mean estimate needs a row-specific true theta reference.
  if (metric_name == "mean_estimate") {
    
    theta_refs <- dat %>%
      distinct(
        theta_plot,
        theta_true
      )
    
    p <- p +
      geom_hline(
        data = theta_refs,
        aes(
          yintercept = theta_true
        ),
        linetype = "dotted",
        linewidth = 0.65,
        color = "grey20"
      )
    
  } else {
    
    p <- add_reference_lines(
      p,
      metric_name
    )
  }
  
  
  p <- add_method_geoms(
    p,
    dat
  ) +
    
    facet_grid(
      rows = vars(
        theta_plot
      ),
      cols = vars(
        phase2_plot
      ),
      labeller = labeller(
        theta_plot = theta_labeller_use,
        phase2_plot = phase2_labeller
      )
    ) +
    
    labs(
      x = "Phase-1 sample size (N)",
      y = spec$y_label
    ) +
    
    theme_simulation
  
  
  p <- add_method_scales(
    p
  )
  
  
  # Fixed y-axis by performance measure.
  if (metric_name != "mean_estimate") {
    
    p <- apply_metric_axis(
      p,
      metric_name
    )
  }
  
  
  # For mean estimate, allow the rows to use the common global theta range.
  if (metric_name == "mean_estimate") {
    
    p <- p +
      scale_y_continuous(
        limits = c(
          -0.025,
          0.33
        )
      )
  }
  
  
  if (is.null(file_prefix)) {
    
    file_prefix <- paste0(
      metric_name,
      "_RA",
      num_code(
        rA,
        2
      ),
      "_RY",
      num_code(
        rY,
        2
      )
    )
  }
  
  
  save_both(
    p,
    filename = file_prefix,
    directory = output_directory,
    width = 9.6,
    height = ifelse(
      length(theta_values) == 3L,
      7.4,
      5.9
    )
  )
  
  
  if (isTRUE(save_plot_data)) {
    
    metric_data_dir <- file.path(
      data_root,
      metric_name
    )
    
    dir.create(
      metric_data_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    write_csv(
      dat %>%
        select(
          scenario_id,
          N,
          phase2_fraction,
          r2_a_marker,
          r2_y_marker,
          theta_true,
          method,
          all_of(
            spec$variable
          )
        ),
      file.path(
        metric_data_dir,
        paste0(
          file_prefix,
          ".csv"
        )
      )
    )
  }
  
  
  invisible(
    p
  )
}


# ==============================================================================
# 14. NULL REJECTION RATE
#
# theta = 0 only.
#
# y-axis:
#   0%-15%
#
# reference lines:
#   solid gray = 4%, 6%
#   dotted     = 5%
# ==============================================================================

make_type1_plot <- function(
    data,
    rA,
    rY,
    output_directory,
    file_prefix = NULL,
    save_plot_data = TRUE
) {
  
  dat <- data %>%
    filter(
      near(
        r2_a_marker,
        rA
      ),
      near(
        r2_y_marker,
        rY
      ),
      near(
        theta_true,
        0
      )
    ) %>%
    mutate(
      type1_error = est__power
    )
  
  
  p <- ggplot(
    dat,
    aes(
      x = N_plot,
      y = type1_error,
      color = method,
      linetype = method,
      group = method
    )
  ) +
    
    geom_hline(
      yintercept = c(
        0.04,
        0.06
      ),
      linetype = "solid",
      linewidth = 0.55,
      color = "grey55"
    ) +
    
    geom_hline(
      yintercept = 0.05,
      linetype = "dotted",
      linewidth = 0.65,
      color = "grey20"
    ) +
    
    facet_grid(
      . ~ phase2_plot,
      labeller = labeller(
        phase2_plot = phase2_labeller
      )
    ) +
    
    scale_y_continuous(
      limits = c(
        0,
        0.15
      ),
      breaks = c(
        0,
        0.04,
        0.05,
        0.06,
        0.10,
        0.15
      ),
      labels = percent_format(
        accuracy = 1
      )
    ) +
    
    labs(
      x = "Phase-1 sample size (N)",
      y = "Null rejection rate"
    ) +
    
    theme_simulation
  
  
  p <- add_method_geoms(
    p,
    dat
  )
  
  p <- add_method_scales(
    p
  )
  
  
  if (is.null(file_prefix)) {
    
    file_prefix <- paste0(
      "type1_RA",
      num_code(
        rA,
        2
      ),
      "_RY",
      num_code(
        rY,
        2
      )
    )
  }
  
  
  save_both(
    p,
    filename = file_prefix,
    directory = output_directory,
    width = 9.6,
    height = 4.6
  )
  
  
  if (isTRUE(save_plot_data)) {
    
    metric_data_dir <- file.path(
      data_root,
      "type1_error"
    )
    
    dir.create(
      metric_data_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    write_csv(
      dat %>%
        select(
          scenario_id,
          N,
          phase2_fraction,
          r2_a_marker,
          r2_y_marker,
          theta_true,
          method,
          type1_error,
          mcse__power
        ),
      file.path(
        metric_data_dir,
        paste0(
          file_prefix,
          ".csv"
        )
      )
    )
  }
  
  
  invisible(
    p
  )
}


# ==============================================================================
# 15. POWER
#
# theta > 0 only.
#
# y-axis:
#   20%-100%
#
# dotted reference:
#   80%
# ==============================================================================

make_power_plot <- function(
    data,
    rA,
    rY,
    output_directory,
    file_prefix = NULL,
    save_plot_data = TRUE
) {
  
  alternative_theta <- theta_levels[
    theta_levels > 0
  ]
  
  
  dat <- data %>%
    filter(
      near(
        r2_a_marker,
        rA
      ),
      near(
        r2_y_marker,
        rY
      ),
      theta_true > 0
    ) %>%
    mutate(
      theta_plot = factor(
        theta_true,
        levels = alternative_theta
      )
    )
  
  
  theta_labeller_alt <- as_labeller(
    setNames(
      paste0(
        'theta == "',
        sprintf("%.2f", alternative_theta),
        '"'
      ),
      as.character(
        alternative_theta
      )
    ),
    label_parsed
  )
  
  
  p <- ggplot(
    dat,
    aes(
      x = N_plot,
      y = est__power,
      color = method,
      linetype = method,
      group = method
    )
  ) +
    
    geom_hline(
      yintercept = 0.80,
      linetype = "dotted",
      linewidth = 0.65,
      color = "grey30"
    ) +
    
    facet_grid(
      rows = vars(
        theta_plot
      ),
      cols = vars(
        phase2_plot
      ),
      labeller = labeller(
        theta_plot = theta_labeller_alt,
        phase2_plot = phase2_labeller
      )
    ) +
    
    scale_y_continuous(
      limits = c(
        0.20,
        1.00
      ),
      breaks = c(
        0.20,
        0.40,
        0.60,
        0.80,
        1.00
      ),
      labels = percent_format(
        accuracy = 1
      )
    ) +
    
    labs(
      x = "Phase-1 sample size (N)",
      y = "Power"
    ) +
    
    theme_simulation
  
  
  p <- add_method_geoms(
    p,
    dat
  )
  
  p <- add_method_scales(
    p
  )
  
  
  if (is.null(file_prefix)) {
    
    file_prefix <- paste0(
      "power_RA",
      num_code(
        rA,
        2
      ),
      "_RY",
      num_code(
        rY,
        2
      )
    )
  }
  
  
  save_both(
    p,
    filename = file_prefix,
    directory = output_directory,
    width = 9.6,
    height = 5.9
  )
  
  
  if (isTRUE(save_plot_data)) {
    
    metric_data_dir <- file.path(
      data_root,
      "power"
    )
    
    dir.create(
      metric_data_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    write_csv(
      dat %>%
        select(
          scenario_id,
          N,
          phase2_fraction,
          r2_a_marker,
          r2_y_marker,
          theta_true,
          method,
          est__power,
          mcse__power
        ),
      file.path(
        metric_data_dir,
        paste0(
          file_prefix,
          ".csv"
        )
      )
    )
  }
  
  
  invisible(
    p
  )
}


# ==============================================================================
# 16. FULL FACTORIAL LIBRARY
#
# For each R_A^2 x R_Y^2 combination, generate all prespecified metric figures.
# ==============================================================================

cat(
  "\nCreating full-factorial figure library...\n"
)


for (metric_name in names(
  metric_specs
)) {
  
  metric_dir <- file.path(
    full_root,
    metric_name
  )
  
  dir.create(
    metric_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  for (rA in rA_levels) {
    
    for (rY in rY_levels) {
      
      make_factorial_metric_plot(
        data = perf,
        rA = rA,
        rY = rY,
        metric_name = metric_name,
        theta_values = theta_levels,
        output_directory = metric_dir
      )
    }
  }
}


# Null rejection rate (theta = 0).

type1_dir <- file.path(
  full_root,
  "type1_error"
)

for (rA in rA_levels) {
  
  for (rY in rY_levels) {
    
    make_type1_plot(
      data = perf,
      rA = rA,
      rY = rY,
      output_directory = type1_dir
    )
  }
}


# Power.

power_dir <- file.path(
  full_root,
  "power"
)

for (rA in rA_levels) {
  
  for (rY in rY_levels) {
    
    make_power_plot(
      data = perf,
      rA = rA,
      rY = rY,
      output_directory = power_dir
    )
  }
}


# ==============================================================================
# 17. FAILURE RATE
#
# Production result: all six methods had zero failures.
# A figure adds no useful information, so stores the failure-rate data only.
# ==============================================================================

write_csv(
  failures,
  file.path(
    data_root,
    "failure_rates_all_scenarios.csv"
  )
)

cat("Failure-rate figure omitted: all methods had zero recorded failures.\n")


# ==============================================================================
# 18. RUNTIME — REDESIGNED
#
# Horizontal dot plot:
#
#   y-axis = method
#   x-axis = median runtime per simulated dataset
#   columns = N = 500 / 800 / 1,500
#   x-axis = log scale
#
# This is much easier to read when JM-MI is orders of magnitude slower.
# ==============================================================================

runtime_plot_data <- runtime %>%
  mutate(
    # Reverse factor levels on the discrete y-axis so the displayed
    # top-to-bottom order matches all other figures:
    # Naive, CCA, FCS-MI, JM-MI, IPW, AIPW.
    method = factor(
      as.character(method),
      levels = rev(method_order)
    ),
    
    N_plot = factor(
      N,
      levels = c(
        500,
        800,
        1500
      ),
      labels = c(
        "N = 500",
        "N = 800",
        "N = 1,500"
      )
    ),
    
    runtime_label = case_when(
      median_seconds >= 60 ~ paste0(
        round(
          median_seconds / 60,
          1
        ),
        " min"
      ),
      
      median_seconds >= 1 ~ paste0(
        round(
          median_seconds,
          1
        ),
        " s"
      ),
      
      TRUE ~ paste0(
        round(
          median_seconds,
          3
        ),
        " s"
      )
    )
  )


write_csv(
  runtime_plot_data,
  file.path(
    data_root,
    "runtime_overview_data.csv"
  )
)


runtime_base <- runtime_plot_data %>%
  filter(method != "AIPW")

runtime_aipw <- runtime_plot_data %>%
  filter(method == "AIPW")

p_runtime <- ggplot(
  runtime_plot_data,
  aes(
    x = median_seconds,
    y = method,
    color = method
  )
) +
  
  geom_point(
    data = runtime_base,
    aes(
      shape = method,
      fill = method
    ),
    size = 3.8,
    stroke = 1.15
  ) +
  
  geom_point(
    data = runtime_aipw,
    shape = 21,
    color = "white",
    fill = "white",
    size = 4.6,
    stroke = 0,
    show.legend = FALSE
  ) +
  
  geom_point(
    data = runtime_aipw,
    aes(
      shape = method
    ),
    size = 3.8,
    stroke = 1.20,
    show.legend = FALSE
  ) +
  
  geom_text(
    aes(
      label = runtime_label
    ),
    hjust = -0.20,
    size = 3.0,
    family = figure_font,
    show.legend = FALSE
  ) +
  
  facet_wrap(
    ~ N_plot,
    nrow = 1
  ) +
  
  scale_y_discrete(
    limits = rev(method_order)
  ) +
  
  scale_x_log10(
    breaks = c(
      0.001,
      0.01,
      0.1,
      1,
      10,
      100,
      1000
    ),
    labels = c(
      "0.001",
      "0.01",
      "0.1",
      "1",
      "10",
      "100",
      "1,000"
    ),
    expand = expansion(
      mult = c(
        0.08,
        0.42
      )
    )
  ) +
  
  labs(
    x = "Median runtime per simulated dataset (seconds, log scale)",
    y = NULL
  ) +
  
  theme_simulation +
  
  theme(
    legend.position = "none"
  )


p_runtime <- p_runtime +
  scale_color_manual(
    values = method_colors,
    breaks = method_order
  ) +
  scale_fill_manual(
    values = setNames(
      rep("white", length(method_order)),
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
  )


save_both(
  p_runtime,
  filename = "runtime_by_method_and_N",
  directory = diagnostic_root,
  width = 10.0,
  height = 5.0
)


# ==============================================================================
# 19. MANUSCRIPT-CANDIDATE REGIMES
#
# Generate comprehensively first; select parsimoniously second.
#
# Weak:
#   R_A^2 = 0.01
#   R_Y^2 = 0.02
#
# MIDUS-like:
#   R_A^2 = 0.05
#   R_Y^2 = 0.05
#
# Strong:
#   R_A^2 = 0.10
#   R_Y^2 = 0.10
# ==============================================================================

candidate_regimes <- tibble(
  regime = c(
    "weak_auxiliary_information",
    "MIDUS_like",
    "strong_auxiliary_information"
  ),
  
  rA = c(
    0.01,
    0.05,
    0.10
  ),
  
  rY = c(
    0.02,
    0.05,
    0.10
  )
)


write_csv(
  candidate_regimes,
  file.path(
    data_root,
    "manuscript_candidate_regimes.csv"
  )
)


alternative_theta <- theta_levels[
  theta_levels > 0
]


for (i in seq_len(
  nrow(
    candidate_regimes
  )
)) {
  
  regime_name <- candidate_regimes$regime[i]
  rA <- candidate_regimes$rA[i]
  rY <- candidate_regimes$rY[i]
  
  
  regime_dir <- file.path(
    candidate_root,
    regime_name
  )
  
  dir.create(
    regime_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # Candidate non-null figures.
  for (metric_name in names(
    metric_specs
  )) {
    
    make_factorial_metric_plot(
      data = perf,
      rA = rA,
      rY = rY,
      metric_name = metric_name,
      theta_values = alternative_theta,
      output_directory = regime_dir,
      file_prefix = paste0(
        metric_name,
        "_",
        regime_name
      ),
      save_plot_data = TRUE
    )
  }
  
  
  # Candidate null rejection-rate figure.
  make_type1_plot(
    data = perf,
    rA = rA,
    rY = rY,
    output_directory = regime_dir,
    file_prefix = paste0(
      "type1_",
      regime_name
    ),
    save_plot_data = TRUE
  )
  
  
  # Candidate power.
  make_power_plot(
    data = perf,
    rA = rA,
    rY = rY,
    output_directory = regime_dir,
    file_prefix = paste0(
      "power_",
      regime_name
    ),
    save_plot_data = TRUE
  )
}


# Copy redesigned runtime candidate.

save_both(
  p_runtime,
  filename = "runtime_candidate",
  directory = candidate_root,
  width = 10.0,
  height = 5.0
)


# ==============================================================================
# 20. FIGURE MANIFEST
# ==============================================================================

figure_files <- list.files(
  figure_root,
  pattern = "\\.(png|pdf)$",
  recursive = TRUE,
  full.names = TRUE
)

figure_manifest <- tibble(
  relative_path = sub(
    paste0(
      "^",
      gsub(
        "\\\\",
        "/",
        figure_root
      ),
      "/?"
    ),
    "",
    gsub(
      "\\\\",
      "/",
      figure_files
    )
  ),
  
  format = tools::file_ext(
    figure_files
  )
)


write_csv(
  figure_manifest,
  file.path(
    data_root,
    "figure_manifest.csv"
  )
)


# ==============================================================================
# 21. README
# ==============================================================================

readme_lines <- c(
  "MIDUS TWO-PHASE SIMULATION — FIGURE LIBRARY",
  "",
  "Key changes:",
  "  - relative precision now uses one common -25% to 450% scale across all scenarios",
  "  - relative precision includes a dotted zero reference line",
  "  - enlarged white mask behind AIPW so connecting lines do not show through",
  "  - common power scale changed to 20%-100% with the 80% reference retained",
  "  - runtime method order explicitly locked to Naive, CCA, FCS-MI, JM-MI, IPW, AIPW",
  "  - Phase-2 sampling fraction displayed with phi_2 to avoid conflict with selection-model eta notation",
  "  - AIPW uses native shape 10 (open circle-plus) in both panels and legend",
  "  - null rejection-rate axis explicitly labels 4%, 5%, and 6%",
  "  - EmpSE and ModSE abbreviations standardized in SE figures",
  "  - runtime methods displayed in the same top-to-bottom order as the manuscript figures",
  "  - Dark2 method colors",
  "  - five white-filled open symbols plus native plotting symbol 10 for AIPW",
  "  - thicker method lines and darker text for manuscript-ready output",
  "  - one-row method legend",
  "  - equally spaced Phase-1 sample-size design levels",
  "  - two-line Phase-2 facet strips using phi_2",
  "  - PDF-safe plotmath theta and phi labels",
  "  - cairo_pdf output",
  "  - fixed y-axis scales by performance measure (bias: -0.01 to 0.025)",
  "  - no shaded acceptable-performance bands",
  "  - reference boundary lines for SE ratio, coverage, and Type I error",
  "  - redesigned horizontal runtime plot with explicit log-scale tick labels",
  "",
  "01_full_factorial/",
  "  Complete R_A^2 x R_Y^2 figure library.",
  "",
  "02_manuscript_candidates/",
  "  Weak, MIDUS-like, and strong auxiliary-information candidate figures.",
  "",
  "03_final_manuscript/",
  "  Intentionally left empty until figures are selected.",
  "",
  "04_diagnostics/",
  "  Redesigned runtime figure; failure rates are stored as data only because all were zero.",
  "",
  "results/figure_data/",
  "  Data underlying figures plus a figure manifest."
)


writeLines(
  readme_lines,
  file.path(
    figure_root,
    "README_FIGURE_LIBRARY.txt"
  )
)


# ==============================================================================
# 22. FINAL MESSAGE
# ==============================================================================

cat(
  "\n",
  paste(
    rep(
      "=",
      72
    ),
    collapse = ""
  ),
  "\n",
  "FIGURE GENERATION COMPLETE\n",
  paste(
    rep(
      "=",
      72
    ),
    collapse = ""
  ),
  "\n\n",
  
  "Full factorial figures:\n  ",
  full_root,
  "\n\n",
  
  "Manuscript candidates:\n  ",
  candidate_root,
  "\n\n",
  
  "Manuscript figures used in the paper are stored under results/manuscript/.\n",
  final_root,
  "\n\n",
  
  "Diagnostics:\n  ",
  diagnostic_root,
  "\n\n",
  
  "Figure data:\n  ",
  data_root,
  "\n\n",
  
  "PNG/PDF files generated: ",
  nrow(
    figure_manifest
  ),
  "\n",
  sep = ""
)
