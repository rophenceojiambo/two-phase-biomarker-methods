################################################################################
# hpc/install_required_packages.R
#
# Run on a COMPUTE NODE after:
#   module load r/4.5.1
# and after R_LIBS_USER has been set to scratch.
################################################################################

lib <- Sys.getenv("R_LIBS_USER")

if (!nzchar(lib)) {
  stop("R_LIBS_USER is not set.")
}

dir.create(
  lib,
  recursive = TRUE,
  showWarnings = FALSE
)

.libPaths(
  unique(
    c(
      lib,
      .libPaths()
    )
  )
)

packages <- c(
  # Core simulation and summaries
  "dplyr",
  "tidyr",
  "purrr",
  "readr",
  "mvtnorm",
  "mice",
  "jomo",
  "broom",
  "ggplot2",
  "scales",
  "rsimsum",
  "generics",
  "tibble",
  # Real-world manuscript outputs
  "patchwork",
  "gtsummary",
  "flextable",
  "officer",
  "Hmisc",
  "smd"
)

installed <- rownames(
  installed.packages()
)

missing <- setdiff(
  packages,
  installed
)

cat(
  "R version:\n",
  R.version.string,
  "\n\nUser library:\n",
  lib,
  "\n\n",
  sep = ""
)

if (length(missing) > 0) {

  cat(
    "Installing:\n",
    paste(
      missing,
      collapse = "\n"
    ),
    "\n\n",
    sep = ""
  )

  install.packages(
    missing,
    repos = "https://cloud.r-project.org",
    lib = lib,
    Ncpus = 4L
  )

} else {
  cat("All required packages are already installed.\n")
}

# Verify all required packages load.
failed <- packages[
  !vapply(
    packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(failed) > 0) {
  stop(
    "These packages still cannot be loaded: ",
    paste(
      failed,
      collapse = ", "
    )
  )
}

cat("\nAll required packages are available.\n\n")
print(
  data.frame(
    package = packages,
    version = vapply(
      packages,
      function(pkg) {
        as.character(
          packageVersion(pkg)
        )
      },
      character(1)
    )
  ),
  row.names = FALSE
)

project_dir <- Sys.getenv(
  "SIM_PROJECT_DIR",
  unset = "."
)

results_dir <- file.path(
  project_dir,
  "results"
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    results_dir,
    "sessionInfo_setup.txt"
  )
)
