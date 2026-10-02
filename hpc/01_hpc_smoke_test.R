################################################################################
# hpc/01_hpc_smoke_test.R
#
# Three representative datasets using final production settings.
# This confirms the Torch R environment and all six analysis methods before
# the large Slurm array is launched.
################################################################################

source("00_config.R")
source("02_dgm.R")
source("03_methods.R")

library(dplyr)
library(readr)

if (!file.exists(calibration_rds)) {
  stop(
    "Calibration object not found: ",
    calibration_rds,
    "\nRun 01_calibrate_midus.R first."
  )
}

calibration <- readRDS(
  calibration_rds
)

smoke_pass_file <- file.path(
  results_dir,
  "HPC_SMOKE_TEST_PASS.txt"
)

if (file.exists(smoke_pass_file)) {
  unlink(smoke_pass_file)
}

scenarios <- data.frame(
  scenario = c(
    "MIDUS-like",
    "Stress",
    "Worst-N benchmark"
  ),
  N = c(
    800L,
    500L,
    1500L
  ),
  phase2_fraction = c(
    0.65,
    0.25,
    0.25
  ),
  r2_a_marker = c(
    0.05,
    0.10,
    0.10
  ),
  r2_y_marker = c(
    0.05,
    0.10,
    0.10
  ),
  theta = c(
    0.15,
    0.15,
    0.15
  ),
  stringsAsFactors = FALSE
)

run_timed <- function(
  method_name,
  expr
) {

  tm <- system.time({

    ans <- tryCatch(
      eval.parent(
        substitute(expr)
      ),
      error = function(e) {
        failed_result(
          method_name,
          conditionMessage(e)
        )
      }
    )
  })

  ans$elapsed_seconds <- unname(
    tm["elapsed"]
  )

  ans
}

out <- list()

for (
  s in seq_len(
    nrow(scenarios)
  )
) {

  sc <- scenarios[
    s,
    ,
    drop = FALSE
  ]

  dgm <- build_dgm_parameters(
    calibration = calibration,
    r2_a_marker = sc$r2_a_marker,
    r2_y_marker = sc$r2_y_marker,
    theta = sc$theta,
    phase2_fraction = sc$phase2_fraction,
    selection = "mar",
    selection_calibration_n = 50000L,
    selection_seed = 101000L + s,
    marker_error = "mvn"
  )

  set.seed(
    102000L + s
  )

  dat <- generate_two_phase_data(
    N = sc$N,
    calibration = calibration,
    dgm_parameters = dgm
  )$observed

  methods <- bind_rows(
    run_timed(
      "Naive",
      fit_naive(
        dat
      )
    ),
    run_timed(
      "CCA",
      fit_cca(
        dat,
        calibration$marker_names
      )
    ),
    run_timed(
      "FCS-MI",
      fit_fcs_mi(
        dat = dat,
        marker_names =
          calibration$marker_names,
        nimp = nimp_primary,
        maxit =
          mice_maxit_primary,
        method =
          mice_method_primary
      )
    ),
    run_timed(
      "JM-MI",
      fit_joint_mi(
        dat = dat,
        marker_names =
          calibration$marker_names,
        nimp = nimp_primary,
        nburn =
          jomo_nburn_primary,
        nbetween =
          jomo_nbetween_primary
      )
    ),
    run_timed(
      "IPW",
      fit_ipw(
        dat,
        calibration$marker_names
      )
    ),
    run_timed(
      "AIPW",
      fit_aipw(
        dat = dat,
        marker_names =
          calibration$marker_names,
        selection_include_y =
          TRUE,
        marker_include_y =
          TRUE,
        method_label =
          "AIPW"
      )
    )
  ) %>%
    mutate(
      scenario = sc$scenario,
      N = sc$N,
      n_phase2 = sum(
        dat$phase2
      ),
      realized_phase2_fraction =
        mean(
          dat$phase2
        ),
      theta_true =
        sc$theta,
      jomo_nburn =
        jomo_nburn_primary,
      jomo_nbetween =
        jomo_nbetween_primary
    )

  out[[s]] <- methods
}

out <- bind_rows(
  out
)

write_csv(
  out,
  file.path(
    results_dir,
    "hpc_smoke_test_results.csv"
  )
)


runtime_projection <- out %>%
  group_by(
    scenario,
    N
  ) %>%
  summarise(
    observed_seconds_one_dataset =
      sum(
        elapsed_seconds,
        na.rm = TRUE
      ),
    projected_hours_100_repetitions =
      100 *
        observed_seconds_one_dataset /
        3600,
    .groups = "drop"
  )

write_csv(
  runtime_projection,
  file.path(
    results_dir,
    "hpc_smoke_runtime_projection.csv"
  )
)

smoke_display <- out %>%
  select(
    scenario,
    method,
    estimate,
    se,
    status,
    elapsed_seconds,
    n_phase2
  )

print(
  as.data.frame(smoke_display),
  row.names = FALSE
)

expected_methods <- c(
  "Naive",
  "CCA",
  "FCS-MI",
  "JM-MI",
  "IPW",
  "AIPW"
)

expected_rows <- nrow(scenarios) * length(expected_methods)

if (nrow(out) != expected_rows) {
  stop(
    "HPC smoke test returned ",
    nrow(out),
    " rows; expected ",
    expected_rows,
    "."
  )
}

if (
  anyNA(out$status) ||
    any(out$status != "ok")
) {
  stop(
    "At least one method failed the HPC smoke test."
  )
}

method_counts <- out %>%
  count(
    scenario,
    method,
    name = "n"
  )

if (
  nrow(method_counts) != expected_rows ||
    any(method_counts$n != 1L) ||
    !setequal(unique(out$method), expected_methods)
) {
  stop(
    "HPC smoke test did not return exactly one result for every ",
    "scenario-method combination."
  )
}

writeLines(
  c(
    "HPC SMOKE TEST: PASS",
    paste0("Completed: ", Sys.time()),
    paste0("nimp: ", nimp_primary),
    paste0("jomo_nburn: ", jomo_nburn_primary),
    paste0("jomo_nbetween: ", jomo_nbetween_primary)
  ),
  smoke_pass_file
)

cat(
  "\nHPC smoke test PASSED.\n",
  "JM-MI settings: nburn=",
  jomo_nburn_primary,
  ", nbetween=",
  jomo_nbetween_primary,
  ", nimp=",
  nimp_primary,
  "\nPASS marker: ",
  smoke_pass_file,
  "\n",
  sep = ""
)
