################################################################################
# hpc/02_check_primary_completion.R
#
# Checks which of the 6,480 expected chunk files are complete.
# No simulation is rerun.
################################################################################

source("00_config.R")

expected <- data.frame(
  task_id = seq_len(
    n_primary_array_tasks
  )
)

expected$scenario_id <- (
  (expected$task_id - 1L) %/%
    chunks_per_scenario
) + 1L

expected$chunk_id <- (
  (expected$task_id - 1L) %%
    chunks_per_scenario
) + 1L

expected$file <- file.path(
  chunk_dir,
  sprintf(
    "scenario_%03d",
    expected$scenario_id
  ),
  sprintf(
    "scenario_%03d_chunk_%02d.rds",
    expected$scenario_id,
    expected$chunk_id
  )
)

expected$complete <- file.exists(
  expected$file
)

missing <- expected[
  !expected$complete,
  ,
  drop = FALSE
]

utils::write.csv(
  expected,
  file.path(
    results_dir,
    "primary_task_completion.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  missing,
  file.path(
    results_dir,
    "primary_missing_tasks.csv"
  ),
  row.names = FALSE
)

writeLines(
  if (
    nrow(
      missing
    ) == 0L
  ) {
    ""
  } else {
    paste(
      missing$task_id,
      collapse = ","
    )
  },
  con = file.path(
    results_dir,
    "primary_missing_array_ids.txt"
  )
)

cat(
  "Completed array tasks: ",
  sum(
    expected$complete
  ),
  " / ",
  nrow(
    expected
  ),
  "\nMissing array tasks: ",
  nrow(
    missing
  ),
  "\n",
  sep = ""
)

if (
  nrow(
    missing
  ) == 0L
) {
  cat(
    "PRIMARY ARRAY COMPLETE.\n"
  )
} else {
  cat(
    "Missing IDs written to:\n",
    file.path(
      results_dir,
      "primary_missing_array_ids.txt"
    ),
    "\n",
    sep = ""
  )
}
