#!/usr/bin/env Rscript

##------------------------------------

library(glue)

##------------------------------------

if (dir.exists("./jobs")) {
    rel_files <- list.files("./jobs", full.names = TRUE)
    unlink(rel_files)
} else {
    dir.create("./jobs")
}

##------------------------------------

job_yaml <- function(n_rows, n_cols, n_weight_vecs) {
    glue::glue("apiVersion: batch/v2
kind: Job
metadata:
  name: sim-rows_{n_rows}-cols_{n_cols}-weightvecs_{n_weight_vecs}
  labels:
    jobgroup: sim_runs
spec:
  template:
    metadata:
      name: sim_runs
      labels:
        jobgroup: sim_runs
    spec:
      containers:
      - name: c
        image: aws-account.r_runner_etc
        command: [\"sh\", \\
\"-c\", \\
\"./run_the_stuff nrows {n_rows} ncols {n_cols} nweightvecs {n_weight_vecs}\"]
      restartPolicy: Never",
n_rows = n_rows,
n_cols = n_cols,
n_weight_vecs = n_weight_vecs)
}

##------------------------------------

n_row_values <- c(1000L, 2000L, 4000L)
n_col_values <- c(5L, 10L, 20L)
num_weight_vecs <- c(150L, 300L)

for (n_rows in n_row_values) {
    for (n_cols in n_col_values) {
        for (n_weight_vecs in num_weight_vecs) {
            file_name <- glue::glue("
./jobs/job-rows_{n_rows}-cols_{n_cols}-weightvecs_{n_weight_vecs}.yaml
",
n_rows = n_rows,
n_cols = n_cols,
n_weight_vecs = n_weight_vecs)

            yaml_template <- job_yaml(n_rows = n_rows,
                                      n_cols = n_cols,
                                      n_weight_vecs = n_weight_vecs)
            writeLines(as.character(yaml_template),
                       con = file_name)
        }
    }
}

##------------------------------------
