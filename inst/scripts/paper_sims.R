#!/usr/bin/env Rscript

library(docopt)

doc_string <- "

Usage:
./paper_sims nrows <n_rows> ncols <n_cols> nweightvecs <num_weight_vec>
./paper_sims [options]

Options:
  -h, --help    display this page
  -x, --usage   show some examples
"


##------------------------------------

library(matchfinder)
library(snowflaker)
library(xgboost)

##------------------------------------

treat_funcs <- paper_treatment_functions(target_mean = 0.425)
mu_funcs <- paper_mean_functions()

n_rows <- arguments[["n_rows"]]
n_cols <- arguments[["n_cols"]]
num_weight_vec <- arguments[["num_weight_vec"]]

##------------------------------------

result_list <- list()
full_frame_list <- list()

f_i <- 1L
m_j <- 1L
for (treat_func_name in names(treat_funcs)) {
    print(paste0(treat_func_name, " ", f_i, "/", length(treat_funcs)))

    mu_list <- list()
    for (mu_func_name in names(mu_funcs)) {
        print(paste0(mu_func_name, " ", m_j, "/", length(mu_funcs)))

        sim_res <- compute_sim_result(
            x_generator = default_x_generator,
            treat_prob_generator = treat_funcs[[treat_func_name]],
            mean_generator = mu_funcs[[mu_func_name]],
            error_generator = default_error_generator,
            n_sink_gen = n_sink_generator(),
            match_method = "with_replacement",
            n_rows = n_rows,
            n_cols = n_cols
            num_weight_vectors = num_weight_vec,
            silent = TRUE)

        mu_list[[mu_func_name]] <- sim_res

        sim_frame <- reshape_list_of_sims(
            list_of_sims = list(sim_res),
            treat_model_name = treat_func_name,
            mu_model_name = mu_func_name,
            n_rows = n_rows,
            n_cols = n_cols,
            num_weight_vectors = num_weight_vec)

        full_frame_list[[length(full_frame_list) + 1L]] <- sim_frame
    }

    result_list[[treat_func_name]] <- mu_list
}

full_frame_output <- do.call(rbind, full_frame_list)

output_list <- list(
    result_list = result_list,
    full_frame_output = full_frame_output)

##------------------------------------
## send to S3

uuid <- paste(sample(c(letters, 0:9), 30), collapse = "")
s3_file <- glue::glue("
sim_res/nrow_{n_rows}/ncol_{n_cols}/nweightvec_{num_weight_vec}/\\
sim_{uuid}.rds
",
n_rows = n_rows,
n_cols = n_cols,
num_weight_vec = num_weight_vec,
uuid = uuid)

temp_file <- tempfile(fileext = ".rds")
saveRDS(output_list, file = temp_file)
on.exit({unlink(temp_file)})
snowflaker::local_to_s3(local_file = temp_file,
                        s3_file = paste0("r_sims/", s3_file),
                        s3_bucket = "via-nyc-data-science")

##------------------------------------
