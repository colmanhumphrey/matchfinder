library(matchfinder)
library(xgboost)

n_row_values <- c(500L, 1000L, 2000L, 4000L)
n_col_values <- c(5L, 10L, 20L)
num_weight_vecs <- 300L

##------------------------------------

treat_funcs <- paper_treatment_functions(treat_mean = 0.425)
mu_funcs <- paper_mean_functions()

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
