library(matchfinder)



if (TRUE) {
    start_time <- Sys.time()
    sim_res <- compute_sim_result(
        treat_prob_generator =
            paper_treatment_functions()[["logistic_treat_prob"]],
        mean_generator =
            paper_mean_functions()[["linear_mu"]],
        num_weight_vectors = 300L,
        n_rows = 4000L,
        n_cols = 20L)
    end_time <- Sys.time()
}
min_diff <- as.numeric(difftime(end_time, start_time, units = "minutes"))

df <- reshape_list_of_sims(list(sim_res),
                     "logistic_treat",
                     "linear_mu",
                     500,
                     5)

##------------------------------------

parallel_sim(x_generator = default_x_generator,
             )
