context("testing simulation.R")


test_that("testing generate_simulation_input", {
    ## using the example functions

    example_list <- generate_simulation_input(n_rows = 1000L,
                                              n_cols = 8L)

    expect_true(setequal(names(example_list),
                         c("x_mat", "treat_vec", "y_vec")))

    expect_equal(dim(example_list[["x_mat"]]),
                 c(1000L, 8L))

    expect_equal(length(example_list[["treat_vec"]]),
                 length(example_list[["y_vec"]]))
    expect_equal(length(example_list[["treat_vec"]]),
                 nrow(example_list[["x_mat"]]))

    expect_true(all(example_list[["treat_vec"]] %in% c(0, 1)))
})


test_that("testing paper_treatment_functions", {
    target_value <- runif(1L, 0.3, 0.7)

    treat_func_list <- paper_treatment_functions(target_mean = target_value)

    x_mat <- default_x_generator(n_rows = 1000L,
                                 n_cols = 4L)

    treat_vectors <- lapply(treat_func_list, function(func){
        func(x_mat)
    })

    tv_means <- unlist(lapply(treat_vectors, mean))

    expect_true(sum(abs(tv_means - target_value)) < 0.00001)

    expect_true(all(0 < unlist(treat_vectors) & unlist(treat_vectors) < 1))

    ##----------------

    expect_true(all(treat_vectors[["constant_treat_prob"]] == target_value))

    sparse_table <- table(treat_vectors[["sparse_treat_prob"]])
    expect_equal(unname(as.integer(sparse_table)), c(500L, 500L))
})



test_that("testing paper_mean_functions", {
    mu_func_list <- paper_mean_functions()

    x_mat <- default_x_generator(n_rows = 1000L,
                                 n_cols = 4L)

    mu_vectors <- lapply(mu_func_list, function(func){
        func(x_mat)
    })

    mv_means <- unlist(lapply(mu_vectors, mean))

    expect_true(sum(abs(mv_means - 0)) < 0.00001)

    ##----------------

    expect_true(all(mu_vectors[["constant_mu"]] == 0))

    sign_table <- table(mu_vectors[["sign_mu"]])
    expect_equal(unname(as.integer(sign_table)), c(500L, 500L))
})


test_that("testing n_sink_generator", {
    expect_error(n_sink_generator(start_frac = -1))
    expect_error(n_sink_generator(end_frac = 100))
    expect_error(n_sink_generator(length_out = 0))

    expect_error(n_sink_generator(start_frac = c(0, 0.2)))
    expect_error(n_sink_generator(end_frac = c(0.9, 1)))
    expect_error(n_sink_generator(length_out = c(10, 100)))

    expect_error(n_sink_generator(start_frac = 0.5,
                                  end_frac = 0.6,
                                  length_out = 1L))
    expect_error(n_sink_generator(start_frac = 0.6,
                                  end_frac = 0.6,
                                  length_out = 1L),
                 NA)

    expect_error(n_sink_generator(start_frac = 0.6,
                                  end_frac = 0.5,
                                  length_out = 2L))
    expect_error(n_sink_generator(start_frac = 0.5,
                                  end_frac = 0.5,
                                  length_out = 2L))
    expect_error(n_sink_generator(start_frac = 0.5,
                                  end_frac = 0.6,
                                  length_out = 2L),
                 NA)

    ##------------------------------------

    n_sink_func <- n_sink_generator(start_frac = 0,
                                    end_frac = 0.6,
                                    length_out = 7L)

    treat_vec <- c(rep(0L, 100),
                   rep(1L, 50))

    n_sink_vector <- n_sink_func(treat_vec)

    expect_equal(n_sink_vector,
                 c(0L, 5L, 10L, 15L, 20L, 25L, 30L))

    treat_vec <- c(rep(0L, 100),
                   rep(1L, 500))

    n_sink_vector <- n_sink_func(treat_vec)

    expect_equal(n_sink_vector,
                 c(0L, 50L, 100L, 150L, 200L, 250L, 300L))
})


test_that("testing compute_sim_result and sim_unwrap", {
    sim_results <- compute_sim_result(
        x_generator = default_x_generator,
        treat_prob_generator =
            paper_treatment_functions()[["logistic_treat_prob"]],
        mean_generator = paper_mean_functions()[["linear_mu"]],
        error_generator = default_error_generator,
        n_sink_gen = n_sink_generator(
            start_frac = 0,
            end_frac = 0.4,
            length_out = 2),
        match_method = "with_replacement",
        n_rows = 100L,
        n_cols = 3L,
        num_weight_vectors = 5L,
        silent = TRUE)

    expect_equal(names(sim_results),
                 c("naive_est", "propensity_results",
                   "mahal_results", "weighted_results"))
    expect_true(length(sim_results[["naive_est"]]) == 1L &&
                abs(sim_results[["naive_est"]] - 1) < 5)

    expect_equal(lengths(sim_results[["propensity_results"]]),
                 c(2L, 2L))

    expect_equal(lengths(sim_results[["mahal_results"]]),
                 c(2L, 2L))

    expect_equal(lengths(sim_results[["weighted_results"]]),
                 c(4L, 4L))

    flat_sims <- reshape_list_of_sims(
        list(sim_results, sim_results),
        treat_model_name = "logistic",
        mu_model_name = "linear",
        n_rows = 100L,
        n_cols = 3L,
        p_cut = c(1, 0.9))

    expect_equal(dim(flat_sims),
                 c(4L, 12L))

    expect_equal(unlist(lapply(flat_sims, typeof)),
                 c(treat_model = "character",
                   mu_model = "character",
                   p_cut = "double",
                   p_brier = "double",
                   raw_brier = "double",
                   n_rows = "integer",
                   n_cols = "integer",
                   n_sinks = "integer",
                   naive_est = "double",
                   propensity_est = "double",
                   mahal_est = "double",
                   weighted_est = "double"))
})
