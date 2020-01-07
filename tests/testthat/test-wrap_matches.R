context("testing wrap_matches.R")


test_that("testing simple_sink_wrap", {
    ##------------------------------------
    ## bipartite

    x_mat <- cbind(
        rnorm(20),
        runif(20))
    treat_vec <- (1L:20L) %in% fixed_sample(1L:20L, size = 10L)
    weight_vec <- c(3, 1)
    dist_mat <- weighted_mahal(x_mat,
                               cov(x_mat),
                               weight_vec = weight_vec,
                               treat_vec = treat_vec)
    wr_zero <- with_replacement_match(dist_mat,
                                      treat_vec = treat_vec)
    worst_three <- rank(wr_zero[["distance"]], ties.method = "first") >= 8L
    wr_three <- lapply(wr_zero, function(x){
        x[!worst_three]
    })

    wr_list <- simple_sink_wrap(wr_zero, n_sinks = c(0L, 3L))

    zero_list <- wr_list[["0"]]
    three_list <- wr_list[["3"]]
    expect_equal(zero_list[["num_sinks"]], 0L)
    expect_equal(three_list[["num_sinks"]], 3L)
    zero_list[["num_sinks"]] <- NULL
    three_list[["num_sinks"]] <- NULL

    expect_equal(zero_list, wr_zero)
    expect_equal(three_list, wr_three)
    ##----------------
    greedy_zero <- greedy_match(dist_mat,
                                treat_vec = treat_vec)
    worst_three <- rank(greedy_zero[["distance"]], ties.method = "first") >= 8L
    greedy_three <- lapply(greedy_zero, function(x){
        x[!worst_three]
    })

    greedy_list <- simple_sink_wrap(greedy_zero, n_sinks = c(0L, 3L))

    zero_list <- greedy_list[["0"]]
    three_list <- greedy_list[["3"]]
    expect_equal(zero_list[["num_sinks"]], 0L)
    expect_equal(three_list[["num_sinks"]], 3L)
    zero_list[["num_sinks"]] <- NULL
    three_list[["num_sinks"]] <- NULL

    expect_equal(zero_list, greedy_zero)
    expect_equal(three_list, greedy_three)

    ##------------------------------------
    ## nonbipartite

    x_mat <- cbind(
        rnorm(20),
        runif(20))
    dist_mat <- weighted_mahal(ranked_x(x_mat, 2),
                               covariance_with_ranks(x_mat, 2))
    diag(dist_mat) <- Inf

    tolerance_vec <- runif(20L)

    ##----------------
    wr_zero <- with_replacement_nbp_match(dist_mat,
                                          tolerance_vec = tolerance_vec,
                                          keep_all = FALSE)
    worst_three <- rank(wr_zero[["distance"]], ties.method = "first") >= 8L
    wr_three <- lapply(wr_zero, function(x){
        x[!worst_three]
    })

    wr_list <- simple_sink_wrap(wr_zero, n_sinks = c(0L, 3L))

    zero_list <- wr_list[["0"]]
    three_list <- wr_list[["3"]]
    expect_equal(zero_list[["num_sinks"]], 0L)
    expect_equal(three_list[["num_sinks"]], 3L)
    zero_list[["num_sinks"]] <- NULL
    three_list[["num_sinks"]] <- NULL

    expect_equal(zero_list, wr_zero)
    expect_equal(three_list, wr_three)
    ##----------------
    greedy_zero <- greedy_nbp_match(dist_mat,
                                    tolerance_vec = tolerance_vec)
    worst_three <- rank(greedy_zero[["distance"]], ties.method = "first") >= 8L
    greedy_three <- lapply(greedy_zero, function(x){
        x[!worst_three]
    })

    greedy_list <- simple_sink_wrap(greedy_zero, n_sinks = c(0L, 3L))

    zero_list <- greedy_list[["0"]]
    three_list <- greedy_list[["3"]]
    expect_equal(zero_list[["num_sinks"]], 0L)
    expect_equal(three_list[["num_sinks"]], 3L)
    zero_list[["num_sinks"]] <- NULL
    three_list[["num_sinks"]] <- NULL

    expect_equal(zero_list, greedy_zero)
    expect_equal(three_list, greedy_three)
})


test_that("testing optimal_sink_wrap", {
    x_mat <- cbind(
        rnorm(20L),
        runif(20L))
    treat_vec <- (1L:20L) %in% fixed_sample(1L:20L, size = 10L)
    dist_mat <- weighted_mahal(ranked_x(x_mat, 2),
                               covariance_with_ranks(x_mat, 2),
                               treat_vec = treat_vec)
    opt_four <- optimal_match(dist_mat,
                              treat_vec = treat_vec,
                              n_sinks = 4,
                              tol_val = 0.0001)
    opt_zero <- optimal_match(dist_mat,
                              treat_vec = treat_vec,
                              n_sinks = 0,
                              tol_val = 0.0001)

    opt_list <- optimal_sink_wrap(dist_mat,
                                  treat_vec = treat_vec,
                                  n_sinks = c(0L, 4L),
                                  tol_val = 0.0001)

    zero_list <- opt_list[["0"]]
    expect_equal(zero_list[["num_sinks"]], 0L)
    zero_list[["num_sinks"]] <- NULL
    expect_equal(zero_list, opt_zero)

    four_list <- opt_list[["4"]]
    expect_equal(four_list[["num_sinks"]], 4L)
    four_list[["num_sinks"]] <- NULL
    expect_equal(four_list, opt_four)

    opt_zero[["num_sinks"]] <- 0L
    expect_equal(optimal_sink_wrap(dist_mat,
                                   treat_vec = treat_vec,
                                   n_sinks = NULL,
                                   tol_val = 0.0001),
                 list("0" = opt_zero))
})


test_that("testing optimal_nbp_sink_wrap", {
    x_mat <- cbind(
        rnorm(20L),
        runif(20L))
    dist_mat <- weighted_mahal(ranked_x(x_mat, 2L),
                               covariance_with_ranks(x_mat, 2L))
    diag(dist_mat) <- Inf

    tolerance_vec <- runif(20L)

    opt_four <- optimal_nbp_match(dist_mat,
                                  tolerance_vec = tolerance_vec,
                                  n_sinks = 4)
    opt_zero <- optimal_nbp_match(dist_mat,
                                  tolerance_vec = tolerance_vec,
                                  n_sinks = 0)

    opt_list <- optimal_nbp_sink_wrap(dist_mat,
                                      tolerance_vec = tolerance_vec,
                                      n_sinks = c(0L, 4L))

    zero_list <- opt_list[["0"]]
    expect_equal(zero_list[["num_sinks"]], 0L)
    zero_list[["num_sinks"]] <- NULL
    expect_equal(zero_list, opt_zero)

    four_list <- opt_list[["4"]]
    expect_equal(four_list[["num_sinks"]], 4L)
    four_list[["num_sinks"]] <- NULL
    expect_equal(four_list, opt_four)

    opt_zero[["num_sinks"]] <- 0L
    expect_equal(optimal_nbp_sink_wrap(dist_mat,
                                       tolerance_vec = tolerance_vec,
                                       n_sinks = NULL),
                 list("0" = opt_zero))
})
