context("testing bipartite_matches")


test_that("testing errors", {
    ## need dist_mat, treat_vec
    expect_error(bipartite_matches())
    expect_error(bipartite_matches(dist_mat = matrix(1:4)))

    ## dist_mat must be matrix
    expect_error(bipartite_matches(dist_mat = 1:100,
                                   treat_vec = rep(c(0, 1), times = c(10, 10))))
    expect_error(bipartite_matches(dist_mat = matrix(1:100, 10, 10),
                                   treat_vec = rep(c(0, 1), times = c(10, 10))),
                 NA)
    expect_error(bipartite_matches(
        dist_mat = matrix(1:100, 10, 10),
        treat_vec = rep(c(0, 1), times = c(10, 10)),
        n_sinks = 0),
        NA)

    ## treat_vec length must be rows + cols of dist_mat, and must be logical
    ## or 0, 1
    expect_error(bipartite_matches(dist_mat = matrix(1:100, 10, 10),
                                   treat_vec = rep(c(0, 1), times = c(9, 9))))
    expect_error(bipartite_matches(dist_mat = matrix(1:100, 10, 10),
                                   treat_vec = rep(c("0", "1"), times = c(10, 10))))
    expect_error(bipartite_matches(dist_mat = matrix(1:100, 10, 10),
                                   treat_vec = rep(c(0, 2), times = c(10, 10))))
    expect_error(bipartite_matches(dist_mat = matrix(1:100, 10, 10),
                                   treat_vec = rep(c(0, 1), times = c(10, 10))),
                 NA)
    expect_error(bipartite_matches(dist_mat = matrix(1:100, 10, 10),
                                   treat_vec = rep(c(FALSE, TRUE), times = c(10, 10))),
                 NA)

    expect_error(bipartite_matches(dist_mat = matrix(1:100, 10, 10),
                                   treat_vec = rep(c(FALSE, TRUE),
                                                   times = c(10, 10)),
                                   match_method = c("greedy", "with_replacement")))
    expect_error(bipartite_matches(dist_mat = matrix(1:100, 10, 10),
                                   treat_vec = rep(c(FALSE, TRUE),
                                                   times = c(10, 10)),
                                   match_method = c("something good")))
    expect_error(bipartite_matches(dist_mat = matrix(1:100, 10, 10),
                                   treat_vec = rep(c(FALSE, TRUE),
                                                   times = c(10, 10)),
                                   match_method = "with_replacement"),
                 NA)

    ## if you set a tol_val, you must be doing optimal
    expect_error(bipartite_matches(dist_mat = matrix(1:100, 10, 10),
                                   treat_vec = rep(c(FALSE, TRUE),
                                                   times = c(10, 10)),
                                   match_method = "with_replacement",
                                   tol_val = 0.01))
    expect_error(bipartite_matches(dist_mat = matrix(1:100, 10, 10),
                                   treat_vec = rep(c(FALSE, TRUE),
                                                   times = c(10, 10)),
                                   match_method = "optimal",
                                   tol_val = 0.01),
                 NA)
})


test_that("testing with replacement", {
    dist_mat <- matrix(1, 10, 10)
    dist_mat[, 1] <- 0.1
    treat_vec <- rep(1, 20)
    treat_vec[c(6:10, 16:20)] <- 0

    wr_matches <- bipartite_matches(dist_mat = dist_mat,
                                    treat_vec = treat_vec,
                                    match_method = "with_replacement")[[1]]
    expect_equal(wr_matches$treat_index_within, 1:10)
    expect_equal(wr_matches$treat_index, c(1:5, 11:15))
    ## they all get matched to the first
    expect_equal(wr_matches$control_index_within, rep(1, 10))
    expect_equal(wr_matches$control_index, rep(6, 10))
    expect_equal(wr_matches$distance, rep(0.1, 10))

    ## slightly more interesting
    dist_mat <- matrix(runif(400, 5, 10), 20, 20)
    best_index <- fixed_sample(1:20, replace = TRUE)
    dist_mat[cbind(1:20, best_index)] <- runif(20, 1, 4)
    treat_vec <- rep(1, 40)
    treat_vec[11:30] <- 0

    wr_matches <- bipartite_matches(dist_mat = dist_mat,
                                    treat_vec = treat_vec,
                                    match_method = "with_replacement")[[1]]
    ## boring:
    expect_equal(wr_matches$treat_index_within, 1:20)
    expect_equal(wr_matches$treat_index, c(1:10, 31:40))

    ## they get matched to the best index
    expect_equal(wr_matches$control_index_within, best_index)
    expect_equal(wr_matches$control_index, c(11:30)[best_index])
    expect_true(sum(wr_matches$distance) >= 1 * 20 &&
                sum(wr_matches$distance) <= 4 * 20)
})


test_that("testing greedy and optimal", {
    ##------------------------------------
    ## first greedy

    dist_mat <- matrix(1, 10, 10)
    dist_mat[, 1] <- 0.1
    treat_vec <- rep(1, 20)
    treat_vec[c(6:10, 16:20)] <- 0

    gr_matches <- bipartite_matches(dist_mat = dist_mat,
                                    treat_vec = treat_vec,
                                    match_method = "greedy")[[1]]
    expect_equal(gr_matches$treat_index_within, 1:10)
    expect_equal(gr_matches$treat_index, c(1:5, 11:15))

    ## who knows who gets picked, but we know it's unique
    expect_equal(sort(gr_matches$control_index_within), 1:10)
    expect_equal(sort(gr_matches$control_index), c(6:10, 16:20))
    ## only one of them gets the nice 0.1
    expect_equal(sort(gr_matches$distance), c(0.1, rep(1, 9)))

    ## slightly more interesting
    dist_mat <- matrix(runif(400, 5, 10), 20, 20)
    ## force best to not be unique
    best_index <- fixed_sample(1:10, size = 20, replace = TRUE)
    dist_mat[cbind(1:20, best_index)] <- 4
    treat_vec <- rep(1, 40)
    treat_vec[11:30] <- 0

    gr_matches <- bipartite_matches(dist_mat = dist_mat,
                                    treat_vec = treat_vec,
                                    match_method = "greedy")[[1]]
    ## boring:
    expect_equal(gr_matches$treat_index_within, 1:20)
    expect_equal(gr_matches$treat_index, c(1:10, 31:40))

    ## not the best!
    expect_true(sum(gr_matches$distance) > 4 * 10 + 5 * 10)

    ## example where greedy isn't optimal
    ## (see optimal_match tests for successful case)

    dist_mat <- rbind(
        c(1, 10, 10),
        c(2, 100, 100)
    )
    treat_vec <- c(1, 0, 0, 0, 1)

    ## rows are randomized, so 50% chance of "messing it up"
    ## 2^(-40) works for me
    gr_list <- list()
    for (j in 1:40) {
        gr_list[[j]] <- bipartite_matches(dist_mat = dist_mat,
                                          treat_vec = treat_vec,
                                          match_method = "greedy")[[1]]
    }
    dist_totals <- lapply(gr_list, function(x) {
        sum(x[["distance"]])
    })
    worst_dist <- max(unlist(dist_totals))
    expect_equal(worst_dist, 101)

    ##------------------------------------
    ## now optimal

    ## first let's revisit the issue above
    optimal_test <- bipartite_matches(dist_mat = dist_mat,
                                      treat_vec = treat_vec,
                                      match_method = "optimal")[[1]]
    expect_equal(sum(optimal_test[["distance"]]), 12)

    ## random numbers
    treat_vec <- rep(1, 20)
    treat_vec[c(6:10, 16:20)] <- 0
    rand_result <- list()
    greedy_stuff <- list()
    for (j in 1:5) {
        rand_mat <- matrix(runif(100), 10, 10)
        rand_result[[j]] <- bipartite_matches(dist_mat = rand_mat,
                                              treat_vec = treat_vec,
                                              match_method = "optimal")[[1]]
        greedy_stuff[[j]] <- bipartite_matches(dist_mat = rand_mat,
                                               treat_vec = treat_vec,
                                              match_method = "greedy")[[1]]
    }

    optimal_wins <- Map(function(x, y){
        sum(x$distance) > sum(y$distance)
    }, greedy_stuff, rand_result)

    ## at least two...
    expect_true(mean(unlist(optimal_wins)) > 0.3)
})


test_that("testing with sinks", {
    default <- bipartite_matches(dist_mat = matrix(1:100, 10, 10),
                                 treat_vec = rep(c(0, 1), times = c(10, 10)),
                                 match_method = "with_replacement")[[1]]
    zero_sinks <- bipartite_matches(
        dist_mat = matrix(1:100, 10, 10),
        treat_vec = rep(c(0, 1), times = c(10, 10)),
        match_method = "with_replacement",
        n_sinks = 0)[[1]]

    ## controversial difference: not type consistent, as default assumes you
    ## don't want a list... but even a single value does

    default[["num_sinks"]] <- 0

    expect_equal(default, zero_sinks)

    ## similar:
    multiple_sinks <- bipartite_matches(
        dist_mat = matrix(1:100, 10, 10),
        treat_vec = rep(c(0, 1), times = c(10, 10)),
        match_method = "with_replacement",
        n_sinks = c(0, 2, 5))

    expect_equal(multiple_sinks[[1]][["distance"]],
                 1:10)
    expect_equal(multiple_sinks[[2]][["distance"]],
                 1:8)
    expect_equal(multiple_sinks[[3]][["distance"]],
                 1:5)

    expect_equal(multiple_sinks[[1]], default)

    ##------------------------------------
    ## and optimal

    treat_vec <- rep(1, 20)
    treat_vec[c(6:10, 16:20)] <- 0
    rand_result <- list()
    for (j in 1:20) {
        rand_mat <- matrix(runif(100), 10, 10)
        rand_result[[j]] <- bipartite_matches(dist_mat = rand_mat,
                              treat_vec = treat_vec,
                              match_method = "optimal",
                              n_sinks = c(0, 4))
    }

    long_worse <- lapply(rand_result, function(x){
        full_distance <- x[[1]][["distance"]]
        short_distance <- x[[2]][["distance"]]

        return(c(sum(sort(full_distance)[1:length(short_distance)]) >
                 sum(short_distance),
                 sum(sort(full_distance)[1:length(short_distance)]) >=
                 sum(short_distance)
                 ))
    })

    strictly_better <- unlist(lapply(long_worse, function(x){
        x[1]
    }))
    at_least_better <- unlist(lapply(long_worse, function(x){
        x[2]
    }))
    expect_true(all(at_least_better))
    expect_true(mean(strictly_better) > 0.4)
})
