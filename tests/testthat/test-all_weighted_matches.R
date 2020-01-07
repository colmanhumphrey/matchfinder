context("testing all weighted matches")


test_that("testing all_bipartite_matches", {
    ##------------------------------------
    ## no caliper/propens

    rows <- 100L
    num_weight_vecs <- 5L
    x_mat <- cbind(rnorm(rows),
                   runif(rows))
    treat_vec <- (1L:rows) %in% fixed_sample(1L:rows,
                                             45L)
    y_vector <- x_mat[, 1] + x_mat[, 2] + treat_vec * 0.3
    cov_x <- covariance_with_ranks(x_mat)

    some_dist_mat = weighted_mahal(x_mat,
                                   cov_x = cov_x,
                                   weight_vec = c(0.66, 0.33),
                                   treat_vec = treat_vec)

    weight_vecs <- generate_random_weights(prior_weights = c(2, 1),
                                           number_vectors = num_weight_vecs,
                                           minimum_weights = c(0.1, 0.1))

    all_wr_matches <- all_bipartite_matches(x_mat = x_mat,
                                            cov_x = cov(x_mat),
                                            weight_list = weight_vecs,
                                            treat_vec = treat_vec,
                                            match_method = "with_replacement",
                                            n_sinks = c(0L, 4L))

    zero_wr_matches <- lapply(all_wr_matches, function(x){
        x[["0"]]
    })
    four_wr_matches <- lapply(all_wr_matches, function(x){
        x[["4"]]
    })

    zero_wr_unique <- unlist(lapply(zero_wr_matches, function(x){
        !any(duplicated(x[["treat_index"]]))
    }))
    zero_wr_dist <- unlist(lapply(zero_wr_matches, function(x){
        sum(x[["distance"]])
    }))
    zero_wr_length <- unlist(lapply(zero_wr_matches, function(x){
        length(x[["treat_index"]])
    }))

    expect_true(all(zero_wr_unique))
    expect_true(all(zero_wr_length == sum(treat_vec)))


    some_match <- bipartite_matches(
        dist_mat = some_dist_mat,
        treat_vec = treat_vec,
        match_method = "with_replacement",
        n_sinks = 0
    )[["0"]]

    expect_true(abs(mean(zero_wr_dist / sum(some_match[["distance"]])) - 1) < 0.2)

    random_dist <- mean(some_dist_mat[
        cbind(1L:nrow(some_dist_mat),
              fixed_sample(1L:ncol(some_dist_mat),
                           nrow(some_dist_mat)))
    ])
    expect_true((mean(some_match[["distance"]]) / random_dist) < 0.5)

    all_optimal_matches <- all_bipartite_matches(x_mat = x_mat,
                                                 cov_x = cov(x_mat),
                                                 weight_list = weight_vecs,
                                                 treat_vec = treat_vec,
                                                 match_method = "optimal",
                                                 n_sinks = c(0L, 4L))

    zero_optimal_matches <- lapply(all_optimal_matches, function(x){
        x[["0"]]
    })
    four_optimal_matches <- lapply(all_optimal_matches, function(x){
        x[["0"]]
    })

    zero_opt_unique <- unlist(lapply(zero_optimal_matches, function(x){
        !any(duplicated(x[["treat_index"]]))
    }))
    zero_opt_dist <- unlist(lapply(zero_optimal_matches, function(x){
        sum(x[["distance"]])
    }))
    zero_opt_length <- unlist(lapply(zero_optimal_matches, function(x){
        length(x[["treat_index"]])
    }))

    expect_true(all(zero_opt_unique))
    expect_true(all(zero_opt_length == sum(treat_vec)))

    some_match <- bipartite_matches(
        dist_mat = some_dist_mat,
        treat_vec = treat_vec,
        match_method = "optimal",
        n_sinks = 0
    )[["0"]]

    expect_true(abs(mean(zero_opt_dist / sum(some_match[["distance"]])) - 1) < 0.3)

    random_dist <- mean(some_dist_mat[
        cbind(1L:nrow(some_dist_mat),
              fixed_sample(1L:ncol(some_dist_mat),
                           nrow(some_dist_mat)))
    ])
    expect_true((mean(some_match[["distance"]]) / random_dist) < 0.5)

    ##----------------

    expect_true(all(zero_opt_dist > zero_wr_dist))

    ##------------------------------------

    zero_wr_briers <- unlist(lapply(zero_wr_matches, function(ml) {
        brier_score_cv(x_mat,
                       match_list = ml)
    }))
    zero_opt_briers <- unlist(lapply(zero_optimal_matches, function(ml) {
        brier_score_cv(x_mat,
                       match_list = ml)
    }))
    four_wr_briers <- unlist(lapply(four_wr_matches, function(ml) {
        brier_score_cv(x_mat,
                       match_list = ml)
    }))
    four_opt_briers <- unlist(lapply(four_optimal_matches, function(ml) {
        brier_score_cv(x_mat,
                       match_list = ml)
    }))

    expect_true(mean(zero_wr_briers) > 0.1 && mean(zero_wr_briers) < 0.5)
    expect_true(mean(four_wr_briers) > 0.1 && mean(four_wr_briers) < 0.5)
    expect_true(mean(zero_opt_briers) > 0.1 && mean(zero_opt_briers) < 0.5)
    expect_true(mean(four_opt_briers) > 0.1 && mean(four_opt_briers) < 0.5)

    zero_wr_sd <- unlist(lapply(1L:num_weight_vecs, function(j) {
        bipartite_match_sd(x_mat,
                           cov_x = cov_x,
                           y_vector = y_vector,
                           match_list = zero_wr_matches[[j]],
                           treat_vec = treat_vec,
                           weight_vec = weight_vecs[[j]])
    }))
    zero_opt_sd <- unlist(lapply(1L:num_weight_vecs, function(j) {
        bipartite_match_sd(x_mat,
                           cov_x = cov_x,
                           y_vector = y_vector,
                           match_list = zero_optimal_matches[[j]],
                           treat_vec = treat_vec,
                           weight_vec = weight_vecs[[j]])
    }))
    four_wr_sd <- unlist(lapply(1L:num_weight_vecs, function(j) {
        bipartite_match_sd(x_mat,
                           cov_x = cov_x,
                           y_vector = y_vector,
                           match_list = four_wr_matches[[j]],
                           treat_vec = treat_vec,
                           weight_vec = weight_vecs[[j]])
    }))
    four_opt_sd <- unlist(lapply(1L:num_weight_vecs, function(j) {
        bipartite_match_sd(x_mat,
                           cov_x = cov_x,
                           y_vector = y_vector,
                           match_list = four_optimal_matches[[j]],
                           treat_vec = treat_vec,
                           weight_vec = weight_vecs[[j]])
    }))

    expect_true(mean(zero_wr_sd) / sd(y_vector) < 0.7)
    expect_true(mean(zero_opt_sd) / sd(y_vector) < 0.7)
    expect_true(mean(four_wr_sd) / sd(y_vector) < 0.7)
    expect_true(mean(four_opt_sd) / sd(y_vector) < 0.7)
})


test_that("testing sink_brier_bipartite_matches", {
    rows <- 100L
    num_weight_vecs <- 5L
    x_mat <- cbind(rnorm(rows),
                   runif(rows))
    treat_vec <- (1L:rows) %in% fixed_sample(1L:rows,
                                             45L)
    y_vector <- x_mat[, 1] + x_mat[, 2] + treat_vec * 0.3
    cov_x <- covariance_with_ranks(x_mat)

    some_dist_mat = weighted_mahal(x_mat,
                                   cov_x = cov_x,
                                   weight_vec = c(0.66, 0.33),
                                   treat_vec = treat_vec)

    weight_vecs <- generate_random_weights(prior_weights = c(2, 1),
                                           number_vectors = num_weight_vecs,
                                           minimum_weights = c(0.1, 0.1))

    all_wr_matches <- all_bipartite_matches(x_mat = x_mat,
                                            cov_x = cov_x,
                                            weight_list = weight_vecs,
                                            treat_vec = treat_vec,
                                            match_method = "with_replacement",
                                            n_sinks = c(0L, 4L))

    sink_brier_wr_matches <- sink_brier_bipartite_matches(
        x_mat = x_mat,
        cov_x = cov_x,
        weight_list = weight_vecs,
        treat_vec = treat_vec,
        match_method = "with_replacement",
        n_sinks = c(0L, 4L),
        silent = TRUE)

    expect_true(all(lengths(sink_brier_wr_matches) == 2L))
    expect_true(all(lengths(
        sink_brier_wr_matches[["matches_by_sinks"]]) == 5L))
    expect_true(all(lengths(
        sink_brier_wr_matches[["briers_by_sinks"]]) == 5L))

    expect_equal(all_wr_matches[[3]][["0"]],
                 sink_brier_wr_matches[["matches_by_sinks"]][["0"]][[3]])
    expect_equal(all_wr_matches[[5]][["1"]],
                 sink_brier_wr_matches[["matches_by_sinks"]][["1"]][[5]])

    brier_scores <- unlist(sink_brier_wr_matches[["briers_by_sinks"]])
    expect_true(all(0 < brier_scores & brier_scores < 1))
    expect_true(all(abs(brier_scores - 0.25) < 0.15))
})


test_that("testing permutation_bipartite_matches", {
    rows <- 40L
    num_weight_vecs <- 3L
    x_mat <- cbind(rnorm(rows),
                   runif(rows))
    treat_vec <- (1L:rows) %in% fixed_sample(1L:rows,
                                             15L)
    y_vector <- x_mat[, 1] + x_mat[, 2] + treat_vec * 0.3
    cov_x <- covariance_with_ranks(x_mat)

    some_dist_mat = weighted_mahal(x_mat,
                                   cov_x = cov_x,
                                   weight_vec = c(0.66, 0.33),
                                   treat_vec = treat_vec)

    weight_vecs <- generate_random_weights(prior_weights = c(2, 1),
                                           number_vectors = num_weight_vecs,
                                           minimum_weights = c(0.1, 0.1))

    sink_brier_wr_matches <- sink_brier_bipartite_matches(
        x_mat = x_mat,
        cov_x = cov_x,
        weight_list = weight_vecs,
        treat_vec = treat_vec,
        match_method = "with_replacement",
        n_sinks = c(0L, 4L),
        silent = TRUE)

    permutation_result <- permutation_bipartite_matches(
        matches_by_sinks = sink_brier_wr_matches[["matches_by_sinks"]],
        briers_by_sinks = sink_brier_wr_matches[["briers_by_sinks"]],
        n_sinks = c(0L, 4L))

    perm_briers <- permutation_result[["permutation_brier_scores"]]
    best_matches <- permutation_result[["best_matches"]]

    expect_equal(lengths(perm_briers), c(3L, 3L))
    expect_equal(lengths(best_matches), c(4L, 4L))
    expect_equal(names(best_matches[[1]]),
                 c("n_sinks", "raw_brier", "permutation_brier", "match_list"))

    expect_true(all(0 <= unlist(perm_briers) & unlist(perm_briers) <= 1))
})
