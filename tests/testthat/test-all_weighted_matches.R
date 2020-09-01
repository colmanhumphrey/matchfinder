context("testing all weighted matches")


test_that("testing all_bipartite_matches", {
    ##------------------------------------
    ## no caliper/propens

    rows <- 200L
    num_weight_vecs <- 5L
    x_mat <- cbind(rnorm(rows),
                   runif(rows))
    treat_vec <- (1L:rows) %in% fixed_sample(1L:rows,
                                             45L)
    y_vector <- x_mat[, 1] + x_mat[, 2] + treat_vec * 0.3 + rnorm(rows)
    cov_x <- covariance_with_ranks(x_mat)

    some_dist_mat <- weighted_mahal(x_mat,
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

    zero_wr_matches <- all_wr_matches[["0"]]
    four_wr_matches <- all_wr_matches[["4"]]

    zero_wr_unique <- unlist(lapply(zero_wr_matches, function(x) {
        !any(duplicated(x[["treat_index"]]))
    }))
    zero_wr_dist <- unlist(lapply(zero_wr_matches, function(x) {
        sum(x[["distance"]])
    }))
    zero_wr_length <- unlist(lapply(zero_wr_matches, function(x) {
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

    expect_true(abs(mean(zero_wr_dist /
                         sum(some_match[["distance"]])) - 1) < 0.2)

    random_dist <- mean(some_dist_mat[
        cbind(seq_len(nrow(some_dist_mat)),
              fixed_sample(seq_len(ncol(some_dist_mat)),
                           nrow(some_dist_mat)))
    ])
    expect_true((mean(some_match[["distance"]]) / random_dist) < 0.5)

    all_optimal_matches <- all_bipartite_matches(x_mat = x_mat,
                                                 cov_x = cov(x_mat),
                                                 weight_list = weight_vecs,
                                                 treat_vec = treat_vec,
                                                 match_method = "optimal",
                                                 n_sinks = c(0L, 4L))

    zero_optimal_matches <- all_optimal_matches[["0"]]
    four_optimal_matches <- all_optimal_matches[["4"]]

    zero_opt_unique <- unlist(lapply(zero_optimal_matches, function(x) {
        !any(duplicated(x[["treat_index"]]))
    }))
    zero_opt_dist <- unlist(lapply(zero_optimal_matches, function(x) {
        sum(x[["distance"]])
    }))
    zero_opt_length <- unlist(lapply(zero_optimal_matches, function(x) {
        length(x[["treat_index"]])
    }))

    expect_true(all(zero_opt_unique))
    expect_true(all(zero_opt_length == sum(treat_vec)))

    some_match <- bipartite_matches(
        dist_mat = some_dist_mat,
        treat_vec = treat_vec,
        match_method = "optimal",
        n_sinks = 0L
    )[["0"]]

    expect_true(abs(mean(zero_opt_dist /
                         sum(some_match[["distance"]])) - 1) < 0.3)

    random_dist <- mean(some_dist_mat[
        cbind(seq_len(nrow(some_dist_mat)),
              fixed_sample(seq_len(ncol(some_dist_mat)),
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

    ## random differences give sqrt(2), so we're below that
    expected_random_sd <- sd(y_vector) * sqrt(2)

    ## opt does great here
    expect_true(mean(zero_opt_sd) < 0.85 * expected_random_sd)
    expect_true(mean(four_opt_sd) < 0.85 * expected_random_sd)

    ## wr not as good on this metric, as expected, since we repeat controls
    expect_true(mean(zero_wr_sd) < 0.95 * expected_random_sd)
    expect_true(mean(four_wr_sd) < 0.95 * expected_random_sd)
})


test_that("testing all_nonbipartite_matches", {
    num_weight_vecs <- 5L
    x_mat <- cbind(
        c(
            4.1, 49.9, 34.3, 22.1, 3.1, 7, 6, 1, 28, 1, 7, 10,
            27, 10, 25, 40, 2, 6, 47, 4
        ),
        c(
            31, 4, 3, 34, 1, 29, 2, 7, 5.1, 3.1, 54, 54,
            15, 6, 11, 15, 7, 16.7, 22.9, 20.3
        )
    )
    cov_x <- covariance_with_ranks(x_mat)

    dist_mat <- weighted_mahal(
        x_mat,
        cov_x
    )
    diag(dist_mat) <- Inf

    ## a random tolerance for now
    tolerance_vec <- runif(nrow(x_mat))

    weight_vecs <- generate_random_weights(
        prior_weights = c(2, 1),
        number_vectors = num_weight_vecs,
        minimum_weights = c(0.1, 0.1)
    )

    wr_res <- with_replacement_nbp_match(dist_mat,
        tolerance_vec = tolerance_vec,
        keep_all = FALSE
    )

    ## we'll nearly always want to explicitly give tolerance info,
    ## so we throw a warning
    expect_warning(
        all_nonbipartite_matches(
            x_mat = x_mat,
            cov_x = cov_x,
            weight_list = weight_vecs,
            match_method = "with_replacement",
            n_sinks = c(0L, 1L, 2L)
        )
    )

    ## It's the same result as passing explicit NULL,
    ## but explicit NULL is assumed intentional
    expect_warning(
        all_nonbipartite_matches(
            x_mat = x_mat,
            cov_x = cov_x,
            tolerance_list = NULL,
            weight_list = weight_vecs,
            match_method = "with_replacement",
            n_sinks = c(0L, 1L, 2L)
        ),
        NA
    )

    tol_list <- gen_tolerance_list(
            tolerance_vec = tolerance_vec
    )

    all_wr_matches <- all_nonbipartite_matches(
        x_mat = x_mat,
        cov_x = cov_x,
        weight_list = weight_vecs,
        tolerance_list = tol_list,
        match_method = "with_replacement",
        n_sinks = c(0L, 1L, 2L)
    )

    ## verify that the treated units are larger in tol
    treat_larger_in_tol <- unlist(lapply(
        all_wr_matches, function(matches_by_sink) {
            !any(unlist(lapply(matches_by_sink, function(match_list) {
                tolerance_check(match_list, tol_list)[["error"]]
            })))
        }
    ))

    expect_true(all(treat_larger_in_tol))

    ##------------------------------------
    ## no caliper/propens

    rows <- 400L
    num_weight_vecs <- 5L
    x_mat <- cbind(rnorm(rows),
                   runif(rows))
    tol_vec <- runif(rows) * 2 + x_mat[, 1] * 0.3
    tol_list <- gen_tolerance_list(
        tolerance_vec = tol_vec,
        tolerance_min = 0.1
    )
    y_vector <- x_mat[, 1] + x_mat[, 2] + tol_vec * 0.3 + rnorm(rows)
    cov_x <- covariance_with_ranks(x_mat)

    some_dist_mat <- weighted_mahal(x_mat,
        cov_x = cov_x,
        weight_vec = c(0.4, 0.6)
    )

    weight_vecs <- generate_random_weights(prior_weights = c(2, 3),
                                           number_vectors = num_weight_vecs,
                                           minimum_weights = c(0.1, 0.1))

    all_wr_matches <- all_nonbipartite_matches(
        x_mat = x_mat,
        cov_x = cov(x_mat),
        weight_list = weight_vecs,
        tolerance_list = tol_list,
        match_method = "with_replacement",
        n_sinks = c(0L, 4L)
    )

    zero_wr_matches <- all_wr_matches[["0"]]
    four_wr_matches <- all_wr_matches[["4"]]

    zero_wr_unique <- unlist(lapply(zero_wr_matches, function(x) {
        !any(duplicated(x[["treat_index"]]))
    }))
    zero_wr_dist <- unlist(lapply(zero_wr_matches, function(x) {
        sum(x[["distance"]])
    }))
    zero_wr_length <- unlist(lapply(zero_wr_matches, function(x) {
        length(x[["treat_index"]])
    }))

    expect_true(all(zero_wr_unique))
    ## should be able to use all
    expect_true(all(zero_wr_length == rows %/% 2L))

    some_match <- nonbipartite_matches(
        dist_mat = some_dist_mat,
        tolerance_list = tol_list,
        match_method = "with_replacement",
        n_sinks = 0L
    )[["0"]]

    expect_true(abs(mean(zero_wr_dist /
                         sum(some_match[["distance"]])) - 1) < 0.2)

    random_dist <- mean(some_dist_mat[
        cbind(seq_len(nrow(some_dist_mat)),
              fixed_sample(seq_len(ncol(some_dist_mat)),
                           nrow(some_dist_mat)))
    ])
    expect_true((mean(some_match[["distance"]]) / random_dist) < 0.5)

    all_optimal_matches <- all_nonbipartite_matches(
        x_mat = x_mat,
        cov_x = cov(x_mat),
        weight_list = weight_vecs,
        tolerance_list = tol_list,
        match_method = "optimal",
        n_sinks = c(0L, 4L)
    )

    zero_optimal_matches <- all_optimal_matches[["0"]]
    four_optimal_matches <- all_optimal_matches[["4"]]

    zero_opt_unique <- unlist(lapply(zero_optimal_matches, function(x) {
        !any(duplicated(x[["treat_index"]]))
    }))
    zero_opt_dist <- unlist(lapply(zero_optimal_matches, function(x) {
        sum(x[["distance"]])
    }))
    zero_opt_length <- unlist(lapply(zero_optimal_matches, function(x) {
        length(x[["treat_index"]])
    }))

    expect_true(all(zero_opt_unique))
    ## should still be able to use all
    expect_true(all(zero_opt_length == rows %/% 2L))

    some_match <- nonbipartite_matches(
        dist_mat = some_dist_mat,
        tolerance_list = tol_list,
        match_method = "optimal",
        n_sinks = 0
    )[["0"]]

    expect_true(abs(mean(zero_opt_dist /
                         sum(some_match[["distance"]])) - 1) < 0.3)

    random_dist <- mean(some_dist_mat[
        cbind(seq_len(nrow(some_dist_mat)),
              fixed_sample(seq_len(ncol(some_dist_mat)),
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
        nonbipartite_match_sd_scaled(
            x_mat,
            cov_x = cov_x,
            y_vector = y_vector,
            match_list = zero_wr_matches[[j]],
            tolerance_list = tol_list,
            weight_vec = weight_vecs[[j]]
        )
    }))
    zero_opt_sd <- unlist(lapply(1L:num_weight_vecs, function(j) {
        nonbipartite_match_sd_scaled(
            x_mat,
            cov_x = cov_x,
            y_vector = y_vector,
            match_list = zero_optimal_matches[[j]],
            tolerance_list = tol_list,
            weight_vec = weight_vecs[[j]]
        )
    }))
    four_wr_sd <- unlist(lapply(1L:num_weight_vecs, function(j) {
        nonbipartite_match_sd_scaled(
            x_mat,
            cov_x = cov_x,
            y_vector = y_vector,
            match_list = four_wr_matches[[j]],
            tolerance_list = tol_list,
            weight_vec = weight_vecs[[j]]
        )
    }))
    four_opt_sd <- unlist(lapply(1L:num_weight_vecs, function(j) {
        nonbipartite_match_sd_scaled(
            x_mat,
            cov_x = cov_x,
            y_vector = y_vector,
            match_list = four_optimal_matches[[j]],
            tolerance_list = tol_list,
            weight_vec = weight_vecs[[j]]
        )
    }))

    expected_random <- y_tolerance_diff_ratio(
        y_vector = y_vector,
        tolerance_list = tol_list
    )
    expected_random_sd <- expected_random[["sd"]]

    ## opt does great here
    expect_true(mean(zero_opt_sd) < 0.95 * expected_random_sd)
    expect_true(mean(four_opt_sd) < 0.95 * expected_random_sd)

    ## wr not as good on this metric, as expected, since we repeat controls
    ## might be worth looking into that more # TODO
})


test_that("testing brier_bipartite_matches", {
    rows <- 100L
    num_weight_vecs <- 5L
    x_mat <- cbind(rnorm(rows),
                   runif(rows))
    treat_vec <- (1L:rows) %in% fixed_sample(1L:rows,
                                             45L)
    y_vector <- x_mat[, 1] + x_mat[, 2] + treat_vec * 0.3
    cov_x <- covariance_with_ranks(x_mat)

    weight_vecs <- generate_random_weights(prior_weights = c(2, 1),
                                           number_vectors = num_weight_vecs,
                                           minimum_weights = c(0.1, 0.1))

    all_wr_matches <- all_bipartite_matches(x_mat = x_mat,
                                            cov_x = cov_x,
                                            weight_list = weight_vecs,
                                            treat_vec = treat_vec,
                                            match_method = "with_replacement",
                                            n_sinks = c(0L, 4L))

    brier_wr_matches <- brier_bipartite_matches(
        x_mat = x_mat,
        cov_x = cov_x,
        weight_list = weight_vecs,
        treat_vec = treat_vec,
        match_method = "with_replacement",
        n_sinks = c(0L, 4L),
        silent = TRUE)

    expect_true(all(lengths(brier_wr_matches) == 2L))
    expect_true(all(lengths(
        brier_wr_matches[["matches_by_sinks"]]) == 5L))
    expect_true(all(lengths(
        brier_wr_matches[["briers_by_sinks"]]) == 5L))

    expect_equal(all_wr_matches[["0"]][[3]],
                 brier_wr_matches[["matches_by_sinks"]][["0"]][[3]])
    expect_equal(all_wr_matches[["1"]][[5]],
                 brier_wr_matches[["matches_by_sinks"]][["1"]][[5]])

    brier_scores <- unlist(brier_wr_matches[["briers_by_sinks"]])
    expect_true(all(0 < brier_scores & brier_scores < 1))
    expect_true(all(abs(brier_scores - 0.25) < 0.15))
})


test_that("testing brier_nonbipartite_matches", {
    treat_effect <- 0.3
    rows <- 200L
    num_weight_vecs <- 5L
    x_mat <- cbind(rnorm(rows),
                   runif(rows))
    tol_vec <- runif(rows) + x_mat[, 1L]
    y_vector <- x_mat[, 1L] + x_mat[, 2L] + tol_vec * treat_effect + rnorm(rows)
    cov_x <- covariance_with_ranks(x_mat)

    tol_list <- gen_tolerance_list(
        tolerance_vec = tol_vec
    )

    weight_vecs <- generate_random_weights(prior_weights = c(2, 1),
                                           number_vectors = num_weight_vecs,
                                           minimum_weights = c(0.1, 0.1))

    all_wr_matches <- all_nonbipartite_matches(
        x_mat = x_mat,
        cov_x = cov_x,
        weight_list = weight_vecs,
        tolerance_list = tol_list,
        match_method = "with_replacement",
        n_sinks = c(0L, 4L)
    )

    brier_wr_matches <- brier_nonbipartite_matches(
        x_mat = x_mat,
        cov_x = cov_x,
        weight_list = weight_vecs,
        tolerance_list = tol_list,
        match_method = "with_replacement",
        n_sinks = c(0L, 4L),
        silent = TRUE
    )

    expect_true(all(lengths(brier_wr_matches) == 2L))
    expect_true(all(lengths(
        brier_wr_matches[["matches_by_sinks"]]) == 5L))
    expect_true(all(lengths(
        brier_wr_matches[["briers_by_sinks"]]) == 5L))

    expect_equal(all_wr_matches[["0"]][[3]],
                 brier_wr_matches[["matches_by_sinks"]][["0"]][[3]])
    expect_equal(all_wr_matches[["1"]][[5]],
                 brier_wr_matches[["matches_by_sinks"]][["1"]][[5]])

    brier_scores <- unlist(brier_wr_matches[["briers_by_sinks"]])
    expect_true(all(0 < brier_scores & brier_scores < 1))
    expect_true(all(abs(brier_scores - 0.25) < 0.15))
})


test_that("testing permutation_matches", {
    ## ------------------------------------
    ## bipartite
    rows <- 40L
    num_weight_vecs <- 3L
    n_sinks_vec <- c(0L, 4L)
    treat_effect <- 0.3

    x_mat <- cbind(rnorm(rows),
                   runif(rows))
    treat_vec <- (1L:rows) %in% fixed_sample(1L:rows,
                                             15L)
    ## don't need it here
    ## y_vector <- tol_vec * treat_effect + <x stuff> + <noise>
    cov_x <- covariance_with_ranks(x_mat)

    weight_vecs <- generate_random_weights(prior_weights = c(2, 1),
                                           number_vectors = num_weight_vecs,
                                           minimum_weights = c(0.1, 0.1))

    brier_wr_matches <- brier_bipartite_matches(
        x_mat = x_mat,
        cov_x = cov_x,
        weight_list = weight_vecs,
        treat_vec = treat_vec,
        match_method = "with_replacement",
        n_sinks = n_sinks_vec,
        silent = TRUE)

    permutation_result <- permutation_matches(
        matches_by_sinks = brier_wr_matches[["matches_by_sinks"]],
        briers_by_sinks = brier_wr_matches[["briers_by_sinks"]],
        x_mat = x_mat,
        n_sinks = n_sinks_vec)

    perm_briers <- permutation_result[["permutation_brier_scores"]]
    best_matches <- permutation_result[["best_matches"]]

    expect_equal(lengths(perm_briers),
                 rep(num_weight_vecs, times = length(n_sinks_vec)))
    expect_equal(length(best_matches), length(n_sinks_vec))
    expect_equal(names(best_matches[[1]]),
                 c("n_sinks", "raw_brier", "permutation_brier", "match_list"))

    expect_true(all(0 <= unlist(perm_briers) & unlist(perm_briers) <= 1))

    ## ----------------
    ## and if we don't want to approximate (very slow!):

    long_permutation_result <- permutation_matches(
        matches_by_sinks = brier_wr_matches[["matches_by_sinks"]],
        briers_by_sinks = brier_wr_matches[["briers_by_sinks"]],
        x_mat = x_mat,
        n_sinks = n_sinks_vec,
        approximate_by_best = FALSE
    )
    long_perm_briers <- long_permutation_result[["permutation_brier_scores"]]

    ## overall should be relatively similar
    approx_scores <- unlist(perm_briers)
    full_scores <- unlist(long_perm_briers)

    abs_perc_diff <- 2 * abs(approx_scores - full_scores) /
        abs(approx_scores + full_scores)

    expect_true(mean(abs_perc_diff) < 0.2)

    ## ------------------------------------
    ## nonbipartite

    treat_effect <- 0.3
    rows <- 200L
    num_weight_vecs <- 5L
    n_sinks_vec <- c(0L, 4L)

    x_mat <- cbind(rnorm(rows),
                   runif(rows))
    tol_vec <- runif(rows) + x_mat[, 1L]
    ## don't need it here
    ## y_vector <- tol_vec * treat_effect + <x stuff> + <noise>
    cov_x <- covariance_with_ranks(x_mat)

    tol_list <- gen_tolerance_list(
        tolerance_vec = tol_vec
    )

    weight_vecs <- generate_random_weights(prior_weights = c(2, 1),
                                           number_vectors = num_weight_vecs,
                                           minimum_weights = c(0.1, 0.1))

    brier_wr_matches <- brier_nonbipartite_matches(
        x_mat = x_mat,
        cov_x = cov_x,
        weight_list = weight_vecs,
        tolerance_list = tol_list,
        match_method = "with_replacement",
        n_sinks = n_sinks_vec,
        silent = TRUE
    )

    permutation_result <- permutation_matches(
        matches_by_sinks = brier_wr_matches[["matches_by_sinks"]],
        briers_by_sinks = brier_wr_matches[["briers_by_sinks"]],
        x_mat = x_mat,
        n_sinks = n_sinks_vec)

    perm_briers <- permutation_result[["permutation_brier_scores"]]
    best_matches <- permutation_result[["best_matches"]]

    expect_equal(lengths(perm_briers),
                 rep(num_weight_vecs, times = length(n_sinks_vec)))
    expect_equal(length(best_matches), length(n_sinks_vec))
    expect_equal(names(best_matches[[1]]),
                 c("n_sinks", "raw_brier", "permutation_brier", "match_list"))

    expect_true(all(0 <= unlist(perm_briers) & unlist(perm_briers) <= 1))
})
