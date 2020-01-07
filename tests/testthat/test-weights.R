context("testing random weight vec generation")


test_that("testing generate_random_weights", {
    ##------------------------------------
    ## first with no hierarchical info

    prior_weights <- 1L:10L

    res_vecs <- generate_random_weights(prior_weights,
                                        number_vectors = 1000L)
    vec_means <- Reduce('+', res_vecs) / length(res_vecs)

    expected_diff <- vec_means - (prior_weights / sum(prior_weights))

    expect_true(sum(abs(expected_diff)) < 0.1)

    ##----------------

    ## add mins

    res_vecs <- generate_random_weights(prior_weights,
                                        number_vectors = 1000L,
                                        minimum_weights = rep(1 / 100, 10L))
    vec_means <- Reduce('+', res_vecs) / length(res_vecs)

    ## we expected this to be balanced also 1:10, but only
    ## beyond the mins
    vec_beyond_mins <- vec_means - 1 / 100
    vec_beyond_mins <- vec_beyond_mins / sum(vec_beyond_mins)

    expected_diff <- vec_beyond_mins - (prior_weights / sum(prior_weights))

    expect_true(sum(abs(expected_diff)) < 0.1)

    ##------------------------------------
    ## hierarchical
    ## not sure why you'd do this full mess

    hier_list <- list(
        list(index = c(1L, 10L), weight = 3),
        list(index = 2L:7L, weight = 2, variance = 10),
        list(index = 8L:9L, weight = 1, variance = 0.001)
    )

    prior_weights <- runif(10)
    prior_weights[c(1, 10)] <- c(1, 5)
    minimum_weights <- runif(10, 0, 0.03)

    res_vecs <- generate_random_weights(prior_weights = prior_weights,
                                        number_vectors = 1000L,
                                        minimum_weights = minimum_weights,
                                        hierarchical_list = hier_list)

    vec_means <- Reduce('+', res_vecs) / length(res_vecs)
    group_means <- lapply(hier_list, function(x){
        vec_means[x[["index"]]]
    })

    expect_true(abs(sum(group_means[[1]]) - 3 / 6) < 0.15)
    expect_true(abs(sum(group_means[[2]]) - 2 / 6) < 0.15)
    expect_true(abs(sum(group_means[[3]]) - 1 / 6) < 0.15)

    group_one_ratio <- unlist(lapply(res_vecs, function(x){
        group_one <- x[hier_list[[1]][["index"]]]
        group_one_minus_min <- group_one -
            minimum_weights[hier_list[[1]][["index"]]]
        group_one_minus_min[2] / group_one_minus_min[1]
    }))

    expect_true(abs(median(group_one_ratio) -
                    prior_weights[10] / prior_weights[1]) < 0.5)
})


test_that("testing hierarchical_random_weights", {
    hier_list <- list(
        list(index = c(1L, 10L), weight = 3),
        list(index = 2L:7L, weight = 2),
        list(index = 8L:9L, weight = 1)
    )

    hier_list_adj <- lapply(hier_list, function(x){
        if (is.null(x[["variance"]])) {
            ## scale invariance
            x[["gamma_alpha"]] <- 1
            x[["gamma_beta"]] <- 1 / x[["weight"]]
        } else {
            x[["gamma_beta"]] <- x[["weight"]] / x[["variance"]]
            x[["gamma_alpha"]] <- x[["weight"]] * x[["gamma_beta"]]
        }
        return(x)
    })

    res_vecs <- hierarchical_random_weights(prior_weights = rep(1, 10L),
                                            number_vectors = 1000L,
                                            minimum_weights = rep(0, 10L),
                                            hierarchical_list = hier_list_adj)
    vec_means <- Reduce('+', res_vecs) / length(res_vecs)

    group_means <- lapply(hier_list, function(x){
        vec_means[x[["index"]]]
    })

    ## they won't be exact, it's a gamma dist, but close
    expect_true(abs(sum(group_means[[1]]) - 3 / 6) < 0.1)
    expect_true(abs(sum(group_means[[2]]) - 2 / 6) < 0.1)
    expect_true(abs(sum(group_means[[3]]) - 1 / 6) < 0.1)

    ##------------------------------------
    ## with low variance

    hier_list <- list(
        list(index = c(1L, 10L), weight = 3, variance = 0.0001),
        list(index = 2L:7L, weight = 2, variance = 0.0001),
        list(index = 8L:9L, weight = 1, variance = 0.0001)
    )

    hier_list_adj <- lapply(hier_list, function(x){
        if (is.null(x[["variance"]])) {
            ## scale invariance
            x[["gamma_alpha"]] <- 1
            x[["gamma_beta"]] <- 1 / x[["weight"]]
        } else {
            x[["gamma_beta"]] <- x[["weight"]] / x[["variance"]]
            x[["gamma_alpha"]] <- x[["weight"]] * x[["gamma_beta"]]
        }
        return(x)
    })

    res_vecs <- hierarchical_random_weights(prior_weights = rep(1, 10L),
                                            number_vectors = 10L,
                                            minimum_weights = rep(0, 10L),
                                            hierarchical_list = hier_list_adj)
    vec_means <- Reduce('+', res_vecs) / length(res_vecs)

    group_means <- lapply(hier_list, function(x){
        vec_means[x[["index"]]]
    })

    ## should now be very close
    expect_true(abs(sum(group_means[[1]]) - 3 / 6) < 0.01)
    expect_true(abs(sum(group_means[[2]]) - 2 / 6) < 0.01)
    expect_true(abs(sum(group_means[[3]]) - 1 / 6) < 0.01)

    ## but not necessarily within group
    max_diff <- unlist(lapply(group_means, function(x){
        max(x) - min(x)
    }))
    expect_true(max(max_diff) > 0.01)

    ##------------------------------------
    ## add prior weights

    hier_list <- list(
        list(index = c(1L, 10L), weight = 3),
        list(index = 2L:7L, weight = 2),
        list(index = 8L:9L, weight = 1)
    )

    prior_weights <- 1:10

    hier_list_adj <- lapply(hier_list, function(x){
        if (is.null(x[["variance"]])) {
            ## scale invariance
            x[["gamma_alpha"]] <- 1
            x[["gamma_beta"]] <- 1 / x[["weight"]]
        } else {
            x[["gamma_beta"]] <- x[["weight"]] / x[["variance"]]
            x[["gamma_alpha"]] <- x[["weight"]] * x[["gamma_beta"]]
        }
        return(x)
    })

    res_vecs <- hierarchical_random_weights(prior_weights = prior_weights,
                                            number_vectors = 10000L,
                                            minimum_weights = rep(0, 10L),
                                            hierarchical_list = hier_list_adj)
    vec_means <- Reduce('+', res_vecs) / length(res_vecs)

    group_means <- lapply(hier_list, function(x){
        vec_means[x[["index"]]]
    })

    ## still should get group means about right
    expect_true(abs(sum(group_means[[1]]) - 3 / 6) < 0.1)
    expect_true(abs(sum(group_means[[2]]) - 2 / 6) < 0.1)
    expect_true(abs(sum(group_means[[3]]) - 1 / 6) < 0.1)

    ## and within group, we should apply the weights
    group_2_rescale <- group_means[[2]] / sum(group_means[[2]])
    diff_to_expected <- group_2_rescale - (2:7) / sum(2:7)
    expect_true(sum(abs(diff_to_expected)) < 0.1)
})
