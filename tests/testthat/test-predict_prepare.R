context("testing predict_prepare")


test_that("testing predict_prepare", {
    x_mat <- cbind(
        1L:10L,
        runif(10L)
    )

    index_list <- list(
        treat_train = c(1L, 3L, 5L),
        control_train = c(2L, 4L, 6L),
        treat_test = c(7L, 9L),
        control_test = c(8L, 10L)
    )

    ##------------------------------------
    ## cross_all

    prep_result <- predict_prepare(x_mat,
                                   index_list,
                                   design = "cross_all")

    expect_equal(prep_result[["y_train"]],
                 c(1L, 1L, 1L, 0L, 0L, 0L))
    expect_equal(prep_result[["y_test"]],
                 c(1L, 1L))

    expect_equal(dim(prep_result[["x_train"]]),
                 c(6L, 4L))

    first_pair <- cbind(x_mat[index_list[["treat_train"]][1], , drop = FALSE],
                        x_mat[index_list[["control_train"]][1], , drop = FALSE])

    expect_equal(prep_result[["x_train"]][1, , drop = FALSE],
                 first_pair)

    first_pair_rev <- first_pair[, c(3, 4, 1, 2), drop = FALSE]
    expect_equal(prep_result[["x_train"]][1 + 3, , drop = FALSE],
                 first_pair_rev)

    ##------------------------------------
    ## cross_random

    prep_result <- predict_prepare(x_mat,
                                   index_list,
                                   design = "cross_random")

    expect_true(abs(sum(prep_result[["y_train"]]) -
                    0.5 * length(prep_result[["y_train"]])) < 1)
    expect_true(abs(sum(prep_result[["y_test"]]) -
                    0.5 * length(prep_result[["y_test"]])) < 1)

    expect_equal(dim(prep_result[["x_train"]]),
                 c(3L, 4L))

    first_pair <- cbind(x_mat[index_list[["treat_train"]][1], , drop = FALSE],
                        x_mat[index_list[["control_train"]][1], , drop = FALSE])
    first_pair_rev <- first_pair[, c(3, 4, 1, 2), drop = FALSE]
    if (prep_result[["y_train"]][1] == 1) {
        rel_pair <- first_pair
    } else {
        rel_pair <- first_pair_rev
    }

    expect_equal(prep_result[["x_train"]][1, , drop = FALSE],
                 rel_pair)

    ##------------------------------------
    ## differences_random

    prep_result <- predict_prepare(x_mat,
                                   index_list,
                                   design = "differences_random")

    expect_true(abs(sum(prep_result[["y_train"]]) -
                    0.5 * length(prep_result[["y_train"]])) < 1)
    expect_true(abs(sum(prep_result[["y_test"]]) -
                    0.5 * length(prep_result[["y_test"]])) < 1)

    expect_equal(dim(prep_result[["x_train"]]),
                 c(3L, 2L))

    first_pair <- x_mat[index_list[["treat_train"]][1], , drop = FALSE] -
        x_mat[index_list[["control_train"]][1], , drop = FALSE]
    first_pair_rev <- -first_pair
    if (prep_result[["y_train"]][1] == 1) {
        rel_pair <- first_pair
    } else {
        rel_pair <- first_pair_rev
    }

    expect_equal(prep_result[["x_train"]][1, , drop = FALSE],
                 rel_pair)

    ##------------------------------------
    ## differences_plain

    prep_result <- predict_prepare(x_mat,
                                   index_list,
                                   design = "differences_plain")

    expect_equal(prep_result[["y_train"]], c(1L, 1L, 1L))
    expect_equal(prep_result[["y_test"]], c(1L, 1L))

    expect_equal(dim(prep_result[["x_train"]]),
                 c(3L, 2L))

    first_pair <- x_mat[index_list[["treat_train"]][1], , drop = FALSE] -
        x_mat[index_list[["control_train"]][1], , drop = FALSE]
    expect_equal(prep_result[["x_train"]][1, , drop = FALSE],
                 first_pair)
})


test_that("testing index_list_from_match", {
    match_list <- list(
        treat_index = c(1L, 2L, 3L, 10L, 4L),
        treat_index_within = c(1L, 2L, 3L, 5L, 4L),
        control_index = c(5L, 6L, 8L, 9L, 7L),
        control_index_within = c(1L, 2L, 4L, 5L, 3L),
        distance = runif(5)
    )
    train_index <- c(TRUE, TRUE, FALSE, TRUE, FALSE)

    res_list <- index_list_from_match(match_list,
                                      train_index)

    expect_equal(res_list[["treat_train"]],
                 match_list[["treat_index"]][train_index])
    expect_equal(res_list[["control_train"]],
                 match_list[["control_index"]][train_index])
    expect_equal(res_list[["treat_test"]],
                 match_list[["treat_index"]][!train_index])
    expect_equal(res_list[["control_test"]],
                 match_list[["control_index"]][!train_index])
})


test_that("testing generate_train_test_split", {
    match_list <- list(
        treat_index = c(1L, 2L, 3L, 10L, 4L),
        treat_index_within = c(1L, 2L, 3L, 5L, 4L),
        control_index = c(5L, 6L, 8L, 9L, 7L),
        control_index_within = c(1L, 2L, 4L, 5L, 3L),
        distance = runif(5)
    )

    res_split <- generate_train_test_split(match_list, 0.6)

    expect_equal(unname(sort(unlist(res_split))),
                 1L:10L)
    expect_equal(unname(lengths(res_split)),
                 c(3L, 3L, 2L, 2L))

    expect_equal(sort(c(res_split[["treat_train"]],
                        res_split[["treat_test"]])),
                 sort(match_list[["treat_index"]]))
})


test_that("testing generate_k_fold_index", {
    match_list <- list(
        treat_index = 1L:20L,
        treat_index_within = 1L:20L,
        control_index = 21L:40L,
        control_index_within = 1L:20L,
        distance = runif(20)
    )

    res_list_list <- generate_k_fold_index(match_list, 5L)

    res_treat_train <- lapply(res_list_list, function(x) {
        x[["treat_train"]]
    })
    res_control_train <- lapply(res_list_list, function(x) {
        x[["control_train"]]
    })
    res_treat_test <- lapply(res_list_list, function(x) {
        x[["treat_test"]]
    })
    res_control_test <- lapply(res_list_list, function(x) {
        x[["control_test"]]
    })

    expect_equal(lengths(res_treat_train), rep(16L, 5L))
    expect_equal(lengths(res_control_train), rep(16L, 5L))
    expect_equal(lengths(res_treat_test), rep(4L, 5L))
    expect_equal(lengths(res_control_test), rep(4L, 5L))

    expect_equal(sort(unlist(res_treat_test)),
                 1L:20L)
    expect_equal(sort(unlist(res_control_test)),
                 21L:40L)
})
