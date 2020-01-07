context("testing match_eval functions")


test_that("testing match_estimate", {
    match_list <- list(
        treat_index = 1L:200L,
        control_index = 201L:400L,
        distance = runif(200L)
    )

    y_vector <- rnorm(400L)

    ##------------------------------------

    expect_equal(match_estimate(match_list = match_list,
                                y_vector = y_vector),
                 mean(y_vector[1L:200L] - y_vector[201L:400L]))

    ##------------------------------------

    expect_error(match_estimate(match_list = match_list,
                                y_vector = y_vector,
                                treat_vec = c(0L, 0L, 1L)))

    expect_error(match_estimate(match_list = match_list,
                                y_vector = y_vector,
                                treat_vec = c(rep(0L, 200L),
                                              rep(1L, 200L))))

    expect_error(match_estimate(match_list = match_list,
                                y_vector = y_vector,
                                treat_vec = c(rep(1L, 200L),
                                              rep(0L, 200L))),
                 NA)
})


test_that("testing brier_score", {
    x_mat <- cbind(
        c(runif(200L, 0, 3),
          runif(200L, 1, 4)),
        rnorm(400L)
    )

    match_list <- list(
        treat_index = 1L:200L,
        treat_index_within = 1L:200L,
        control_index = 201L:400L,
        control_index_within = 1L:200L,
        distance = runif(200L)
    )

    index_list <- index_list_from_match(match_list,
                                        c(rep(TRUE, 120L),
                                          rep(FALSE, 80L)))

    train_test_list <- predict_prepare(x_mat,
                                       index_list,
                                       "cross_all")

    brier_res <- brier_score(train_test_list,
                             match_predict_xgb,
                             avg = TRUE)
    expect_true(brier_res < 0.25 && brier_res > 0.01)

    ##------------------------------------

    ## now just random
    x_mat <- cbind(
        runif(400L),
        rnorm(400L)
    )

    match_list <- list(
        treat_index = 1L:200L,
        treat_index_within = 1L:200L,
        control_index = 201L:400L,
        control_index_within = 1L:200L,
        distance = runif(200L)
    )

    index_list <- index_list_from_match(match_list,
                                        c(rep(TRUE, 120L),
                                          rep(FALSE, 80L)))

    train_test_list <- predict_prepare(x_mat,
                                       index_list,
                                       "cross_all")

    brier_res <- brier_score(train_test_list,
                             match_predict_xgb,
                             avg = FALSE)

    ## on average will be larger than 0.25, and we have 80
    expect_true(abs(brier_res - 0.25 * 80L) < 0.15 * 80L)
})


test_that("testing brier_score_split", {
    x_mat <- cbind(
        runif(400L),
        rnorm(400L)
    )

    match_list <- list(
        treat_index = 1L:200L,
        treat_index_within = 1L:200L,
        control_index = 201L:400L,
        control_index_within = 1L:200L,
        distance = runif(200L)
    )

    brier_res <- brier_score_split(x_mat,
                                   match_list,
                                   design = "cross_all",
                                   train_fraction = 0.7,
                                   match_predict_xgb)
    expect_true(abs(brier_res - 0.25) < 0.1)
})


test_that("testing brier_score_cv", {
    x_mat <- cbind(
        runif(400L),
        rnorm(400L)
    )

    match_list <- list(
        treat_index = 1L:200L,
        treat_index_within = 1L:200L,
        control_index = 201L:400L,
        control_index_within = 1L:200L,
        distance = runif(200L)
    )

    brier_res <- brier_score_cv(x_mat,
                                match_list,
                                design = "cross_all",
                                num_folds = 5,
                                match_predict_xgb)

    ## little tighter than with split
    expect_true(abs(brier_res - 0.28) < 0.06)
})


test_that("testing permutation_brier", {
    x_mat <- cbind(
        c(runif(200L, 0, 2),
          runif(200L, 1, 3)),
        rnorm(400L)
    )

    match_list <- list(
        treat_index = 1L:200L,
        treat_index_within = 1L:200L,
        control_index = 201L:400L,
        control_index_within = 1L:200L,
        distance = runif(200L)
    )

    expect_error(permutation_brier(x_mat,
                                   match_list,
                                   design = "cross_all",
                                   use_cv = TRUE,
                                   num_permutations = 20L,
                                   match_predict_function = match_predict_xgb,
                                   train_fraction = 0.7))
    expect_error(permutation_brier(x_mat,
                                   match_list,
                                   design = "cross_all",
                                   use_cv = FALSE,
                                   num_permutations = 20L,
                                   match_predict_function = match_predict_xgb,
                                   num_folds = 5))
    ## see below for non-error versions

    ##------------------------------------
    ## cv

    ## easy to predict:
    brier_res <- brier_score_cv(x_mat,
                                match_list,
                                design = "cross_all",
                                num_folds = 5,
                                match_predict_xgb)
    expect_true(brier_res < 0.15)

    ## now get the dist:
    perm_res <- permutation_brier(x_mat,
                                  match_list,
                                  design = "cross_all",
                                  use_cv = TRUE,
                                  num_permutations = 8L,
                                  match_predict_function = match_predict_xgb,
                                  num_folds = 3L)
    expect_true(abs(mean(perm_res) - 0.28) < 0.04)

    ##------------------------------------
    ## split

    ## still easy to predict:
    brier_res <- brier_score_split(x_mat,
                                   match_list,
                                   design = "cross_all",
                                   train_fraction = 0.7,
                                   match_predict_xgb)
    expect_true(brier_res < 0.15)

    ## now get the dist:
    perm_res <- permutation_brier(x_mat,
                                  match_list,
                                  design = "cross_all",
                                  use_cv = FALSE,
                                  num_permutations = 8L,
                                  match_predict_function = match_predict_xgb,
                                  train_fraction = 0.7)

    expect_true(abs(mean(perm_res) - 0.28) < 0.05)
})
