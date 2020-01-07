context("testing propensity")


test_that("testing gen_propensity_list", {
    expect_error(gen_propensity_list(oos_propensity = 5L))
    expect_error(gen_propensity_list(oos_propensity = c(TRUE, FALSE)))
    expect_error(gen_propensity_list(oos_propensity = FALSE,
                                     n_folds = 5L))
    expect_error(gen_propensity_list(oos_propensity = TRUE,
                                     n_folds = 5L),
                 NA)

    expect_equal(names(gen_propensity_list()),
                 c("propensity_function", "oos_propensity", "n_folds"))

    expect_true(gen_propensity_list(n_folds = 10L)[["oos_propensity"]])
    expect_false(gen_propensity_list(n_folds = NULL)[["oos_propensity"]])
})


test_that("testing propensity_score", {
    n_rows <- 1000L
    n_cols <- 10L

    x_mat <- matrix(rnorm(n_rows * n_cols), ncol = n_cols)
    treat_prob <- expit((apply(x_mat, 1L, sum) / 10) *
                        abs(apply(x_mat, 1L, prod))^(1/10))

    treat_vec <- rbinom(length(treat_prob), 1L, treat_prob)

    in_sample <- propensity_score(
        x_mat,
        treat_vec,
        propensity_list = gen_propensity_list(oos_propensity = FALSE))
    out_sample_5 <- propensity_score(
        x_mat,
        treat_vec,
        propensity_list = gen_propensity_list(oos_propensity = TRUE,
                                              n_folds = 5L))
    out_sample_10 <- propensity_score(
        x_mat,
        treat_vec,
        propensity_list = gen_propensity_list(oos_propensity = TRUE,
                                              n_folds = 10L))

    brier_vec <- c(calc_brier(in_sample, treat_vec),
                   calc_brier(out_sample_5, treat_vec),
                   calc_brier(out_sample_10, treat_vec))
    expect_equal(rank(brier_vec)[1], 1L)

    expect_true(brier_vec[2] + brier_vec[3] < 0.55)

    ##------------------------------------

    x_mat <- matrix(rnorm(n_rows * n_cols), ncol = n_cols)
    treat_prob <- expit(ifelse(x_mat[, 1] > 0,
                               x_mat[, 2] * 3,
                               x_mat[, 3] * 3))

    treat_vec <- rbinom(length(treat_prob), 1L, treat_prob)

    in_sample <- propensity_score(
        x_mat,
        treat_vec,
        propensity_list = gen_propensity_list(oos_propensity = FALSE))
    out_sample_5 <- propensity_score(
        x_mat,
        treat_vec,
        propensity_list = gen_propensity_list(oos_propensity = TRUE,
                                              n_folds = 5L))
    out_sample_10 <- propensity_score(
        x_mat,
        treat_vec,
        propensity_list = gen_propensity_list(oos_propensity = TRUE,
                                              n_folds = 10L))
    brier_vec <- c(calc_brier(in_sample, treat_vec),
                   calc_brier(out_sample_5, treat_vec),
                   calc_brier(out_sample_10, treat_vec))
    expect_equal(rank(brier_vec)[1], 1L)

    expect_true(brier_vec[2] + brier_vec[3] < 0.3)
})


test_that("propensity_score_linear", {
    n_rows <- 1000L
    n_cols <- 10L

    train_test_list <- list(
        x_train = matrix(rnorm(n_rows * n_cols), ncol = n_cols),
        x_test = matrix(rnorm(n_rows * n_cols), ncol = n_cols),
        y_train = rbinom(n_rows, 1, prob = 0.5),
        y_test = rbinom(n_rows, 1, prob = 0.5)
    )

    glm_pred <- propensity_score_linear(train_test_list,
                                        use_lm_approx = FALSE)
    lm_pred <- propensity_score_linear(train_test_list,
                                       use_lm_approx = FALSE)

    expect_equal(length(glm_pred), n_rows)
    expect_equal(length(lm_pred), n_rows)

    expect_true(all(glm_pred >= 0 & glm_pred <= 1))
    expect_true(all(lm_pred >= 0 & lm_pred <= 1))


    expect_true(abs(mean(glm_pred) - 0.5) + abs(mean(lm_pred) - 0.5) < 0.1)
    expect_true(calc_brier(glm_pred, train_test_list[["y_test"]]) < 0.33 &&
                calc_brier(glm_pred, train_test_list[["y_test"]]) > 0.2)
    expect_true(calc_brier(lm_pred, train_test_list[["y_test"]]) < 0.33 &&
                calc_brier(lm_pred, train_test_list[["y_test"]]) > 0.2)

    ##----------------

    x_mat <- matrix(rnorm(n_rows * n_cols), ncol = n_cols)
    coef_vec <- runif(n_cols, -1, 1)
    treat_prob <- expit(x_mat %*% coef_vec)

    treat_vec <- rbinom(length(treat_prob), 1L, treat_prob)
    train_ind <- sample(c(TRUE, FALSE), n_rows,
                        replace = TRUE, prob = c(0.7, 0.3))

    train_test_list <- list(
        x_train = x_mat[train_ind, ],
        x_test = x_mat[!train_ind, ],
        y_train = treat_vec[train_ind],
        y_test = treat_vec[!train_ind]
    )

    glm_pred <- propensity_score_linear(train_test_list,
                                        use_lm_approx = FALSE)
    lm_pred <- propensity_score_linear(train_test_list,
                                       use_lm_approx = TRUE)

    expect_true(calc_brier(glm_pred, train_test_list[["y_test"]]) < 0.2)
    expect_true(calc_brier(lm_pred, train_test_list[["y_test"]]) < 0.2)

    ##------

    treat_prob <- x_mat[, 1]^2 / (x_mat[, 1]^2 + (0.5) * x_mat[, 2]^2)
    treat_vec <- rbinom(length(treat_prob), 1L, treat_prob)
    train_ind <- sample(c(TRUE, FALSE), n_rows,
                        replace = TRUE, prob = c(0.7, 0.3))

    train_test_list <- list(
        x_train = x_mat[train_ind, ],
        x_test = x_mat[!train_ind, ],
        y_train = treat_vec[train_ind],
        y_test = treat_vec[!train_ind]
    )

    glm_pred <- propensity_score_linear(train_test_list,
                                        use_lm_approx = FALSE)
    lm_pred <- propensity_score_linear(train_test_list,
                                       use_lm_approx = TRUE)

    expect_true(calc_brier(glm_pred, train_test_list[["y_test"]]) < 0.3)
    expect_true(calc_brier(lm_pred, train_test_list[["y_test"]]) < 0.3)

})


test_that("propensity_score_xgb", {
    n_rows <- 1000L
    n_cols <- 10L

    train_test_list <- list(
        x_train = matrix(rnorm(n_rows * n_cols), ncol = n_cols),
        x_test = matrix(rnorm(n_rows * n_cols), ncol = n_cols),
        y_train = rbinom(n_rows, 1, prob = 0.5),
        y_test = rbinom(n_rows, 1, prob = 0.5)
    )

    xgb_pred <- propensity_score_xgb(train_test_list)

    expect_equal(length(xgb_pred), n_rows)

    expect_true(all(xgb_pred >= 0 & xgb_pred <= 1))


    expect_true(abs(mean(xgb_pred) - 0.5) < 0.05)
    expect_true(calc_brier(xgb_pred, train_test_list[["y_test"]]) < 0.33 &&
                calc_brier(xgb_pred, train_test_list[["y_test"]]) > 0.2)

    ##----------------

    x_mat <- matrix(rnorm(n_rows * n_cols), ncol = n_cols)
    coef_vec <- runif(n_cols, -1, 1)
    treat_prob <- expit(x_mat %*% coef_vec)

    treat_vec <- rbinom(length(treat_prob), 1L, treat_prob)
    train_ind <- sample(c(TRUE, FALSE), n_rows,
                        replace = TRUE, prob = c(0.7, 0.3))

    train_test_list <- list(
        x_train = x_mat[train_ind, ],
        x_test = x_mat[!train_ind, ],
        y_train = treat_vec[train_ind],
        y_test = treat_vec[!train_ind]
    )

    xgb_pred <- propensity_score_xgb(train_test_list)

    expect_true(calc_brier(xgb_pred, train_test_list[["y_test"]]) < 0.25)

    ##------

    treat_prob <- x_mat[, 1]^2 / (x_mat[, 1]^2 + (0.5) * x_mat[, 2]^2)
    treat_vec <- rbinom(length(treat_prob), 1L, treat_prob)
    train_ind <- sample(c(TRUE, FALSE), n_rows,
                        replace = TRUE, prob = c(0.7, 0.3))

    train_test_list <- list(
        x_train = x_mat[train_ind, ],
        x_test = x_mat[!train_ind, ],
        y_train = treat_vec[train_ind],
        y_test = treat_vec[!train_ind]
    )

    xgb_pred <- propensity_score_xgb(train_test_list)

    expect_true(calc_brier(xgb_pred, train_test_list[["y_test"]]) < 0.2)
})
