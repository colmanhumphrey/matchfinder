context("testing match_predict")


test_that("testing match_predict functions", {
    x_mat <- cbind(rnorm(400), runif(400))
    treat_ind <- rank(rowSums(x_mat) + runif(400, 0, 2),
                      ties.method = "random") <= 200

    ## fake matched pairs, but should be predictable:
    index_list <- list(
        treat_train = which(treat_ind)[1:120],
        control_train = which(!treat_ind)[1:120],
        treat_test = which(treat_ind)[121:200],
        control_test = which(!treat_ind)[121:200]
    )

    ##------------------------------------
    ## xgb

    train_test_list <- predict_prepare(x_mat,
                                       index_list,
                                       "cross_all")
    test_pred_func <- match_predict_xgb(nrounds = 10,
                                        nthread = 1,
                                        params = list(eta = 0.3, max.depth = 5))
    test_pred <- test_pred_func(train_test_list)

    expect_equal(length(test_pred), length(train_test_list[["y_test"]]))
    ## should get most correct!
    test_pred <- ifelse(train_test_list[["y_test"]] == 1,
                        test_pred, 1 - test_pred)
    expect_true(mean(test_pred > 0.5) > 0.7)

    ##------------------------------------
    ## linear

    train_test_list <- predict_prepare(x_mat,
                                       index_list,
                                       "differences_random")
    test_pred_func <- match_predict_linear(use_linear_lm = FALSE)
    test_pred <- test_pred_func(train_test_list)

    expect_equal(length(test_pred), length(train_test_list[["y_test"]]))

    ## should get most correct!
    test_pred <- ifelse(train_test_list[["y_test"]] == 1,
                        test_pred, 1 - test_pred)
    expect_true(mean(test_pred > 0.5) > 0.7)
})
