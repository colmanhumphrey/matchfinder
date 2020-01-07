context("utils for simulation work")


test_that("testing x_from_cov", {
    random_mat <- matrix(runif(25L), 5L, 5L)
    random_cor <- cov2cor(random_mat %*% t(random_mat))

    ##----------------

    x_500 <- x_from_cov(random_cor,
                        n_rows = 500L)

    expect_equal(dim(x_500), c(500L, ncol(random_cor)))
    expect_true(mean(abs(cov(x_500) - random_cor)) < 0.25)

    some_p_values <- lapply(1:100, function(j){
        shapiro.test(x_from_cov(random_cor, n_rows = 500L)[, 1])$p.value
    })

    expect_true(mean(unlist(some_p_values) < 0.05) < 0.15)

    ##----------------

    x_50000 <- x_from_cov(random_cor,
                          n_rows = 50000L)

    expect_equal(dim(x_50000), c(50000L, ncol(random_cor)))
    expect_true(mean(abs(cov(x_50000) - random_cor)) < 0.025)

    some_p_values <- lapply(1:25, function(j){
        shapiro.test(sample(x_from_cov(random_cor,
                                       n_rows = 50000L)[, 1],
                            size = 5000L))$p.value
    })

    expect_true(mean(unlist(some_p_values) < 0.05) < 0.3)
})


test_that("testing target_mean_expit", {
    x_mat <- default_x_generator(n_rows = 500L,
                                 n_cols = 10L)

    coef_vec <- rnorm(10L)

    linear_vec <- x_mat %*% coef_vec

    target_mean <- mean(expit(linear_vec + 0.3))

    expit_res <- target_mean_expit(target_mean = target_mean,
                                   linear_vector = linear_vec)

    expect_true(mean(abs(expit_res - expit(linear_vec + 0.3))) < 0.00001)
})


test_that("testing example_treat_prob_generator", {
    x_mat <- matrix(rnorm(1000L), 100L, 10L)

    treat_vec <- example_treat_prob_generator(x_mat)

    expect_equal(length(treat_vec),
                 nrow(x_mat))

    expect_true(all(0 < treat_vec & treat_vec < 1))
})


test_that("testing example_mean_generator", {
    x_mat <- matrix(rnorm(1000L), 100L, 10L)

    mean_vec <- example_mean_generator(x_mat)

    expect_equal(length(mean_vec),
                 nrow(x_mat))
})


test_that("testing default_x_generator", {
    arb_rows <- sample(5000L:10000L, 1L)
    arb_cols <- sample(5L:25L, 1L)

    expect_equal(
        dim(default_x_generator(n_rows = arb_rows,
                            n_cols = arb_cols)),
        c(arb_rows, arb_cols))
})


test_that("testing default_error_generator", {
    arb_rows <- sample(5000L:10000L, 1L)
    some_error_vec <- default_error_generator(n_rows = arb_rows)

    expect_equal(length(some_error_vec), arb_rows)

    expect_true(abs(mean(sign(some_error_vec))) < 0.03)
})
