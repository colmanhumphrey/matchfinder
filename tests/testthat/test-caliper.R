context("testing caliper adding")


test_that("testing gen_caliper_list", {
    expect_equal(NULL,
                 gen_caliper_list())

    expect_equal(NULL,
                 gen_caliper_list(caliper_vec = NULL))

    expect_error(gen_caliper_list(caliper_max = 3))
    expect_error(gen_caliper_list(caliper_vec = c(1, 2, 10)))
    expect_error(gen_caliper_list(caliper_vec = c(1, 2, 10),
                                  caliper_max = c(3, 10)))
    expect_error(gen_caliper_list(caliper_vec = c(1, 2, 10),
                                  caliper_max = c(3),
                                  continuous_mult = c(20, 1)))

    caliper_vec <- runif(10)
    caliper_max <- 0.2
    continuous_mult <- 10

    expect_equal(
        list(caliper_vec = caliper_vec,
             caliper_max = caliper_max,
             continuous_mult = continuous_mult),
        gen_caliper_list(caliper_vec, caliper_max, continuous_mult)
    )

    expect_equal(
        list(caliper_vec = caliper_vec,
             caliper_max = caliper_max,
             continuous_mult = 100),
        gen_caliper_list(caliper_vec, caliper_max)
    )
})


test_that("failure modes", {
    ##------------------------------------
    ## need caliper_list
    expect_error(create_caliper())

    ##------------------------------------
    ## can't have treat_vec different length than caliper_vec
    expect_error(create_caliper(gen_caliper_list(caliper_vec = c(1, 2, 3),
                                                 caliper_max = 1),
                                treat_vec = c(1, 1, 0, 0)))
    expect_error(create_caliper(gen_caliper_list(caliper_vec = c(1, 2, 3),
                                                 caliper_max = 1),
                                treat_vec = c(1, 1, 0)),
                 NA)
})


test_that("success modes", {
    ##------------------------------------
    ## with Inf mult (NULL same behaviour), the default
    expect_equal(create_caliper(gen_caliper_list(caliper_vec = c(1, 2, 3),
                                                 caliper_max = 1,
                                                 continuous_mult = Inf)),
                 matrix(c(0, 0, Inf, 0, 0, 0, Inf, 0, 0),
                        nrow = 3))
    expect_equal(create_caliper(gen_caliper_list(caliper_vec = c(1, 2, 3),
                                                 caliper_max = 1,
                                                 continuous_mult = NULL),
                                treat_vec = c(1, 1, 0)),
                 matrix(c(Inf, 0), nrow = 2))

    ##------------------------------------
    ## with cont. mult

    result_mat <- create_caliper(
        gen_caliper_list(caliper_vec = 1:10,
                         caliper_max = 3,
                         continuous_mult = 7),
        treat_vec = rep(c(TRUE, FALSE), times = c(5, 5)))

    expect_equal(result_mat[cbind(
        c(1, 2, 3), c(3, 4, 5))],
        rep(28, 3))
})
