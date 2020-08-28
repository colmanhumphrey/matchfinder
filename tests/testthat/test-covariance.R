context("covariance.R testing")


test_that("testing covariance_with_ranks", {
    n <- 5
    rand_mat <- matrix(runif(n^2, -1, 1), n, n)
    sigma <- rand_mat %*% t(rand_mat)
    mean_vec <- rnorm(n, 1, 2)

    rows <- 100000
    z_mat <- matrix(rnorm(rows * n), ncol = n)
    x_mat <- z_mat %*% t(rand_mat) + matrix(mean_vec,
                                            nrow = rows,
                                            ncol = n,
                                            byrow = TRUE)

    ## should be just sigma, i.e. just calls plain cov(x_mat)
    x_cov <- covariance_with_ranks(x_mat)

    expect_true(sum(abs(x_cov - sigma)) < 1)

    ##------------------------------------
    ## with ranks

    n <- 5
    rand_mat <- matrix(runif(n^2, -1, 1), n, n)
    sigma <- rand_mat %*% t(rand_mat)
    mean_vec <- rnorm(n, 1, 2)

    rows <- 100000
    z_mat <- matrix(rnorm(rows * n), ncol = n)
    x_mat <- z_mat %*% t(rand_mat) + matrix(mean_vec,
                                            nrow = rows,
                                            ncol = n,
                                            byrow = TRUE)
    ## the values have some unknown rank-preserving function applied
    x_mat[, 2] <- exp(x_mat[, 2])
    x_mat[, 3] <- x_mat[, 3]^3
    x_mat[, 4] <- exp(x_mat[, 4] * 2 + 3)

    colnames(x_mat) <- as.character(1:5)

    x_cov_rank <- covariance_with_ranks(x_mat,
                                        rank_cols = c("2", "3", "4"))
    x_cov_raw <- covariance_with_ranks(x_mat)

    expect_true(sum(abs(x_cov_rank - sigma)) < 100)

    ## a 20x (usually much bigger) improvement over the raw version:
    expect_true(sum(abs(x_cov_rank - sigma)) <
                (sum(abs(x_cov_raw - sigma)) / 20))

    ##------------------------------------
    ## error modes

    ## must give a matrix
    expect_error(covariance_with_ranks())
    expect_error(covariance_with_ranks(x_mat = c(1, 2, 3)))

    ## if rank_cols given, x_mat must have colnames
    expect_error(covariance_with_ranks(x_mat = matrix(1:4, 2, 2),
                                       rank_cols = "some_col"))
    ## if x_mat has colnames, rank_cols must be a subset
    expect_error(covariance_with_ranks(
        x_mat = matrix(1:4, 2, 2,
                       dimnames = list(c(), c("a", "b"))),
        rank_cols = "some_col"))
    ## no error here
    expect_error(covariance_with_ranks(
        x_mat = matrix(1:4, 2, 2,
                       dimnames = list(c(), c("a", "b"))),
        rank_cols = "b"),
        NA)

    ##------------------------------------
})
