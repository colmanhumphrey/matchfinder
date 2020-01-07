context("mahal.R testing")


test_that("testing weighted_mahal", {
    cols <- 5
    half_rows <- 200
    rows <- half_rows * 2

    x_mat <- matrix(rnorm(rows * cols), nrow = rows)
    treat_vec <- rep(c(TRUE, FALSE), times = c(half_rows, half_rows))
    half_diff <- x_mat[treat_vec, ] - x_mat[!treat_vec, ]

    cov_x <- covariance_with_ranks(x_mat)

    ##----------------
    ## in this simple case, the mahal dist is approx
    ## just the summed squared diffs

    res <- weighted_mahal(x_mat,
                          cov_x = cov_x,
                          treat_vec = treat_vec)
    sqrt_mahals <- sqrt(apply(half_diff^2, 1, sum)) / cols
    expect_true(mean(abs(sqrt_mahals - diag(res))) < 0.1)

    res_sq <- weighted_mahal(x_mat,
                             cov_x = cov_x,
                             treat_vec = treat_vec,
                             sqrt_mahal = FALSE)
    mahals <- apply(half_diff^2, 1, sum) / cols
    expect_true(mean(abs(sqrt_mahals - diag(res))) < 0.1)

    ##----------------
    ## now let's put all the weight on the first two columns

    weight_vec <- c(rep(2, 2), rep(0.01, cols - 2))
    w_res <- weighted_mahal(x_mat,
                            cov_x = cov_x,
                            treat_vec = treat_vec,
                            weight_vec = weight_vec)

    w_sqrt_mahals <- sqrt(apply(half_diff[, 1:2]^2, 1, sum)) * 2

    expect_true(mean(abs(w_sqrt_mahals - diag(w_res))) < 0.5)

    ##------------------------------------
    ## test with ranks

    cols <- 3
    half_rows <- 1000
    rows <- half_rows * 2

    r_mat <- matrix(rnorm(rows * cols), nrow = rows)
    r_mat[, 1] <- exp(r_mat[, 1])
    r_mat[, 2] <- r_mat[, 2]^3
    colnames(r_mat) <- c("first", "second", "third")

    cov_r <- covariance_with_ranks(r_mat, rank_cols = c("first", "second"))

    treat_vec <- rep(c(TRUE, FALSE), times = c(half_rows, half_rows))
    r_res <- weighted_mahal(x_mat = ranked_x(r_mat, rank_cols = c("first", "second")),
                            cov_x = cov_r,
                            weight_vec = c(1, 1, 0.01),
                            treat_vec = treat_vec,
                            sqrt_mahal = TRUE)

    first_diff <- rank(r_mat[, 1])[treat_vec] - rank(r_mat[, 1])[!treat_vec]
    second_diff <- rank(r_mat[, 2])[treat_vec] - rank(r_mat[, 2])[!treat_vec]

    approx_diff <- sqrt((first_diff / rows)^2 + (second_diff / rows)^2)

    expect_true(mean(abs((rank(approx_diff) - rank(diag(r_res))))) < 20)

    ##------------------------------------
    ## partial index is fairly simple

    cols <- 5
    half_rows <- 200
    rows <- half_rows * 2

    x_mat <- matrix(rnorm(rows * cols), nrow = rows)
    partial_list <- list(c(1, 10, 100),
                         201:220)
    some_diffs <- x_mat[partial_list[[1]], ] - x_mat[partial_list[[2]][1:3], ]

    cov_x <- covariance_with_ranks(x_mat)

    partial_res <- weighted_mahal(x_mat,
                                  cov_x = cov_x,
                                  partial_index = partial_list)

    sqrt_mahals <- sqrt(apply(some_diffs^2, 1, sum)) / cols
    expect_true(mean(abs(sqrt_mahals - partial_res[cbind(1:3, 1:3)])) < 0.5)

    ##------------------------------------
})
