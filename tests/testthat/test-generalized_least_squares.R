test_that("testing gls_svd_lm", {
    ## ------------------------------------
    ## set up base data
    n_rows <- 300L
    beta_vec <- c(0.0, 0.3, -0.4)
    y_sigma <- 1.5

    ## We don't use intercepts in `gls_action` or `gls_svd_lm`,
    ## that's why we add one manually
    ## else it would be tricky to specify that you're not using one
    x_mat <- cbind(
        rep(1L, n_rows),
        rnorm(n_rows),
        runif(n_rows)
    )

    y_fit <- c(x_mat %*% beta_vec)

    ## ------------------------------------
    ## first, we'll check that we get the
    ## same result if we have an identity `rel_var_mat`
    y_vector <- y_fit + rnorm(n_rows, sd = y_sigma)

    plain_res <- gls_svd_lm(
        y_vector = y_vector,
        x_mat = x_mat,
        rel_var_mat = diag(rep(1.0, n_rows))
    )

    plain_lm <- lm(y_vector ~ x_mat + 0)

    expect_true(inherits(plain_res, "lm"))

    expect_true(sum(abs(sigma(plain_res) - sigma(plain_lm))) < 0.00001)
    expect_true(sum(abs(vcov(plain_res) - vcov(plain_lm))) < 0.00001)
    expect_true(sum(abs(coefficients(plain_res) -
                        coefficients(plain_lm))) < 0.00001)

    expect_true(abs((sigma(plain_lm) - y_sigma) / y_sigma) < 0.2)

    ## ------------------------------------
    ## can we replicate with weights?

    weight_vec <- runif(n_rows, min = 1.0, max = 3.0)
    rel_var_mat <- diag(weight_vec)
    sqrt_var_mat <- diag(sqrt(weight_vec))

    plain_noise <- rnorm(n_rows, sd = y_sigma)
    rel_noise <- c(sqrt_var_mat %*% plain_noise)

    y_vector_w <- y_fit + rel_noise

    weight_res <- gls_svd_lm(
        y_vector = y_vector_w,
        x_mat = x_mat,
        rel_var_mat = rel_var_mat
    )

    ## Ok why does this work?
    ## in GLS, we are solving something like:
    ## (y - Xb)' inv(Sigma) (y - Xb)
    ## In weighted reg, we solve
    ## (y - Xb)' W (y - Xb)
    ## so we should have W prop to inv(Sigma)
    ## here they're both diag, so only the diag matters
    ## or conceptually: we're adding variance
    weight_lm <- lm(y_vector_w ~ x_mat + 0,
                    weights = 1 / weight_vec)

    expect_true(inherits(weight_res, "lm"))

    expect_true(sum(abs(sigma(weight_res) - sigma(weight_lm))) < 0.00001)
    expect_true(sum(abs(vcov(weight_res) - vcov(weight_lm))) < 0.00001)
    expect_true(sum(abs(coefficients(weight_res) -
                        coefficients(weight_lm))) < 0.00001)

    expect_true(abs((sigma(weight_res) - y_sigma) / y_sigma) < 0.2)

    ## ------------------------------------
    ## and finally a complex case

    rand_mat <- matrix(
        runif(n_rows * n_rows, -0.1, 0.1),
        nrow = n_rows
    )

    sym_eigen <- eigen(t(rand_mat) %*% rand_mat, symmetric = TRUE)
    rand_diag <- runif(n_rows, 1, 3)

    half_cov <- sym_eigen$vectors %*%
        diag(sqrt(rand_diag))
    full_cov <- half_cov %*% t(half_cov)

    rand_cov <- sym_eigen$vectors %*%
        diag(rand_diag) %*%
        t(sym_eigen$vectors)

    noise_vec <- c(half_cov %*% rnorm(n_rows, sd = y_sigma))

    y_vector_r <- y_fit + noise_vec

    rand_reg <- gls_svd_lm(
        y_vector = y_vector_r,
        x_mat = x_mat,
        rel_var_mat = full_cov
    )

    ## we'll compare against a "simple approx":
    ## we'll maintain the diag, but no correlations
    simple_rand_lm <- lm(y_vector_r ~ x_mat + 0,
                         weights = 1 / rand_diag)

    expect_true(sum(abs(sigma(rand_reg) - sigma(simple_rand_lm))) < 0.5)

    ## let's compare coefs... and ignore intercept
    coef_rand <- coefficients(rand_reg)[2:3]
    coef_lm <- coefficients(simple_rand_lm)[2:3]

    mean_diff <- mean(abs(coef_rand - coef_lm))

    expect_true(mean_diff < 0.5)

    expect_true(abs((sigma(rand_reg) -
                     y_sigma) / y_sigma) < 0.2)
})

test_that("testing gls_action", {
    ## ------------------------------------
    ## set up base data
    n_rows <- 300L
    beta_vec <- c(0.0, 0.3, -0.4)
    y_sigma <- 1.5

    ## We don't use intercepts in `gls_action` or `gls_svd_lm`,
    ## that's why we add one manually
    ## else it would be tricky to specify that you're not using one
    x_mat <- cbind(
        rep(1L, n_rows),
        rnorm(n_rows),
        runif(n_rows)
    )

    y_fit <- c(x_mat %*% beta_vec)

    ## ------------------------------------
    ## can repeat a lot of the above
    ## ------------------------------------

    ## first, we'll check that we get the
    ## same result if we have an identity `rel_var_mat`
    y_vector <- y_fit + rnorm(n_rows, sd = y_sigma)

    plain_res <- gls_action(
        y_vector = y_vector,
        x_mat = x_mat,
        rel_var_mat = diag(rep(1.0, n_rows))
    )

    plain_lm <- lm(y_vector ~ x_mat + 0)

    expect_true(sum(abs(plain_res[["sigma_sq_est"]] -
                        sigma(plain_lm)^2)) < 0.00001)
    expect_true(sum(abs(plain_res[["var_beta_gls"]] -
                        vcov(plain_lm))) < 0.00001)
    expect_true(sum(abs(plain_res[["beta_gls"]] -
                        coefficients(plain_lm))) < 0.00001)

    ## pointless test initially really, just testing
    ## `lm`, but more useful later, so stay consistent!
    expect_true(abs((sqrt(plain_res[["sigma_sq_est"]]) -
                     y_sigma) / y_sigma) < 0.2)

    ## ------------------------------------
    ## replicate with weights

    weight_vec <- runif(n_rows, min = 1.0, max = 3.0)
    rel_var_mat <- diag(weight_vec)
    sqrt_var_mat <- diag(sqrt(weight_vec))

    plain_noise <- rnorm(n_rows, sd = y_sigma)
    rel_noise <- c(sqrt_var_mat %*% plain_noise)

    y_vector_w <- y_fit + rel_noise

    weight_res <- gls_action(
        y_vector = y_vector_w,
        x_mat = x_mat,
        rel_var_mat = rel_var_mat
    )

    ## Ok why does this work?
    ## in GLS, we are solving something like:
    ## (y - Xb)' inv(Sigma) (y - Xb)
    ## In weighted reg, we solve
    ## (y - Xb)' W (y - Xb)
    ## so we should have W prop to inv(Sigma)
    ## here they're both diag, so only the diag matters
    ## or conceptually: we're adding variance
    weight_lm <- lm(y_vector_w ~ x_mat + 0,
                    weights = 1 / weight_vec)

    expect_true(sum(abs(weight_res[["sigma_sq_est"]] -
                        sigma(weight_lm)^2)) < 0.00001)
    expect_true(sum(abs(weight_res[["var_beta_gls"]] -
                        vcov(weight_lm))) < 0.00001)
    expect_true(sum(abs(weight_res[["beta_gls"]] -
                        coefficients(weight_lm))) < 0.00001)

    expect_true(abs((sqrt(weight_res[["sigma_sq_est"]]) -
                     y_sigma) / y_sigma) < 0.2)

    ## ------------------------------------
    ## and finally a complex case
    ## can now also compare to `gls_svd_lm`

    rand_mat <- matrix(
        runif(n_rows * n_rows, -0.1, 0.1),
        nrow = n_rows
    )

    sym_eigen <- eigen(t(rand_mat) %*% rand_mat, symmetric = TRUE)
    rand_diag <- runif(n_rows, 1, 3)

    half_cov <- sym_eigen$vectors %*%
        diag(sqrt(rand_diag))
    full_cov <- half_cov %*% t(half_cov)

    rand_cov <- sym_eigen$vectors %*%
        diag(rand_diag) %*%
        t(sym_eigen$vectors)

    noise_vec <- c(half_cov %*% rnorm(n_rows, sd = y_sigma))

    y_vector_r <- y_fit + noise_vec

    rand_action <- gls_action(
        y_vector = y_vector_r,
        x_mat = x_mat,
        rel_var_mat = full_cov
    )

    rand_svd <- gls_svd_lm(
        y_vector = y_vector_r,
        x_mat = x_mat,
        rel_var_mat = full_cov
    )

    ## we'll compare against a "simple approx":
    ## we'll maintain the diag, but no correlations
    simple_rand_lm <- lm(y_vector_r ~ x_mat + 0,
                         weights = 1 / rand_diag)

    expect_true(sum(abs(rand_action[["sigma_sq_est"]] -
                        sigma(simple_rand_lm)^2)) < 0.5)

    ## let's compare coefs... and ignore intercept
    coef_rand <- rand_action[["beta_gls"]][2:3]
    coef_lm <- coefficients(simple_rand_lm)[2:3]

    mean_diff <- mean(abs(coef_rand - coef_lm))

    expect_true(mean_diff < 0.5)

    expect_true(abs((sqrt(rand_action[["sigma_sq_est"]]) -
                     y_sigma) / y_sigma) < 0.2)

    ## compare svd and action

    expect_true(
        sum(abs(rand_action[["beta_gls"]] -
                coefficients(rand_svd))) < 0.00001
    )
    expect_true(
        sum(abs(rand_action[["var_beta_gls"]] -
                vcov(rand_svd))) < 0.00001
    )
    expect_true(sum(abs(rand_action[["sigma_sq_est"]] -
                        sigma(rand_svd)^2)) < 0.00001)
})

test_that("testing gls_known_relative_variance", {
    ## only really have to test that it catches errors,
    ## `gls_action` does the hard work

    ## ------------------------------------
    ## set up some example data
    n_rows <- 300L
    beta_vec <- c(0.0, 0.3, -0.4)
    y_sigma <- 1.5

    x_mat <- cbind(
        rep(1L, n_rows),
        rnorm(n_rows),
        runif(n_rows)
    )

    y_fit <- c(x_mat %*% beta_vec)

    y_vector <- y_fit + rnorm(n_rows, sd = y_sigma)

    base_var_mat <- diag(rep(1.0, n_rows))

    ## ------------------------------------
    ## testing y

    ## can't be matrix
    expect_error(
        gls_known_relative_variance(
            y_vector = as.matrix(y_vector),
            x_mat = x_mat,
            rel_var_mat = base_var_mat
        )
    )

    ## needs to match var dim
    expect_error(
        gls_known_relative_variance(
            y_vector = y_vector[1L:10L],
            x_mat = x_mat,
            rel_var_mat = base_var_mat
        )
    )

    expect_error(
        gls_known_relative_variance(
            y_vector = y_vector,
            x_mat = x_mat,
            rel_var_mat = base_var_mat
        ),
        NA
    )

    ## ------------------------------------
    ## testing x_mat

    ## dims should match y
    expect_error(
        gls_known_relative_variance(
            y_vector = y_vector,
            x_mat = x_mat[1L:10L, ],
            rel_var_mat = base_var_mat
        )
    )

    ## if vector, length should match y
    expect_error(
        gls_known_relative_variance(
            y_vector = y_vector,
            x_mat = x_mat[1L:10L, 1L],
            rel_var_mat = base_var_mat
        )
    )
    expect_error(
        gls_known_relative_variance(
            y_vector = y_vector,
            x_mat = x_mat[, 1L],
            rel_var_mat = base_var_mat
        ),
        NA
    )

    expect_error(
        gls_known_relative_variance(
            y_vector = y_vector,
            x_mat = x_mat,
            rel_var_mat = base_var_mat
        ),
        NA
    )

    ## ------------------------------------
    ## testing rel_var_mat

    ## should be a matrix
    expect_error(
        gls_known_relative_variance(
            y_vector = y_vector,
            x_mat = x_mat,
            rel_var_mat = c(base_var_mat)
        )
    )
    expect_error(
        gls_known_relative_variance(
            y_vector = y_vector,
            x_mat = x_mat,
            rel_var_mat = base_var_mat[, 1L]
        )
    )

    ## needs to be square
    expect_error(
        gls_known_relative_variance(
            y_vector = y_vector,
            x_mat = x_mat,
            rel_var_mat = base_var_mat[1L:10L, ]
        )
    )

    ## needs to match y length
    expect_error(
        gls_known_relative_variance(
            y_vector = y_vector,
            x_mat = x_mat,
            rel_var_mat = base_var_mat[1L:10L, 1L:10L]
        )
    )
    expect_error(
        gls_known_relative_variance(
            y_vector = y_vector[1L:10L],
            x_mat = x_mat[1L:10L, ],
            rel_var_mat = base_var_mat[1L:10L, 1L:10L]
        ),
        NA
    )

    expect_error(
        gls_known_relative_variance(
            y_vector = y_vector,
            x_mat = x_mat,
            rel_var_mat = base_var_mat
        ),
        NA
    )
})
