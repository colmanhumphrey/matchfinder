context("testing standard dev calcs for matches")


test_that("testing gen_[non]bipartite_repeated_variance", {
    cols <- 5L
    quarter_rows <- 100L
    half_rows <- quarter_rows * 2L
    rows <- half_rows * 2L

    x_mat <- matrix(rnorm(rows * cols), nrow = rows)
    treat_vec <- rep(c(TRUE, FALSE), times = c(half_rows, half_rows))

    cov_x <- covariance_with_ranks(x_mat)
    y_vector <- runif(rows) + treat_vec * 0.2

    expect_equal(gen_bipartite_repeated_variance(
        x_mat = x_mat,
        cov_x = cov_x,
        y_vector = y_vector,
        control_index = (1L:half_rows) + half_rows,
        treat_vec = treat_vec),
        0)

    ##----------------

    sub_lengths <- c(2L, 4L, 10L, 20L, 50L, 100L)
    just_80 <- rep(80L, 10L)

    control_lists <- lapply(sub_lengths, function(x){
        unique_size_sub((1L:half_rows) + half_rows, x)
    })
    control_80s <- lapply(just_80, function(x){
        unique_size_sub((1L:half_rows) + half_rows, x)
    })

    rep_variances <- unlist(lapply(control_lists, function(cl) {
        gen_bipartite_repeated_variance(x_mat = x_mat,
                                        cov_x = cov_x,
                                        y_vector = y_vector,
                                        treat_vec = treat_vec,
                                        control_index = cl)
    }))
    rep_80s <- unlist(lapply(control_80s, function(cl) {
        gen_bipartite_repeated_variance(x_mat = x_mat,
                                        cov_x = cov_x,
                                        y_vector = y_vector,
                                        treat_vec = treat_vec,
                                        control_index = cl)
    }))

    res_80 <- mean(sqrt(rep_80s))
    scale_80s <- sqrt((half_rows - 80L) / 80L) / res_80

    should_be_close <- sqrt((half_rows - sub_lengths) / sub_lengths) /
        scale_80s

    ratios <- sqrt(rep_variances) / should_be_close

    adjusted_diffs <- abs(ratios - 1) *
        sqrt(half_rows) / sqrt(half_rows - sub_lengths)

    expect_true(mean(adjusted_diffs) < 0.5)

    ##------------------------------------
    ## again but with tol vec

    x_mat <- matrix(rnorm(rows * cols), nrow = rows)
    match_list <- list(
        treat_index = 1L:half_rows,
        control_index = (1L:half_rows) + half_rows
    )
    tol_vec <- (rows:1L) * (10 / rows)

    cov_x <- covariance_with_ranks(x_mat)
    y_vector <- runif(rows) + tol_vec / 10

    expect_equal(gen_nonbipartite_repeated_variance(
        x_mat = x_mat,
        cov_x = cov_x,
        y_vector = y_vector,
        control_index = match_list[["control_index"]],
        tolerance_list = gen_tolerance_list(
            tolerance_vec = tol_vec,
            tolerance_min = 5,
            tolerance_max = 100
        )),
        0)

    sub_lengths <- c(10L, 50L, 100L, 200L, 400L)

    ## sub from full list
    control_lists <- lapply(sub_lengths, function(x){
        unique_size_sub(1L:rows, x)
    })
    control_lists <- lapply(control_lists, function(cl) {
        same_ind <- cl == match_list[["treat_index"]]
        while (any(same_ind)) {
            cl[same_ind] <- fixed_sample(1L:rows, sum(same_ind))
            same_ind <- cl == match_list[["treat_index"]]
        }
        return(cl)
    })

    test_variances <- unlist(lapply(control_lists, function(cl) {
        gen_nonbipartite_repeated_variance(
            x_mat = x_mat,
            cov_x = cov_x,
            y_vector = y_vector,
            control_index = cl,
            tolerance_list = gen_tolerance_list(
                tolerance_vec = tol_vec,
                tolerance_min = 5,
                tolerance_max = 100
            ))
    }))

    decreased <- test_variances[1L:(length(test_variances) - 1L)] >
        test_variances[2L:length(test_variances)]

    expect_true(mean(decreased) > 0.7)
})


test_that("testing bipartite_match_sd", {
    ##------------------------------------
    ## first, using cases with no repeats
    ## this should be very boring

    cols <- 5L
    quarter_rows <- 250L
    half_rows <- quarter_rows * 2L
    rows <- half_rows * 2L

    x_mat <- matrix(rnorm(rows * cols), nrow = rows)
    treat_vec <- rep(c(TRUE, FALSE), times = c(half_rows, half_rows))

    cov_x <- covariance_with_ranks(x_mat)
    y_vector <- runif(rows) + treat_vec * 0.2

    y_vector <- runif(half_rows * 2L)

    ## won't bother with stuff we won't use
    match_list <- list(
        treat_index = 1L:half_rows,
        control_index = (1L:half_rows) + half_rows
    )

    boring_sd <- sd(y_vector[match_list[["treat_index"]]] -
                    y_vector[match_list[["control_index"]]])

    match_sd <- bipartite_match_sd(x_mat = x_mat,
                                   cov_x = cov_x,
                                   y_vector = y_vector,
                                   match_list = match_list)
    expect_equal(match_sd, boring_sd)

    ## actually don't even need x
    match_sd_no_x <- bipartite_match_sd(x_mat = NULL,
                                        cov_x = NULL,
                                        y_vector = y_vector,
                                        match_list = match_list)
    expect_equal(match_sd_no_x, boring_sd)

    ##------------------------------------
    ## now with repeats

    sub_lengths <- c(2L, 10L, 50L, quarter_rows, half_rows - 1L)

    control_lists <- lapply(sub_lengths, function(x){
        unique_size_sub((1L:half_rows) + half_rows, x)
    })
    match_lists <- lapply(control_lists, function(cl) {
        list(
            treat_index = 1L:half_rows,
            control_index = cl
        )
    })

    ## need treat_vec
    expect_error(gen_bipartite_repeated_variance(
        x_mat = x_mat,
        cov_x = cov_x,
        y_vector = y_vector,
        control_index = control_lists[[2]]))

    boring_variances <- unlist(lapply(control_lists, function(cl) {
        var(y_vector[1L:half_rows] - y_vector[cl])
    }))

    rep_variances <- unlist(lapply(control_lists, function(cl) {
        gen_bipartite_repeated_variance(x_mat = x_mat,
                                        cov_x = cov_x,
                                        y_vector = y_vector,
                                        control_index = cl,
                                        treat_vec = treat_vec)
    }))

    full_variances <- unlist(lapply(match_lists, function(ml) {
        bipartite_match_sd(x_mat,
                           cov_x = cov_x,
                           y_vector = y_vector,
                           match_list = ml,
                           treat_vec = treat_vec)
    }))^2

    should_be_zero <-
        full_variances - rep_variances - boring_variances

    expect_true(sum(abs(should_be_zero)) < 0.00000001)
})


test_that("testing nonbipartite_match_sd_scaled", {
    ##------------------------------------
    ## first, using cases with no repeats
    ## this should be very boring

    cols <- 5L
    quarter_rows <- 250L
    half_rows <- quarter_rows * 2L
    rows <- half_rows * 2L

    x_mat <- matrix(rnorm(rows * cols), nrow = rows)

    cov_x <- covariance_with_ranks(x_mat)
    treat_vec <- rep(c(TRUE, FALSE), each = half_rows)
    y_vector <- runif(rows) + treat_vec * 0.2

    y_vector <- runif(half_rows * 2L)

    ## won't bother with stuff we won't use
    match_list <- list(
        treat_index = 1L:half_rows,
        control_index = (1L:half_rows) + half_rows
    )

    tol_vec <- rep(rpois(half_rows, 10), times = 2) +
        (((1L:rows) %in% match_list[["treat_index"]]) *
         sample(c(1L, 2L), rows, replace = TRUE))

    boring_sd <- sd((y_vector[match_list[["treat_index"]]] -
                     y_vector[match_list[["control_index"]]]) /
                    (tol_vec[match_list[["treat_index"]]] -
                     tol_vec[match_list[["control_index"]]]))

    match_sd <- nonbipartite_match_sd_scaled(
        x_mat = x_mat,
        cov_x = cov_x,
        y_vector = y_vector,
        match_list = match_list,
        tolerance_list = gen_tolerance_list(
            tolerance_vec = tol_vec,
            tolerance_min = 1))
    expect_equal(match_sd, boring_sd)

    ## actually don't even need x
    match_sd_no_x <- nonbipartite_match_sd_scaled(
        x_mat = NULL,
        cov_x = NULL,
        y_vector = y_vector,
        match_list = match_list,
        tolerance_list = gen_tolerance_list(
            tolerance_vec = tol_vec,
            tolerance_min = 1))
    expect_equal(match_sd_no_x, boring_sd)

    ##------------------------------------
    ## now with repeats

    sub_lengths <- c(rep(50L, 4L),
                     rep(quarter_rows, 4L))

    control_lists <- lapply(sub_lengths, function(x){
        unique_size_sub(1L:rows, x)
    })
    control_lists <- lapply(control_lists, function(cl) {
        bad_ind <- cl == match_list[["treat_index"]] |
            tol_vec[cl] == tol_vec[match_list[["treat_index"]]]
        while (any(bad_ind)) {
            cl[bad_ind] <- fixed_sample(1L:rows, sum(bad_ind))
            bad_ind <- cl == match_list[["treat_index"]] |
                tol_vec[cl] == tol_vec[match_list[["treat_index"]]]
        }
        return(cl)
    })
    match_lists <- lapply(control_lists, function(cl) {
        list(
            treat_index = 1L:half_rows,
            control_index = cl
        )
    })

    boring_variances <- unlist(lapply(control_lists, function(cl) {
        var((y_vector[1L:half_rows] - y_vector[cl]) /
            (tol_vec[1L:half_rows] - tol_vec[cl]))
    }))

    rep_variances <- unlist(lapply(control_lists, function(cl) {
        gen_nonbipartite_repeated_variance(
            x_mat = x_mat,
            cov_x = cov_x,
            y_vector = y_vector,
            control_index = cl,
            tolerance_list = gen_tolerance_list(
                tolerance_vec = tol_vec,
                tolerance_min = 1))
    }))

    full_variances <- unlist(lapply(match_lists, function(ml) {
        nonbipartite_match_sd_scaled(x_mat,
                                     cov_x = cov_x,
                                     y_vector = y_vector,
                                     match_list = ml,
                                     tolerance_list = gen_tolerance_list(
                                         tolerance_vec = tol_vec,
                                         tolerance_min = 1))
    }))^2

    boring_avg <- colMeans(matrix(boring_variances, ncol = 2))
    rep_avg <- colMeans(matrix(rep_variances, ncol = 2))
    full_avg <- colMeans(matrix(full_variances, ncol = 2))

    close_to_zero <- full_avg - rep_avg - boring_avg

    expect_true(sum(abs(close_to_zero)) < 0.5)
})
