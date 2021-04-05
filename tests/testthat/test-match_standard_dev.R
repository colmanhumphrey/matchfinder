context("testing standard dev calcs for matches")

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

    match_se <- regression_eval(match_list, y_vector)[["standard_error"]]
    match_sd <- match_se * sqrt(length(match_list[["treat_index"]]))

    expect_equal(match_sd, boring_sd)

    ##------------------------------------
    ## now with different number of repeats

    sub_lengths <- c(1L, 2L, 4L, 10L, 50L, quarter_rows, half_rows - 1L)

    se_ratios <- lapply(seq_len(10L), function(j_iter) {
        control_lists <- lapply(sub_lengths, function(x) {
            unique_size_sub((1L:half_rows) + half_rows, x)
        })
        match_lists <- lapply(control_lists, function(cl) {
            list(
                treat_index = 1L:half_rows,
                control_index = cl
            )
        })

        naive_ses <- unlist(lapply(control_lists, function(cl) {
            sd(y_vector[1L:half_rows] - y_vector[cl]) / sqrt(length(cl))
        }))

        match_ses <- unlist(lapply(match_lists, function(match_list) {
            regression_eval(match_list, y_vector = y_vector)[["standard_error"]]
        }))

        match_ses / naive_ses
    })
    se_ratio_mat <- do.call(rbind, se_ratios)

    ## this should be something like sqrt(500 / sub_lengths) worse
    ## or... sqrt(500 / sub_lengths + 1) / sqrt(2)?
    ## TODO: need to work this out in greater detail
    approx_expect <- sqrt(500 / sub_lengths + 1) / sqrt(2)
    naive_expect <- sqrt(500 / sub_lengths)
    ## hmm
    approx_weight <- (sub_lengths - 1) / (sub_lengths + 1)
    combined_expect <- approx_expect * approx_weight +
        naive_expect * (1 - approx_weight)
    ratio_diffs <- apply(se_ratio_mat, 2, mean) / combined_expect
    dist_from_one <- abs(ratio_diffs - 1)

    expect_true(mean(dist_from_one) < 0.05)
})

test_that("testing approx_ratio_sd", {
    ##------------------------------------
    ## value testing

    n <- 1000L

    some_numer <- runif(n) + rnorm(n)
    denom_list <- list(
        neg_cor_denom = runif(n, 9, 18) - 0.3 * some_numer,
        zero_cor_denom = abs(rt(n, 10) + rt(n, 20) + 1) + 2,
        low_cor_denom = 5 * rbeta(n, 3, 5) + 0.05 * some_numer + 3,
        high_cor_denom = 10 * rbeta(n, 3, 5) + 0.6 * some_numer + 3
    )

    ratio_sd <- lapply(denom_list, function(denom) {
        sd(some_numer / denom)
    })

    approx_ratios <- lapply(denom_list, function(denom) {
        approx_ratio_sd(
            mean(some_numer),
            sd(some_numer),
            mean(denom),
            sd(denom),
            cov(some_numer, denom)
        )
    })

    sym_mean_error <- lapply(names(denom_list), function(name) {
        abs(ratio_sd[[name]] - approx_ratios[[name]]) /
            (2 * (abs(ratio_sd[[name]]) + abs(approx_ratios[[name]])))
    })

    expect_true(mean(unlist(sym_mean_error)) < 0.05)

    ##------------------------------------
    ## error testing

    ## negative sds
    expect_error(
        approx_ratio_sd(
            numerator_mean = 5,
            numerator_sd = -0.1,
            denominator_mean = 10,
            denominator_sd = 1,
            covariance = 0
        )
    )

    ## denom must have pos mean
    expect_error(
        approx_ratio_sd(
            numerator_mean = 5,
            numerator_sd = 2,
            denominator_mean = -10,
            denominator_sd = 1,
            covariance = 0
        )
    )

    ## impossible cov
    expect_error(
        approx_ratio_sd(
            numerator_mean = 5,
            numerator_sd = 2,
            denominator_mean = 10,
            denominator_sd = 1,
            covariance = 100
        )
    )

    ## if denom seems to be awfully close to giving negs...:
    ## (ratio of sd to mean)
    expect_warning(
        approx_ratio_sd(
            numerator_mean = 5,
            numerator_sd = 2,
            denominator_mean = 1,
            denominator_sd = 100,
            covariance = 0
        )
    )
})

test_that("testing y_tolerance_diff_ratio_sd", {
    ##------------------------------------
    ## easy tolerance version
    treat_effect <- 0.1
    opt_test_list <- lapply(1L:200L, function(j) {
        rows <- 1000L
        some_noise <- runif(rows, -1, 1)
        x_vec_1 <- runif(rows, -2, 2) + some_noise * 0.5
        x_vec_2 <- runif(rows, -2, 2) + some_noise * 0.5
        x_mat <- cbind(x_vec_1, x_vec_2)

        tol_vec <- (1:rows) / rows + x_vec_1 * 0.1 + x_vec_2 * 0.3 -
            some_noise * 0.05
        tol_list <- gen_tolerance_list(tolerance_vec = tol_vec)

        dist_mat <- weighted_mahal(
            x_mat = ranked_x(x_mat, rank_cols = NULL),
            cov_x = covariance_with_ranks(x_mat, rank_cols = NULL)
        )

        y_vector <- treat_effect * tol_vec + 0.5 * x_vec_1 + 0.8 * x_vec_2 +
            0.3 * runif(rows) + some_noise * 0.2

        direct_reg <- lm(y_vector ~ tol_vec + x_vec_1 + x_vec_2)

        ## fit optimal
        nonbimatch <- nonbipartite_matches(
            dist_mat = dist_mat,
            tolerance_list = tol_list,
            match_method = "optimal",
            n_sinks = 0L
        )[["0"]]

        nonbi_eval <- regression_eval(
            match_list = nonbimatch,
            y_vector = y_vector,
            tolerance_list = tol_list
        )

        y_diffs <- y_vector[nonbimatch[["treat_index"]]] -
            y_vector[nonbimatch[["control_index"]]]
        tol_vec <- tol_list[["tolerance_vec"]]
        tol_diffs <- tol_vec[nonbimatch[["treat_index"]]] -
            tol_vec[nonbimatch[["control_index"]]]
        nonbi_mean_naive <- mean(y_diffs / tol_diffs)

        random_match_res <- y_tolerance_diff_ratio(
            y_vector = y_vector,
            tolerance_list = tol_list
        )
        return(list(
            match_se = nonbi_eval[["standard_error"]],
            match_mean = nonbi_eval[["estimate"]],
            match_mean_naive = nonbi_mean_naive,
            random_se_reg = random_match_res[["se"]],
            random_mean_reg = random_match_res[["mean"]],
            direct_reg_mean = unname(coef(direct_reg)["tol_vec"])
        ))
    })

    res_frame <- as.data.frame(
        do.call(rbind, lapply(opt_test_list, unlist)))

    ## OK this test wasn't really supposed to be about the estimates,
    ## but why not
    mean_percent_treat_error <- function(x) {
        sqrt(mean(abs(x - treat_effect)^2)) / treat_effect
    }
    perc_errors <- list(
        match = list(
            reg = mean_percent_treat_error(res_frame[["match_mean"]]),
            naive = mean_percent_treat_error(res_frame[["match_mean_naive"]])
        ),
        random_reg = mean_percent_treat_error(res_frame[["random_mean_reg"]]),
        actual_reg = mean_percent_treat_error(res_frame[["direct_reg_mean"]])
    )
    res_means <- colMeans(res_frame)
    perc_bias <- list(
        match = list(
            reg = mean_percent_treat_error(res_means["match_mean"]),
            naive = mean_percent_treat_error(res_means["match_mean_naive"])
        ),
        random_reg = mean_percent_treat_error(res_means["random_mean_reg"]),
        actual_reg = mean_percent_treat_error(res_means["direct_reg_mean"])
    )

    ## within 50% over the 20 runs - usually better,
    ## avg 30%
    expect_true(perc_errors[["match"]][["reg"]] < 0.50)
    ## beats the random
    expect_true(perc_errors[["match"]][["reg"]] <
                perc_errors[["random_reg"]])
    ## beats regression!
    expect_true(perc_errors[["match"]][["reg"]] <
                perc_errors[["actual_reg"]])

    ## same results for "bias",
    ## avg 20%, within 40%
    expect_true(perc_bias[["match"]][["reg"]] < 0.40)
    ## beats the random
    expect_true(perc_bias[["match"]][["reg"]] <
                perc_bias[["random_reg"]])
    ## beats regression!
    expect_true(perc_bias[["match"]][["reg"]] <
                perc_bias[["actual_reg"]])

    ## now comparing the standard deviations
    ## big decrease
    expect_true(res_means["match_se"] <
                0.4 * res_means["random_se_reg"])

    ## ------------------------------------
    ## simpler overall test, more tricky tolerance

    treat_effect <- 0.1
    opt_test_list <- lapply(1L:20L, function(j) {
        print(j)
        rows <- 300L
        some_noise <- runif(rows, -1, 1)
        x_vec <- runif(rows, -2, 2)
        x_mat <- matrix(x_vec, ncol = 1L)

        tol_vec <- 1:rows
        tol_list <- gen_tolerance_list(
            tolerance_vec = tol_vec,
            tolerance_min = rows / 3,
            tolerance_max = 2 * rows / 3
        )

        dist_mat <- weighted_mahal(
            x_mat = ranked_x(x_mat, rank_cols = NULL),
            cov_x = covariance_with_ranks(x_mat, rank_cols = NULL)
        )

        y_vector <- treat_effect * tol_vec + 1.2 * x_vec + 0.3 * rnorm(rows)

        direct_reg <- lm(y_vector ~ tol_vec + x_vec)

        ## fit optimal
        nonbimatch <- nonbipartite_matches(
            dist_mat = dist_mat,
            tolerance_list = tol_list,
            match_method = "optimal",
            n_sinks = 0L
        )[["0"]]
        nonbi_eval <- regression_eval(
            match_list = nonbimatch,
            y_vector = y_vector,
            tolerance_list = tol_list
        )

        y_diffs <- y_vector[nonbimatch[["treat_index"]]] -
            y_vector[nonbimatch[["control_index"]]]
        tol_vec <- tol_list[["tolerance_vec"]]
        tol_diffs <- tol_vec[nonbimatch[["treat_index"]]] -
            tol_vec[nonbimatch[["control_index"]]]
        nonbi_mean_naive <- mean(y_diffs / tol_diffs)

        random_match_res <- y_tolerance_diff_ratio(
            y_vector = y_vector,
            tolerance_list = tol_list
        )
        return(list(
            match_se = nonbi_eval[["standard_error"]],
            match_mean = nonbi_eval[["estimate"]],
            match_mean_naive = nonbi_mean_naive,
            random_se_reg = random_match_res[["se"]],
            random_mean_reg = random_match_res[["mean"]],
            direct_reg_mean = unname(coef(direct_reg)["tol_vec"])
        ))
    })

    res_frame <- as.data.frame(
        do.call(rbind, lapply(opt_test_list, unlist)))

    ## OK this test wasn't really supposed to be about the estimates,
    ## but why not
    mean_percent_treat_error <- function(x) {
        sqrt(mean(abs(x - treat_effect)^2)) / treat_effect
    }
    perc_errors <- list(
        match = list(
            reg = mean_percent_treat_error(res_frame[["match_mean"]]),
            naive = mean_percent_treat_error(res_frame[["match_mean_naive"]])
        ),
        random_reg = mean_percent_treat_error(res_frame[["random_mean_reg"]]),
        actual_reg = mean_percent_treat_error(res_frame[["direct_reg_mean"]])
    )
    res_means <- colMeans(res_frame)
    perc_bias <- list(
        match = list(
            reg = mean_percent_treat_error(res_means["match_mean"]),
            naive = mean_percent_treat_error(res_means["match_mean_naive"])
        ),
        random_reg = mean_percent_treat_error(res_means["random_mean_reg"]),
        actual_reg = mean_percent_treat_error(res_means["direct_reg_mean"])
    )

    ## within 1% over the 20 runs - usually better
    expect_true(perc_errors[["match"]][["reg"]] < 0.01)
    ## beats the random
    expect_true(perc_errors[["match"]][["reg"]] <
                perc_errors[["random_reg"]])
    ## here the regression is better as the setup
    ## is perfectly suited to it

    ## NOT same results for "bias"
    ## random is totally unbiased!
    expect_true(perc_bias[["match"]][["reg"]] < 0.003)

    ## now comparing the standard deviations
    expect_true(res_means["match_se"] <
                0.5 * res_means["random_se_reg"])
})

test_that("testing tol_random_sample", {
    rows <- 500L
    tol_vec <- runif(rows)

    easy_tol <- gen_tolerance_list(tol_vec)
    medium_tol <- gen_tolerance_list(
        tol_vec,
        tolerance_min = 0.1,
        tolerance_max = 0.9
    )
    hard_tol <- gen_tolerance_list(
        tol_vec,
        tolerance_min = 0.3,
        tolerance_max = 0.6
    )

    ##------------------------------------

    easy_pairs <- tol_random_sample(easy_tol)
    expect_equal(unname(lengths(easy_pairs)), rep(rows %/% 2L, 2L))

    expect_false(tolerance_check(match_list = easy_pairs,
                                 tolerance_list = easy_tol)[["error"]])

    ##------------------------------------

    medium_pair_res <- lapply(1L:100L, function(j) {
        tol_random_sample(medium_tol)
    })

    num_pairs <- unlist(lapply(medium_pair_res, function(x) {
        length(x[["treat_index"]])
    }))

    ## ballpark low 30s
    expect_true(sum(num_pairs == (rows %/% 2L)) > 5L)
    expect_true(sum(num_pairs == (rows %/% 2L)) < 80L)

    any_result <- sample(medium_pair_res, 1L)[[1L]]

    expect_false(tolerance_check(match_list = any_result,
                                 tolerance_list = medium_tol)[["error"]])

    ##------------------------------------

    hard_pair_res <- lapply(1L:100L, function(j) {
        tol_random_sample(hard_tol)
    })

    num_pairs <- unlist(lapply(hard_pair_res, function(x) {
        length(x[["treat_index"]])
    }))

    ## ballpark none
    expect_true(!any(num_pairs == (rows %/% 2L)))
    ## but shouldn't be too bad as long as rows isn't tiny
    expect_true(min(num_pairs) > (rows %/% 2L) * 0.7)

    any_result <- sample(hard_pair_res, 1L)[[1L]]

    expect_false(tolerance_check(match_list = any_result,
                                 tolerance_list = hard_tol)[["error"]])
})
