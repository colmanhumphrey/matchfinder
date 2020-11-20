test_that("regression_eval for simple cases", {
    half_rows <- 50L
    n_rows <- 2 * half_rows
    y_vector <- rnorm(n_rows)
    match_list <- list(
        treat_index = 1L:half_rows,
        control_index = 1L:half_rows + half_rows
    )

    y_diffs <- y_vector[match_list[["treat_index"]]] -
        y_vector[match_list[["control_index"]]]

    ## ------------------------------------
    ## bipartite

    eval_reg <- regression_eval(match_list, y_vector)

    expected_lm <- lm(y_diffs ~ 1)
    expected_res <- summary(expected_lm)[["coefficients"]][1L, ]

    expect_true(
        abs(eval_reg[["estimate"]] - expected_res[1L]) < 0.00001
    )
    expect_true(
        abs(eval_reg[["standard_error"]] - expected_res[2L]) < 0.00001
    )

    ## ------------------------------------
    ## nonbipartite

    ## fake out the diffs
    tol_diffs <- runif(half_rows, 0.3, 1.3)
    tol_vec <- c(tol_diffs, rep(0.0, half_rows))

    tol_list <- gen_tolerance_list(tolerance_vec = tol_vec)

    eval_nb_reg <- regression_eval(match_list, y_vector, tol_list)

    expected_lm <- lm(y_diffs ~ tol_diffs + 0)
    expected_res <- summary(expected_lm)[["coefficients"]][1L, ]

    expect_true(
        abs(eval_nb_reg[["estimate"]] - expected_res[1L]) < 0.00001
    )
    expect_true(
        abs(eval_nb_reg[["standard_error"]] - expected_res[2L]) < 0.00001
    )
})

test_that("regression_eval error testing", {
    half_rows <- 50L
    n_rows <- 2 * half_rows
    y_vector <- rnorm(n_rows)
    match_list <- list(
        treat_index = 1L:half_rows,
        control_index = 1L:half_rows + half_rows
    )

    ## base case is good
    expect_error(
        regression_eval(match_list, y_vector),
        NA
    )

    y_with_nas <- c(y_vector[1L:(n_rows - 2L)], NA, NA)
    expect_error(
        regression_eval(match_list, y_with_nas)
    )

    match_list_too_much <- match_list
    match_list_too_much[["treat_index"]][half_rows] <- n_rows + 10L
    expect_error(
        regression_eval(match_list_too_much, y_vector)
    )
})

test_that("regression_eval with repeated controls", {
    num_neg <- 20L
    num_half <- num_neg * 3L

    y_diffs <- c(
        sample(seq(0.5, 1.5, length.out = num_neg * 2L)),
        sample(seq(-1.5, -0.5, length.out = num_neg))
    )

    ## we want no noise here because we want to test
    ## sharing controls, but not having the control value
    ## differ
    noise_vector <- rep(0, length(y_diffs))
    y_vector <- c(y_diffs + noise_vector, noise_vector)

    match_list <- list(
        treat_index = seq_len(length(y_diffs)),
        control_index = seq_len(length(y_diffs)) + length(y_diffs)
    )

    base_result <- regression_eval(match_list, y_vector)

    one_controls <- seq_len(num_neg * 2L) + num_half
    neg_controls <- seq_len(num_neg * 1L) + num_half + num_neg * 2L

    ## ------------------------------------
    ## this will seem a bit fuzzy only because
    ## we're not doing any matching on anything useful,
    ## So we'll try to get fun results anyway.

    ## a comparison case: throw data away (but "balanced")
    ## stopping at 10 because that's where `repeat` has to stop
    throw_results <- lapply(seq_len(10L), function(toss) {
        throw_some <- lapply(seq_len(20L), function(j) {
            match_list_throw <- match_list
            ## toss four of the ones, two of the neg ones
            match_list_throw[["treat_index"]] <- c(
                sample(seq_len(num_neg * 2L), size = num_neg * 2L - toss * 2L),
                num_neg * 2L + sample(seq_len(num_neg), size = num_neg - toss)
            )
            match_list_throw[["control_index"]] <-
                match_list_throw[["treat_index"]] + num_half

            regression_eval(match_list_throw, y_vector)
        })
        res_df <- do.call(rbind.data.frame, throw_some)
        res_df[["n"]] <- toss
        res_df
    })

    agg_func <- function(res_df) {
        data.frame(
            min_est = min(res_df$estimate),
            mean_est = mean(res_df$estimate),
            max_est = max(res_df$estimate),
            min_se = min(res_df$standard_error),
            mean_se = mean(res_df$standard_error),
            max_se = max(res_df$standard_error),
            n = res_df$n[1]
        )
    }

    throw_agg <- do.call(rbind.data.frame, lapply(throw_results, agg_func))

    ## use some controls twice
    ## can only go to ten here, else too few
    some_twice <- lapply(seq_len(10L), function(two_count) {
        double_control <- lapply(seq_len(20L), function(j) {
            match_two_control <- match_list
            retain_one <- sample(seq_len(num_neg * 2L),
                                 size = num_neg * 2L - two_count * 2L)
            retain_neg <- num_neg * 2L +
                sample(seq_len(num_neg), size = num_neg - two_count)

            ## why the extra sample? to ensure we don't always pair
            ## the "first" set of values with the "last"
            match_two_control[["control_index"]] <- num_half + c(
                c(sample(c(retain_one, sample(retain_one, two_count * 2L))),
                  sample(c(retain_neg, sample(retain_neg, two_count)))
                  )
            )

            ## almost used this:
            ## but this is a huge problem: this will actually
            ## always pair the first treat with the 23rd, the second
            ## with the 24th... etc
            ## so e.g. if the controls are constant say,
            ## you'll actually just get the same result repeated
            ## and if they're not, you'll falsely reduce the between-run
            ## variance. Tricky!
            ## match_two_control[["control_index"]] <- num_half + c(
            ##     c(retain_one, retain_one[(two_count * 2L):1L],
            ##       retain_neg, retain_neg[two_count:1L])
            ## )

            regression_eval(match_two_control, y_vector)
        })
        res_df <- do.call(rbind.data.frame, double_control)
        res_df[["n"]] <- two_count
        res_df
    })
    repeat_agg <- do.call(rbind.data.frame, lapply(some_twice, agg_func))

    ## ----------------

    perc_diff <- function(value, target) {
        abs(value - target) / target
    }

    expect_true(mean(perc_diff(throw_agg$mean_est,
                               base_result$estimate)) < 0.05)
    expect_true(mean(perc_diff(repeat_agg$mean_est,
                               base_result$estimate)) < 0.05)

    throw_se_diff <- perc_diff(throw_agg$mean_se,
                               base_result$standard_error)
    repeat_se_diff <- perc_diff(repeat_agg$mean_se,
                                base_result$standard_error)

    ## I think under uniformity should bounce around 1/3
    rank_diffs <- function(values) {
        diffs <- rank(values) - seq_len(length(values))
        mean(abs(diffs)) / length(values)
    }
    expect_true(rank_diffs(throw_se_diff) < 0.1)
    expect_equal(rank_diffs(repeat_se_diff) < 0.1)

    ## probably all true
    expect_true(mean(repeat_agg$mean_se <
                     throw_agg$mean_se[seq_len(nrow(repeat_agg))]) > 0.8)
})
