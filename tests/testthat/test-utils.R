context("testing utilities")


test_that("test min_different_rank", {
    ##------------------------------------
    ## with no ties or any other issues, should be simple
    rows <- 20L
    dist_mat <- matrix(rnorm(rows * rows), nrow = rows)
    min_inds <- min_different_rank(dist_mat)
    min_from_apply <- apply(dist_mat, 1, which.min)
    equal_to_ind <- min_from_apply == 1L:rows
    expect_equal(min_inds[!equal_to_ind],
                 min_from_apply[!equal_to_ind])

    ##------------------------------------
    dist_mat <- rbind(
        c(0, 100, 200),
        c(0, -5, 3),
        c(1, 3, -1)
    )
    min_inds <- min_different_rank(dist_mat)
    expect_equal(min_inds, c(2L, 1L, 1L))
    ## compare:
    expect_equal(apply(dist_mat, 1, which.min),
                 c(1L, 2L, 3L))
})


test_that("test min_blocked_rank", {
    ##------------------------------------
    ## with no ties or any other issues, should be simple
    rows <- 20L
    cols <- 50L
    blocked_ind <- (1L:20L) * 2L
    dist_mat <- matrix(rnorm(rows * cols), nrow = rows)
    min_inds <- min_blocked_rank(dist_mat, blocked_ind)

    min_from_apply <- apply(dist_mat, 1, which.min)
    equal_to_ind <- min_from_apply == blocked_ind
    expect_equal(min_inds[!equal_to_ind],
                 min_from_apply[!equal_to_ind])

    ##------------------------------------
    dist_mat <- rbind(
        c(0, 100, 200),
        c(0, 2, 3),
        c(1, 3, 4)
    )
    min_inds <- min_blocked_rank(dist_mat, c(1L, 1L, 1L))
    expect_equal(min_inds, c(2L, 2L, 2L))
    ## compare:
    expect_equal(apply(dist_mat, 1, which.min),
                 c(1L, 1L, 1L))
})


test_that("test fixed_sample", {
    some_vector <- sample(1:100, 50)
    expect_equal(sort(fixed_sample(some_vector)),
                 sort(some_vector))

    ## and the purpose:
    expect_equal(fixed_sample(10L),
                 10L)
    expect_true(length(fixed_sample(10L)) == 1L)
    ## and in contrast:
    expect_true(length(sample(10L)) == 10L)
})


test_that("test unique_size_sub", {
    base_vector <- sample(1L:20000L, 500L)

    random_400 <- unique_size_sub(base_vector, 400L)

    expect_equal(length(random_400), 500L)
    expect_equal(length(unique(random_400)), 400L)
})


test_that("test is_tf", {
    expect_true(is_tf(TRUE))
    expect_true(is_tf(FALSE))

    expect_false(is_tf("hello"))
    expect_false(is_tf(1))

    expect_false(is_tf(c(TRUE, FALSE)))
})


test_that("testing fold_indexing", {
    nine_five <- fold_indexing(9, 5)
    expect_equal(1L:9L, sort(unlist(nine_five)))

    expect_true(all(lengths(nine_five) %in% c(1L, 2L)))

    ##----------------

    ten_five <- fold_indexing(10, 5)
    expect_equal(1L:10L, sort(unlist(ten_five)))

    expect_true(all(lengths(ten_five) == 2L))
})


test_that("testing calc_brier", {
    boring_res <- calc_brier(rep(0.5, 10L),
                             rbinom(10, 1, runif(1)),
                             avg = TRUE)
    expect_equal(boring_res, 0.25)

    boring_sum_res <- calc_brier(rep(0.5, 10L),
                             rbinom(10, 1, runif(1)),
                             avg = FALSE)
    expect_equal(boring_sum_res, 2.5)

    ##------------------------------------

    rand_res <- calc_brier(runif(1000),
                           rbinom(1000, 1, runif(1)),
                           avg = TRUE)
    expect_true(abs(rand_res - 1 / 3) < 0.04)

    ##------------------------------------

    pred_vec <- runif(1000)
    outcome_vec <- ifelse(pred_vec > 0.5, 1, 0)

    ## should be solid
    expect_true(calc_brier(pred_vec, outcome_vec, avg = TRUE) < 0.15)
})


test_that("testing swap_pairs", {
    match_list <- list(
        treat_index = 1L:20L,
        treat_index_within = 1L:20L,
        control_index = 21L:40L,
        control_index_within = 1L:20L,
        distance = runif(20)
    )

    ##------------------------------------

    expect_equal(swap_pairs(match_list, rep(FALSE, 20L)),
                 match_list)

    expect_equal(swap_pairs(match_list,
                            sample(c(TRUE, FALSE), 20L, TRUE))[["distance"]],
                 match_list[["distance"]])

    ##------------------------------------

    swaps <- (1L:20L) %in% sample(1L:20L, size = 10L)

    res_list <- swap_pairs(match_list, swaps)

    expect_equal(res_list[["treat_index"]][!swaps],
                 match_list[["treat_index"]][!swaps])
    expect_equal(res_list[["treat_index"]][swaps],
                 match_list[["control_index"]][swaps])
})


test_that("testing rank_integer_index", {
    x_mat <- matrix(runif(30), ncol = 3)

    expect_equal(rank_integer_index(NULL, ),
                 vector("integer", 0L))

    expect_equal(rank_integer_index(c(TRUE, FALSE, TRUE), x_mat),
                 c(1L, 3L))
    expect_error(rank_integer_index(c(TRUE, FALSE), x_mat))

    expect_equal(rank_integer_index(c(2L, 3L), x_mat),
                 c(2L, 3L))

    expect_error(rank_integer_index(c("aa", "bb"), x_mat))

    colnames(x_mat) <- c("aa", "bb", "cc")
    expect_equal(rank_integer_index(c("aa", "bb"), x_mat),
                 c(1L, 2L))

    expect_error(rank_integer_index(c("aa", "zz"), x_mat))
})


test_that("testing ranked_x", {
    x_mat <- matrix(1:1000, ncol = 10)

    expect_error(ranked_x(x_mat, rank_cols = c("col1", "col3")))
    expect_error(ranked_x(x_mat, rank_cols = c(1, 3)),
                 NA)

    ## need all ten cols to use logical
    expect_error(ranked_x(x_mat,
                          rank_cols = c(TRUE, FALSE, TRUE, rep(FALSE, 6L))))
    expect_error(ranked_x(x_mat,
                          rank_cols = c(TRUE, FALSE, TRUE, rep(FALSE, 7L))),
                 NA)

    expect_equal(ranked_x(x_mat,
                          rank_cols = c(TRUE, FALSE, TRUE, rep(FALSE, 7L))),
                 ranked_x(x_mat, rank_cols = c(1, 3)))

    colnames(x_mat) <- paste0("col", 1L:10L)
    expect_error(ranked_x(x_mat, rank_cols = c("col1", "col3")),
                 NA)
    expect_equal(ranked_x(x_mat, rank_cols = c(1, 3)),
                 ranked_x(x_mat, rank_cols = c("col1", "col3")))

    ##------------------------------------

    x_mat <- cbind((1:10)^3,
                   rnorm(10))

    expect_false(isTRUE(all.equal(x_mat,
                                  ranked_x(x_mat, c(1, 2)))))
    expect_equal(rank(ranked_x(x_mat, c(1, 2))[, 1]),
                 1L:10L)
})


test_that("testing near_given_match", {
    sorted_vec <- sort(runif(1000L))
    given_ind <- fixed_sample(1L:1000L, 300L, replace = FALSE)
    given_value <- sorted_vec[given_ind]

    above_ind <- pmin(given_ind + 1L, 1000L)
    below_ind <- pmax(given_ind - 1L, 1L)
    above_value <- sorted_vec[above_ind]
    above_value[given_ind == 1000L] <- Inf
    below_value <- sorted_vec[below_ind]
    below_value[given_ind == 1L] <- -Inf

    closest_ind <- ifelse(above_value - given_value > given_value - below_value,
                          below_ind, above_ind)

    new_order <- fixed_sample(1L:1000L)
    new_value_vec <- sorted_vec[new_order]

    new_given <- match(given_ind, new_order)
    new_closest <- match(closest_ind, new_order)

    expect_equal(new_closest,
                 near_given_match(sorted_vec[new_order],
                                  new_given))

    ##------------------------------------

    random_vec <- runif(500L)
    random_given_lgl <- (1L:500L) %in% fixed_sample(1L:500L, 30L)
    random_given <- which(random_given_lgl)

    dist_mat <- abs(outer(random_vec[random_given_lgl],
                          random_vec,
                          "-"))
    dist_mat[cbind(seq_len(length(random_given)),
                   random_given)] <- Inf

    min_dist_ind <- apply(dist_mat, 1, function(x) {
        which(rank(x, ties.method = "random") == 1)
    })

    ## this is must faster because we just have one vector,
    ## here about 50 times so
    closest_ind <- near_given_match(random_vec, random_given)

    expect_equal(min_dist_ind, closest_ind)
})


test_that("testing expit", {
    x <- runif(10)
    expect_equal(expit(x),
                 (exp(x) / (1 + exp(x))))
})


test_that("testing sym_mat", {
    expect_error(sym_mat())

    expect_equal(sym_mat(runif(1, -10000000, 1000000)),
                 matrix(0, 1, 1))

    diag_mat <- diag(x = c(rnorm(50), runif(50)))
    expect_equal(sym_mat(diag_mat),
                 matrix(0, 100, 100))

    expect_error(sym_mat(matrix(runif(10), 5, 2)))

    some_matrix <- rbind(
        c(8, 2, 4),
        c(-3, 3, 3),
        c(0, 2, 2)
    )
    res_matrix <- rbind(
        c(0, 2 - 3, 4 + 0),
        c(-3 + 2, 0, 3 + 2),
        c(0 + 4, 2 + 3, 0)
    )
    expect_equal(sym_mat(some_matrix),
                 res_matrix)
})


test_that("testing binary_search", {
    root_3 <- binary_search(3,
                            monotone_function = function(x) {
                                x^2
                            },
                            init_bounds = c(0, 2))
    root_3_again <- binary_search(3,
                                  monotone_function = function(x) {
                                      x^2
                                  })

    expect_true(abs(root_3 - sqrt(3)) < 1e-5)
    expect_true(abs(root_3_again - sqrt(3)) < 1e-5)

    expect_error(binary_search(3,
                               monotone_function = function(x) {
                                   x^2
                               },
                               max_iters = 10L))

    ##------------------------------------

    mono_func <- function(x) {
        - x / 3 + 50
    }

    mono_sol <- binary_search(target_value = -600,
                              monotone_function = mono_func,
                              error_gap = 1e-7)

    expect_true(abs(mono_func(mono_sol) + 600) < 1e-6)
    expect_false(abs(mono_func(mono_sol) + 600) < 1e-9)

    mono_sol <- binary_search(target_value = 1200,
                              monotone_function = mono_func,
                              error_gap = 1e-7)

    expect_true(abs(mono_func(mono_sol) - 1200) < 1e-6)
    expect_false(abs(mono_func(mono_sol) - 1200) < 1e-9)
})


test_that("testing tolerance_check", {
    control_tol <- runif(1000L)
    treat_tol <- control_tol + runif(1000L, 1, 3)
    tol_vec <- c(control_tol, treat_tol)

    match_list <- list(
        treat_index = 1001L:2000L,
        control_index = 1L:1000L,
        distance = runif(1000L)
    )

    nothing_wrong <- tolerance_check(
        match_list,
        gen_tolerance_list(
            tolerance_vec = tol_vec
        )
    )
    expect_false(nothing_wrong[["error"]])

    ## ------------------------------------

    min_violations <- tolerance_check(
        match_list,
        gen_tolerance_list(
            tolerance_vec = tol_vec,
            tolerance_min = 2
        )
    )

    expect_true(min_violations[["error"]])
    expect_true(grepl("have difference below or at",
                      min_violations[["message"]]))

    wrong_dir <- tolerance_check(
        match_list,
        gen_tolerance_list(
            tolerance_vec = -tol_vec
        )
    )

    expect_true(wrong_dir[["error"]])
    expect_true(grepl("all pairs have treatment tolerance",
                      wrong_dir[["message"]]))
    expect_true(grepl("did the treatment and control get switched?",
                      wrong_dir[["message"]]))

    max_violations <- tolerance_check(
        match_list,
        gen_tolerance_list(
            tolerance_vec = tol_vec,
            tolerance_max = 2
        )
    )

    expect_true(max_violations[["error"]])
    expect_true(grepl("have difference above",
                      max_violations[["message"]]))

    random_tol <- tolerance_check(
        match_list,
        gen_tolerance_list(
            tolerance_vec = sample(tol_vec),
            tolerance_min = 2
        )
    )
    expect_true(random_tol[["error"]])
    expect_true(grepl("treatment tolerance with lower value than",
                      random_tol[["message"]]))
    expect_true(grepl("a further \\d{1,3} pairs have min constraint violated",
                      random_tol[["message"]]))
})
