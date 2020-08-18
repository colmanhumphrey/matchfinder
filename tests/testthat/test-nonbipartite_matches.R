context("testing nonbipartite_matches")

test_that("testing gen_tolerance_list", {
    expect_equal(NULL,
                 gen_tolerance_list())

    expect_equal(NULL,
                 gen_tolerance_list(tolerance_vec = NULL))

    expect_error(gen_tolerance_list(tolerance_min = 3))
    expect_error(gen_tolerance_list(tolerance_vec = c(1, 2, 10),
                                    tolerance_min = c(3, 10)))
    expect_error(gen_tolerance_list(tolerance_vec = c(1, 2, 10),
                                    tolerance_min = c(3, 10)))
    expect_error(gen_tolerance_list(tolerance_vec = c(1, 2, 10),
                                    tolerance_min = c(3),
                                    tolerance_max = c(20, 1)))

    tolerance_vec <- runif(10)
    tolerance_min <- 0.2
    tolerance_max <- 10

    expect_equal(
        list(tolerance_vec = tolerance_vec,
             tolerance_min = tolerance_min,
             tolerance_max = tolerance_max),
        gen_tolerance_list(tolerance_vec, tolerance_min, tolerance_max)
    )

    expect_equal(
        list(tolerance_vec = tolerance_vec,
             tolerance_min = tolerance_min,
             tolerance_max = NULL),
        gen_tolerance_list(tolerance_vec, tolerance_min)
    )

    expect_equal(
        list(tolerance_vec = tolerance_vec,
             tolerance_min = 0,
             tolerance_max = NULL),
        gen_tolerance_list(tolerance_vec)
    )

    ##------------------------------------
    ## testing tolerance_to_caliper_list
    tol_list <- gen_tolerance_list(tolerance_vec = c(1, 2, 10),
                                   tolerance_min = 2,
                                   tolerance_max = 3)
    cal_min_list <- tolerance_to_caliper_list(tol_list, use_min = TRUE)
    cal_max_list <- tolerance_to_caliper_list(tol_list, use_min = FALSE)

    expect_equal(
        list(
            caliper_vec = c(1, 2, 10),
            caliper_max = 2,
            continuous_mult = 1),
        cal_min_list)
    expect_equal(
        list(
            caliper_vec = c(1, 2, 10),
            caliper_max = 3,
            continuous_mult = 1),
        cal_max_list)
})


test_that("testing errors", {
    ## need dist_mat, must be symmetric
    expect_error(nonbipartite_matches())
    expect_error(nonbipartite_matches(dist_mat = matrix(1:4)))
    expect_error(nonbipartite_matches(
        dist_mat = matrix(c(1, 2, 2, 1), 2, 2),
        tolerance_list = gen_tolerance_list()),
        NA)

    ## dist_mat must be matrix
    expect_error(nonbipartite_matches(
        dist_mat = 1:100,
        tolerance_list = gen_tolerance_list(rep(c(0, 1), times = c(5, 5)))))
    expect_error(nonbipartite_matches(
        dist_mat = sym_mat(matrix(runif(100), 10, 10)),
        tolerance_list = gen_tolerance_list(rep(c(0, 1), times = c(5, 5)))),
        NA)
    expect_error(nonbipartite_matches(
        dist_mat = sym_mat(matrix(runif(100), 10, 10)),
        tolerance_list = gen_tolerance_list(rep(c(0, 1), times = c(5, 5))),
        n_sinks = 0),
        NA)

    ## tolerance_vec length must be = rows = cols of dist_mat,
    ## and must be logical or numeric
    expect_error(nonbipartite_matches(
        dist_mat = sym_mat(matrix(runif(100), 10, 10)),
        tolerance_list = gen_tolerance_list(rep(c(0, 1), times = c(4, 4)))))
    expect_error(
        nonbipartite_matches(
            dist_mat = sym_mat(matrix(runif(100), 10, 10)),
            tolerance_list = gen_tolerance_list(rep(c("0", "1"), times = c(5, 5)))))
    expect_error(
        nonbipartite_matches(
            dist_mat = sym_mat(matrix(runif(100), 10, 10)),
            tolerance_list = gen_tolerance_list(rep(c(0, 2), times = c(5, 5)))),
        NA)
    expect_error(
        nonbipartite_matches(
            dist_mat = sym_mat(matrix(runif(100), 10, 10)),
            tolerance_list = gen_tolerance_list(rep(c(FALSE, TRUE),
                                                    times = c(5, 5)))),
        NA)

    expect_error(nonbipartite_matches(
        dist_mat = sym_mat(matrix(runif(100), 10, 10)),
        tolerance_list = gen_tolerance_list(rep(c(FALSE, TRUE),
                                                times = c(5, 5))),
        match_method = "best_match"))
    expect_error(
        nonbipartite_matches(
            dist_mat = sym_mat(matrix(runif(100), 10, 10)),
            tolerance_list = gen_tolerance_list(rep(c(FALSE, TRUE),
                                                    times = c(5, 5))),
            match_method = c("optimal", "greedy")))
    expect_error(nonbipartite_matches(
        dist_mat = sym_mat(matrix(runif(100), 10, 10)),
        tolerance_list = gen_tolerance_list(rep(c(FALSE, TRUE),
                                                times = c(5, 5))),
        match_method = "greedy"),
        NA)
})


test_that("testing nonbipartite_matches", {
    x_mat <- cbind(
        c(0.1, 2.4, 0.1, -0.9, 2.3, 1.8, -0.9, 0.5, 2.9, -0.8,
          2.5, 2, -1, -1.4, -2.3, 2.1, 2.4, 1.4, -0.5, -0.8),
        c(3.3, 3.7, 1.5, 1.7, 3.9, 3.7, 3.2, 2.7, 3.4, 2.3, 3.7,
          2.9, 1.4, 1.6, 2, 2.7, 2.6, 2.3, 3.3, 3.2)
    )
    dist_mat <- weighted_mahal(ranked_x(x_mat, 2),
                               covariance_with_ranks(x_mat, 2))
    diag(dist_mat) <- Inf

    tolerance_vec <- rep((1L:10L) * 3L, each = 2L)

    wr_res <- nonbipartite_matches(
        dist_mat,
        tolerance_list = gen_tolerance_list(tolerance_vec,
                                            tolerance_min = 2),
        match_method = "with_replacement",
        n_sinks = c(0L, 1L, 2L))
    greedy_res <- nonbipartite_matches(
        dist_mat,
        tolerance_list = gen_tolerance_list(tolerance_vec,
                                            tolerance_min = 2),
        match_method = "greedy",
        n_sinks = c(0L, 1L, 2L)
    )
    opt_res <- nonbipartite_matches(
        dist_mat,
        tolerance_list = gen_tolerance_list(tolerance_vec,
                                            tolerance_min = 2),
        match_method = "optimal",
        n_sinks = c(0L, 2L, 4L)
    )

    all_res_list <- list(
        "with_replacement" = wr_res,
        "greedy" = greedy_res,
        "optimal" = opt_res
    )

    is_ordered <- unlist(lapply(all_res_list, function(x) {
        ranks <- rank(unlist(lapply(x, function(y) {
            length(y[["treat_index"]])
        })), ties.method = "first") # first to penalise accidental ties
        isTRUE(all.equal(unname(ranks), c(3L, 2L, 1L)))
    }))
    expect_true(all(is_ordered))

    ##----------------

    zero_list <- lapply(all_res_list, function(x) {
        x[["0"]]
    })
    ## should get all matched up
    zero_len <- unlist(lapply(zero_list, function(x) {
        length(x[["treat_index"]])
    }))
    expect_true(all(zero_len == 10L))

    ## order of dists should be wr < optimal < greedy
    zero_dists <- unlist(lapply(zero_list, function(x) {
        sum(x[["distance"]])
    }))
    expect_true(zero_dists["with_replacement"] < zero_dists["optimal"])
    expect_true(zero_dists["optimal"] < zero_dists["greedy"])

    ##----------------

    eight_list <- lapply(all_res_list, function(x) {
        x[[3]]
    })
    eight_len <- unlist(lapply(eight_list, function(x) {
        length(x[["treat_index"]])
    }))
    expect_true(all(eight_len == 8L))

    ## order of dists should be wr < optimal < greedy
    eight_dists <- unlist(lapply(eight_list, function(x) {
        sum(x[["distance"]])
    }))
    expect_true(eight_dists["with_replacement"] < eight_dists["optimal"])
    expect_true(eight_dists["optimal"] < eight_dists["greedy"])
})


test_that("testing with_replacement_nbp_match", {
    x_mat <- cbind(
        c(4.1, 49.9, 34.3, 22.1, 3.1, 7, 6, 1, 28, 1, 7, 10,
          27, 10, 25, 40, 2, 6, 47, 4),
        c(31, 4, 3, 34, 1, 29, 2, 7, 5.1, 3.1, 54, 54,
          15, 6, 11, 15, 7, 16.7, 22.9, 20.3))

    dist_mat <- weighted_mahal(x_mat,
                               cov(x_mat))
    diag(dist_mat) <- Inf
    tolerance_vec <- runif(nrow(x_mat))

    wr_res <- with_replacement_nbp_match(dist_mat,
                                         tolerance_vec = tolerance_vec,
                                         keep_all = FALSE)
    greedy_res <- greedy_nbp_match(dist_mat,
                                   tolerance_vec = tolerance_vec)
    expect_true(sum(wr_res[["distance"]]) <
                sum(greedy_res[["distance"]]))

    treated_units <- wr_res[["treat_index"]]
    control_units <- wr_res[["control_index"]]

    expect_true(!any(duplicated(treated_units)))
    expect_true(all(treated_units != control_units))
    expect_true(all(tolerance_vec[treated_units] >
                    tolerance_vec[control_units]))

    ## tolerance of course blocks some matches:
    min_vals <- apply(dist_mat, 1, min)
    expect_true(sum(min_vals[rank(min_vals, ties.method = "first") <= 10L]) <
                sum(wr_res[["distance"]]))

    ##------------------------------------

    wr_res_all <- with_replacement_nbp_match(dist_mat,
                                             tolerance_vec = tolerance_vec,
                                             keep_all = TRUE)
    expect_equal(length(wr_res_all[["treat_index"]]),
                 nrow(dist_mat))
    treated_units <- wr_res_all[["treat_index"]]
    control_units <- wr_res_all[["control_index"]]

    ## at least one unit will have no real match
    inf_ind <- wr_res_all[["distance"]] == Inf
    expect_true(any(inf_ind))

    expect_true(!any(duplicated(treated_units[!inf_ind])))
    expect_true(all(treated_units[!inf_ind] != control_units[!inf_ind]))
    expect_true(all(tolerance_vec[treated_units[!inf_ind]] >
                    tolerance_vec[control_units[!inf_ind]]))
})


test_that("testing greedy_nbp_match", {
    x_mat <- matrix(rnorm(100L), 25L, 4L)
    weight_vec <- c(1, 3, 3, 2)
    dist_mat <- weighted_mahal(x_mat,
                               cov(x_mat),
                               weight_vec = weight_vec)
    diag(dist_mat) <- Inf
    tolerance_vec <- runif(25L)

    greedy_res <- greedy_nbp_match(dist_mat,
                                   tolerance_vec)
    treated_units <- greedy_res[["treat_index"]]
    control_units <- greedy_res[["control_index"]]

    expect_true(!any(duplicated(c(treated_units, control_units))))
    expect_true(all(c(treated_units, control_units) %in% 1L:25L))

    expect_true(all(tolerance_vec[treated_units] >
                    tolerance_vec[control_units]))
})


test_that("testing optimal_nbp_match", {
    x_mat <- cbind(c(0, 0, 0, 1, 1, 1),
                   c(0, 3, 2, 1, 0, 2))
    tolerance_vec <- c(0, 5, 0, 10, 20, 10)
    dist_mat <- weighted_mahal(x_mat,
                               covariance_with_ranks(x_mat),
                               treat_vec = NULL)
    ## simulate tol mins
    dist_mat[cbind(c(1, 3, 4, 6),
                   c(3, 1, 6, 4))] <- Inf
    diag(dist_mat) <- Inf
    nbp_result <- optimal_nbp_match(dist_mat, tolerance_vec = tolerance_vec)

    expect_equal(sort(nbp_result[["treat_index"]]),
                 c(2, 5, 6))
    expect_equal(sort(nbp_result[["control_index"]]),
                 c(1, 3, 4))
    expect_true(abs(sum(sort(nbp_result[["distance"]]) -
                        c(0.5, 0.5, 1.5))) < 0.5)

    ##----------------
    ## greedy res is somewhat random
    opt_dist <- sum(nbp_result[["distance"]])

    greedy_res_list <- lapply(1:50, function(j) {
        greedy_nbp_match(dist_mat, tolerance_vec = tolerance_vec)
    })
    dist_vec <- unlist(lapply(greedy_res_list, function(x) {
        sum(x[["distance"]])
    }))
    num_pairs <- unlist(lapply(greedy_res_list, function(x) {
        length(x[["distance"]])
    }))

    ## we will likely sometimes only get two pairs, because
    ## if we are left with the 1,3 or 4,6 pairs, they can't be matched
    expect_false(all(num_pairs == 3L))

    ## and we'll sometimes do sub-optimally
    dist_with_three <- dist_vec[num_pairs == 3L]
    expect_true(any(dist_with_three > opt_dist))

    ##------------------------------------
    ## with sinks
    ## this shows that including sinks is better than
    ## later throwing away the worst pairs

    x_mat <- cbind(
        c(4, 49, 34, 22, 3, 7, 6, 1, 28, 1, 7, 10,
          27, 10, 25, 40, 2, 6, 47, 4),
        c(31, 4, 3, 34, 1, 29, 2, 7, 5, 3, 54, 54,
          15, 6, 11, 15, 7, 16, 22, 20))
    dist_mat <- weighted_mahal(ranked_x(x_mat, 2),
                               covariance_with_ranks(x_mat, 2))
    diag(dist_mat) <- Inf

    opt_res <- optimal_nbp_match(dist_mat, n_sinks = 4)
    opt_full <- optimal_nbp_match(dist_mat, n_sinks = 0)

    ## 4 sinks here means we lose two pairs
    expect_equal(length(opt_res[["treat_index"]]),
                 8L)
    expect_equal(length(opt_full[["treat_index"]]),
                 10L)

    full_dist <- opt_full[["distance"]]
    keep_8 <- rank(full_dist, ties.method = "first") <= 8L
    opt_8_from_full <- lapply(opt_full, function(x) {
        x[keep_8]
    })

    expect_true(sum(opt_8_from_full[["distance"]]) >
                sum(opt_res[["distance"]]))
})


test_that("testing add_nbp_sinks", {
    if (!requireNamespace("nbpMatching")) {
        skip("nbpMatching not installed, skipping some nbp tests")
    }

    ##----------------
    ## even case:
    sym_matrix <- sym_mat(matrix(rnorm(36L), 6L, 6L))

    distance_add_0 <- add_nbp_sinks(sym_matrix, 0L)
    expect_true(sum(abs(distance_add_0 - sym_matrix)) == 0)

    expect_warning(add_nbp_sinks(sym_matrix, 1L))
    distance_add_1 <- suppressWarnings(add_nbp_sinks(sym_matrix, 1L))
    expect_equal(sum(distance_add_1[distance_add_1 != Inf]), sum(sym_matrix))

    distance_add_4 <- add_nbp_sinks(sym_matrix, 4L)
    expect_equal(dim(distance_add_4), dim(sym_matrix) + c(4L, 4L))

    ##----------------
    ## odd case:
    sym_matrix <- sym_mat(matrix(rnorm(49L), 7L, 7L))

    expect_warning(add_nbp_sinks(sym_matrix, 0L))
    distance_add_0 <- suppressWarnings(add_nbp_sinks(sym_matrix, 0L))
    expect_equal(sum(distance_add_0), sum(sym_matrix))
    expect_equal(dim(distance_add_0), dim(sym_matrix) + c(1L, 1L))


    distance_add_1 <- add_nbp_sinks(sym_matrix, 1L)
    expect_equal(unname(distance_add_1), unname(distance_add_0))

    distance_add_5 <- add_nbp_sinks(sym_matrix, 5L)
    expect_equal(dim(distance_add_5), dim(sym_matrix) + c(5L, 5L))
})


test_that("testing fix_nbp_match", {
    if (!requireNamespace("nbpMatching")) {
        skip("nbpMatching not installed, skipping some nbp tests")
    }

    ##----------------

    nbp_opt_input <- add_nbp_sinks(sym_mat(matrix(runif(16), 4, 4)),
                                   0)

    nbp_opt_output <- nbpMatching::nonbimatch(nbp_opt_input)
    fixed_output <- fix_nbp_match(nbp_opt_output, nrow_match = 4L)

    expect_equal(names(fixed_output),
                 c("treat_index", "control_index", "distance"))

    ##----------------

    sym_matrix <- sym_mat(matrix(runif(16), 4, 4))
    sym_matrix[cbind(c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4),
                  c(3, 4, 2, 3, 4, 1, 1, 2, 1, 2))] <- Inf


    nbp_opt_input <- add_nbp_sinks(sym_matrix, 0)

    nbp_opt_output <- nbpMatching::nonbimatch(nbp_opt_input)
    fixed_output <- fix_nbp_match(nbp_opt_output, nrow_match = 4L,
                                  tolerance_vec = c(1, 1, 10, 9))

    expect_equal(fixed_output[["treat_index"]], 3L)
    expect_equal(fixed_output[["control_index"]], 4L)

    ##----------------

    nbp_opt_input <- add_nbp_sinks(sym_mat(matrix(runif(16), 4, 4)),
                                   n_sinks = 2)

    nbp_opt_output <- nbpMatching::nonbimatch(nbp_opt_input)
    fixed_output <- fix_nbp_match(nbp_opt_output, nrow_match = 4L)

    expect_true(length(fixed_output[["treat_index"]]) == 1L)
})


test_that("testing reorder_nbp", {
    treat_index <- fixed_sample(1L:100L, 40L)
    control_index <- fixed_sample((1L:100L)[!((1L:100L) %in% treat_index)],
                                  40L)
    match_list <- list(
        treat_index = treat_index,
        control_index = control_index,
        distance = runif(40L)
    )

    reordered <- reorder_nbp(match_list)
    expect_equal(setdiff(reordered[["treat_index"]], treat_index),
                 integer(0))
    expect_equal(setdiff(reordered[["control_index"]], control_index),
                 integer(0))
    expect_equal(setdiff(reordered[["distance"]], match_list[["distance"]]),
                 numeric(0))

    expect_equal(rank(reordered[["treat_index"]]),
                 seq_len(length(treat_index)))

    expect_equal(reordered[["control_index"]],
                 control_index[order(treat_index)])

    ##------------------------------------

    treat_index <- fixed_sample(1L:100L, 40L)
    control_index <- fixed_sample((1L:100L)[!((1L:100L) %in% treat_index)],
                                  40L)
    match_list <- list(
        treat_index = treat_index,
        control_index = control_index,
        distance = runif(40L)
    )

    tolerance_vec <- 1L:100L

    reordered_both <- reorder_nbp(match_list,
                                  tolerance_vec)

    expect_equal(setdiff(c(match_list[["treat_index"]],
                          match_list[["control_index"]]),
                        c(reordered_both[["treat_index"]],
                          reordered_both[["control_index"]])),
                integer(0))
    expect_equal(setdiff(reordered_both[["distance"]],
                         match_list[["distance"]]),
                 numeric(0))

    ## because tolerance_vec is increasing
    expect_true(all(reordered_both[["treat_index"]] >
                    reordered_both[["control_index"]]))
})
