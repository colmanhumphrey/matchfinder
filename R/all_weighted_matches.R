#' idea is to complete a search for a given N,
#' or vector N
#'
#' @param x_mat input/design matrix (already rank-adjusted etc)
#' @param cov_x the (potentially rank-adjusted) covariance of \code{x_mat}.
#'   This means it's possible that \code{cov(x_mat)} is not equal to
#'   \code{cov_x}; see \code{covariance_with_ranks} for more details.
#' @param weight_list list of weight vectors. See `generate_random_weights` to
#'   automatically generate a reasonable set of vectors.
#' @param treat_vec Logical (or 1/0) vector, indicating treatment (or control).
#' @param n_sinks how many potential matches to not bother with
#'   NOTE: you can do this as a vector, but not for optimal matching.
#' @param caliper_list Optional, see \code{gen_caliper_list}. Provide
#'   this to force matches that are close on some metric.
#' @param tol_val For optimal matches, you can set a tolerance
#'   to be within optimality of, which can be zero for perfect optimality.
#'   Default 1e-4 is reasonable in many cases.
#'
#' @export
all_bipartite_matches <- function(x_mat,
                                  cov_x,
                                  weight_list,
                                  treat_vec,
                                  match_method = c(
                                      "with_replacement",
                                      "optimal",
                                      "greedy"
                                  ),
                                  n_sinks = 0L,
                                  caliper_list = gen_caliper_list(),
                                  propensity_list = match_propensity_list(NULL),
                                  sqrt_mahal = TRUE,
                                  tol_val = NULL) {
    if (is.null(n_sinks)) {
        n_sinks <- 0L
    }

    if (!is.null(propensity_list)) {
        if (!is.null(caliper_list)) {
            stop(
                "don't use both `caliper_list` and `propensity_list`: ",
                " If you do want both, create the combined caliper separately"
            )
        }

        ## in case of logical
        treat_vec <- treat_vec * 1L

        ## generate propensity score
        prop_list_names <- c(
            "propensity_function",
            "oos_propensity",
            "n_folds"
        )
        prop_score <- propensity_score(
            x_mat = x_mat,
            treat_vec = treat_vec,
            propensity_list = propensity_list[prop_list_names]
        )
        caliper_list <- gen_caliper_list(
            caliper_vec = prop_score,
            caliper_max = sd(prop_score) * propensity_list[["caliper_sd_mult"]],
            continuous_mult = propensity_list[["continuous_mult"]]
        )
    }

    if (!is.null(caliper_list)) {
        caliper_dist_mat <- create_caliper(caliper_list,
            treat_vec = treat_vec
        )
    }

    by_weight_list <- lapply(weight_list, function(weight_vec) {
        w_dist_mat <- weighted_mahal(x_mat,
                                     cov_x = cov_x,
                                     weight_vec = weight_vec,
                                     treat_vec = treat_vec,
                                     sqrt_mahal = sqrt_mahal
                                     )

        if (!is.null(caliper_list)) {
            w_dist_mat <- w_dist_mat + caliper_dist_mat
        }

        bipartite_matches(
            dist_mat = w_dist_mat,
            treat_vec = treat_vec,
            match_method = match_method,
            n_sinks = n_sinks,
            tol_val = tol_val,
            weight_vec = weight_vec
        )
    })

    ## more natural to group by sink value
    setNames(lapply(n_sinks, function(x) {
        lapply(by_weight_list, function(y) {
            y[[as.character(x)]]
        })
    }), n_sinks)
}


#' Computes all matches, then gets the brier scores for each. Reorder by
#' number of sinks.
#'
#' @inheritParams all_bipartite_matches
#' @inheritParams brier_score_cv
#' @return List of matches within sink values,
#'  and brier scores for each.
#' @author Colman Humphrey
#'
#' @export
brier_bipartite_matches <- function(x_mat,
                                    cov_x,
                                    weight_list,
                                    treat_vec,
                                    match_method = c(
                                        "with_replacement",
                                        "optimal",
                                        "greedy"
                                    ),
                                    n_sinks = 0L,
                                    caliper_list = gen_caliper_list(),
                                    propensity_list =
                                        match_propensity_list(NULL),
                                    sqrt_mahal = TRUE,
                                    tol_val = NULL,
                                    design = "cross_all",
                                    num_folds = 5,
                                    match_predict_function =
                                        match_predict_xgb(),
                                    silent = !interactive()) {
    if (is.null(n_sinks)) {
        n_sinks <- 0L
    }

    ## generate all matches: one per weight vector per n_sink value
    all_matches <- all_bipartite_matches(
        x_mat = x_mat,
        cov_x = cov_x,
        weight_list = weight_list,
        treat_vec = treat_vec,
        match_method = match_method,
        n_sinks = n_sinks,
        caliper_list = caliper_list,
        propensity_list = propensity_list,
        sqrt_mahal = sqrt_mahal,
        tol_val = tol_val
    )

    if (!silent) {
        message("getting briers")
    }

    ## get all brier scores for all results
    briers_by_sinks <- lapply(all_matches, function(all_by_sink) {
        if (!silent) {
            print(all_by_sink[[1]]["num_sinks"])
        }
        unlist(lapply(all_by_sink, function(indiv_match_list) {
            brier_score_cv(
                x_mat = x_mat,
                match_list = indiv_match_list,
                design = design,
                num_folds = num_folds,
                match_predict_function = match_predict_function
            )
        }))
    })

    list(
        matches_by_sinks = all_matches,
        briers_by_sinks = briers_by_sinks
    )
}


#' Takes matches and their Brier scores, and computes
#' permutation Brier scores and the best matches
#'
#' Can work for bipartite or non-bipartite matches,
#' the permutation is just over the labels.
#' @param matches_by_sinks List by number of sinks, each a list
#'   of match results (a match list), for each weight vector.
#' @param briers_by_sinks List by number of sinks, each a vector
#'   of Brier results. Basically a number for each match in
#'   \code{matches_by_sinks}.
#' @param x_mat Typical input matrix
#' @param n_sinks Vector of number of sinks - probably could get this
#'   directly from \code{matches_by_sinks}, but nice to be explicit.
#' @param approximate_by_best Logical, default \code{TRUE}. Only
#'   compute one permutation distribution, using the best result by brier
#'   score to do so. Useful because it changes little, but saves a ton of
#'   time.
#' @param silent Do you want to suppress message output? Default
#'   \code{!interactive()}.
#' @return Returns a list of two lists. The first is vectors
#'   of permutation Brier scores (one per match). The second is the
#'   best match at each sink value, along with some extra info about that
#'   match.
#' @author Colman Humphrey
#'
#' @export
permutation_matches <- function(matches_by_sinks,
                                briers_by_sinks,
                                x_mat,
                                n_sinks = 0L,
                                approximate_by_best = TRUE,
                                silent = !interactive()) {
    if (is.null(n_sinks)) {
        n_sinks <- 0L
    }

    if (length(n_sinks) != length(matches_by_sinks)) {
        stop("`n_sinks` is length ", length(n_sinks),
             " but `matches_by_sinks` is length ", length(matches_by_sinks))
    }

    if (length(n_sinks) != length(briers_by_sinks)) {
        stop("`n_sinks` is length ", length(n_sinks),
             " but `briers_by_sinks` is length ", length(briers_by_sinks))
    }

    best_brier_inds <- lapply(briers_by_sinks, function(x) {
        which(rank(x, ties.method = "first") == length(x))
    })

    if (!silent) {
        message("running permutations, will be a little slow")
    }

    permutation_briers <- lapply(seq_len(length(n_sinks)), function(j) {
        if (!silent) {
            message(paste0("Running permutations for ",
                           n_sinks[j], " sinks"))
        }
        if (approximate_by_best) {
            ## we'll just use the best to save time
            best_brier_ind <- best_brier_inds[[j]]
            return(permutation_brier(
                x_mat,
                match_list = matches_by_sinks[[j]][[best_brier_ind]]
            ))
        } else {
            return(lapply(matches_by_sinks[[j]], function(match_list) {
                permutation_brier(
                    x_mat,
                    match_list = match_list
                )
            }))
        }
    })

    ## compute the permutation score for each match
    permutation_brier_scores <- setNames(lapply(seq_len(length(n_sinks)),
                                                function(j) {
        if (approximate_by_best) {
            permutation_vec <- permutation_briers[[j]]
            return(unlist(lapply(briers_by_sinks[[j]], function(x) {
                mean(x <= permutation_vec)
            })))
        } else {
            return(unlist(lapply(
                seq_len(length(briers_by_sinks[[j]])), function(k) {
                    permutation_vec <- permutation_briers[[j]][[k]]
                    mean(briers_by_sinks[[j]][[k]] <= permutation_vec)
            })))
        }
    }), n_sinks)

    ## now that we're doing one-sided brier,
    ## the lowest value will just be the best
    ## so will highest brier

    best_matches <- setNames(lapply(seq_len(length(n_sinks)), function(j) {
        best_brier_ind <- if (approximate_by_best) {
                              best_brier_inds[[j]]
                          } else {
                              which.min(permutation_brier_scores[[j]])
                          }

        if (approximate_by_best) {
            ## because they all used the same vector -
            ## so the best must also have the best score here
            stopifnot(
                permutation_brier_scores[[j]][best_brier_ind] ==
                min(permutation_brier_scores[[j]])
            )
        }

        list(
            n_sinks = n_sinks[j],
            raw_brier = briers_by_sinks[[j]][best_brier_ind],
            permutation_brier = permutation_brier_scores[[j]][best_brier_ind],
            match_list = matches_by_sinks[[j]][[best_brier_ind]]
        )
    }), n_sinks)

    list(
        permutation_brier_scores = permutation_brier_scores,
        best_matches = best_matches
    )
}


#' Same as \code{all_bipartite_matches} but for nbp matching
#'
#' @inheritParams all_bipartite_matches
#' @inheritParams nonbipartite_matches
#' @param tolerance_list See \code{gen_tolerance_list}
#'
#' @export
all_nonbipartite_matches <- function(x_mat,
                                     cov_x,
                                     weight_list,
                                     tolerance_list = gen_tolerance_list(),
                                     match_method = c(
                                         "with_replacement",
                                         "optimal",
                                         "greedy"
                                     ),
                                     n_sinks = 0L,
                                     caliper_list = gen_caliper_list(),
                                     propensity_list =
                                         match_propensity_list(NULL),
                                     sqrt_mahal = TRUE,
                                     keep_all_with_replacement = FALSE) {
    if (is.null(n_sinks)) {
        n_sinks <- 0L
    }

    if (!is.null(propensity_list)) {
        if (!is.null(caliper_list)) {
            stop(
                "don't use both `caliper_list` and `propensity_list`: ",
                " If you do want both, create the combined caliper separately"
            )
        }

        ## generate propensity score
        prop_list_names <- c(
            "propensity_function",
            "oos_propensity",
            "n_folds"
        )
        prop_score <- propensity_score(
            x_mat = x_mat,
            treat_vec = tolerance_list[["tolerance_vec"]],
            propensity_list = propensity_list[prop_list_names]
        )
        caliper_list <- gen_caliper_list(
            caliper_vec = prop_score,
            caliper_max = sd(prop_score) * propensity_list[["caliper_sd_mult"]],
            continuous_mult = propensity_list[["continuous_mult"]]
        )
    }

    if (!is.null(caliper_list)) {
        caliper_dist_mat <- create_caliper(caliper_list)
    }

    if (missing(tolerance_list)) {
        warning("assuming that all pairs are matchable, ",
                "not neccesarily correct. If this is what you want, ",
                "you can silence this warning by explicitly supplying ",
                "`tolerance_vec = NULL`",
                call. = FALSE
                )
    }

    by_weight_list <- lapply(weight_list, function(weight_vec) {
        w_dist_mat <- weighted_mahal(x_mat,
            cov_x = cov_x,
            weight_vec = weight_vec,
            sqrt_mahal = sqrt_mahal
        )

        if (!is.null(caliper_list)) {
            w_dist_mat <- w_dist_mat + caliper_dist_mat
        }

        nonbipartite_matches(
            dist_mat = w_dist_mat,
            tolerance_list = tolerance_list,
            match_method = match_method,
            n_sinks = n_sinks,
            keep_all_with_replacement = keep_all_with_replacement,
            weight_vec = weight_vec
        )
    })

    ## more natural to group by sink value
    setNames(lapply(n_sinks, function(x) {
        lapply(by_weight_list, function(y) {
            y[[as.character(x)]]
        })
    }), n_sinks)
}


#' Computes all matches, then gets the brier scores for each. Reorder by
#' number of sinks.
#'
#' @inheritParams all_nonbipartite_matches
#' @inheritParams brier_score_cv
#' @return List of matches within sink values,
#'  and brier scores for each.
#' @author Colman Humphrey
#'
#' @export
brier_nonbipartite_matches <- function(x_mat,
                                       cov_x,
                                       weight_list,
                                       tolerance_list = gen_tolerance_list(),
                                       match_method = c(
                                           "with_replacement",
                                           "optimal",
                                           "greedy"
                                       ),
                                       n_sinks = 0L,
                                       caliper_list = gen_caliper_list(),
                                       propensity_list =
                                           match_propensity_list(NULL),
                                       sqrt_mahal = TRUE,
                                       keep_all_with_replacement = FALSE,
                                       design = "cross_all",
                                       num_folds = 5,
                                       match_predict_function =
                                           match_predict_xgb(),
                                       silent = !interactive()) {
    if (is.null(n_sinks)) {
        n_sinks <- 0L
    }

    ## generate all matches: one per weight vector per n_sink value
    all_matches <- all_nonbipartite_matches(
        x_mat = x_mat,
        cov_x = cov_x,
        weight_list = weight_list,
        tolerance_list = tolerance_list,
        match_method = match_method,
        n_sinks = n_sinks,
        caliper_list = caliper_list,
        propensity_list = propensity_list,
        sqrt_mahal = sqrt_mahal,
        keep_all_with_replacement = keep_all_with_replacement
    )

    if (!silent) {
        message("getting briers")
    }

    ## get all brier scores for all results
    briers_by_sinks <- lapply(all_matches, function(all_by_sink) {
        if (!silent) {
            print(all_by_sink[[1]]["num_sinks"])
        }
        unlist(lapply(all_by_sink, function(indiv_match_list) {
            brier_score_cv(
                x_mat = x_mat,
                match_list = indiv_match_list,
                design = design,
                num_folds = num_folds,
                match_predict_function = match_predict_function
            )
        }))
    })

    list(
        matches_by_sinks = all_matches,
        briers_by_sinks = briers_by_sinks
    )
}
