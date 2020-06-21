#' idea is to complete a search for a given N,
#' or vector N
#'
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
                                  match_method = c("with_replacement",
                                                   "optimal",
                                                   "greedy"),
                                  n_sinks = 0,
                                  caliper_list = gen_caliper_list(),
                                  propensity_list = match_propensity_list(NULL),
                                  sqrt_mahal = TRUE,
                                  tol_val = NULL) {
    if (!is.null(propensity_list)) {
        if (!is.null(caliper_list)) {
            stop("don't use both `caliper_list` and `propensity_list`: ",
                 " If you do want both, create the combined caliper separately")
        }

        ## in case of logical
        treat_vec <- treat_vec * 1L

        ## generate propensity score
        prop_list_names <- c("propensity_function",
                             "oos_propensity",
                             "n_folds")
        prop_score <- propensity_score(
            x_mat = x_mat,
            treat_vec = treat_vec,
            propensity_list = propensity_list[prop_list_names])
        caliper_list <- gen_caliper_list(
            caliper_vec = prop_score,
            caliper_max = sd(prop_score) * propensity_list[["caliper_sd_mult"]],
            continuous_mult = propensity_list[["continuous_mult"]])
    }

    if (!is.null(caliper_list)) {
        caliper_dist_mat <- create_caliper(caliper_list,
                                           treat_vec = treat_vec)
    }

    lapply(weight_list, function(weight_vec) {
        w_dist_mat <- weighted_mahal(x_mat,
                                     cov_x = cov_x,
                                     weight_vec = weight_vec,
                                     treat_vec = treat_vec,
                                     sqrt_mahal = sqrt_mahal)

        if(!is.null(caliper_list)){
            w_dist_mat <- w_dist_mat + caliper_dist_mat
        }

        bipartite_matches(dist_mat = w_dist_mat,
                          treat_vec = treat_vec,
                          match_method = match_method,
                          n_sinks = n_sinks,
                          tol_val = tol_val)
    })
}


#' Computes all matches, then gets the brier scores for each. Reorder by
#' number of sinks.
#'
#' @inheritParams all_bipartite_matches
#' @return List of matches within sink values,
#'  and brier scores for each.
#' @author Colman Humphrey
#'
#' @export
sink_brier_bipartite_matches <- function(x_mat,
                                         cov_x,
                                         weight_list,
                                         treat_vec,
                                         match_method = c("with_replacement",
                                                          "optimal",
                                                          "greedy"),
                                         n_sinks = 0,
                                         caliper_list = gen_caliper_list(),
                                         propensity_list = match_propensity_list(NULL),
                                         sqrt_mahal = TRUE,
                                         silent = !interactive(),
                                         tol_val = NULL) {

    ## generate all matches: one per weight vector per n_sink value
    all_matches <- all_bipartite_matches(
        x_mat = x_mat,
        cov_x = covariance_with_ranks(x_mat),
        weight_list = weight_list,
        treat_vec = treat_vec,
        match_method = match_method,
        n_sinks = n_sinks)

    ## reorder: by sink instead of by weight vector
    all_by_sinks <- setNames(lapply(n_sinks, function(x) {
        lapply(all_matches, function(y) {
            y[[as.character(x)]]
        })
    }), n_sinks)

    if (!silent) {
        message("getting briers")
    }

    ## get all brier scores for all results
    briers_by_sinks <- lapply(all_by_sinks, function(x){
        if (!silent) {
            print(x[[1]]["num_sinks"])
        }
        unlist(lapply(x, function(y) {
            brier_score_cv(x_mat,
                           y)
        }))
    })

    list(
        matches_by_sinks = all_by_sinks,
        briers_by_sinks = briers_by_sinks
    )
}


#' Takes matches and their Brier scores, and computes
#' permutation Brier scores and the best matches
#'
#' @param matches_by_sinks List by number of sinks, each a list
#'   of match results (a match list), for each weight vector.
#' @param briers_by_sinks List by number of sinks, each a vector
#'   of Brier results. Basically a number for each match in
#'   \code{matches_by_sinks}.
#' @param x_mat Typical input matrix
#' @param n_sinks Vector of number of sinks - probably could get this
#'   directly from \code{matches_by_sinks}, but nice to be explicit.
#' @param silent Do you want to suppress message output? Default
#'   \code{!interactive()}.
#' @return Returns a list of two lists. The first is vectors
#'   of permutation Brier scores (one per match). The second is the
#'   best match at each sink value, along with some extra info about that
#'   match.
#' @author Colman Humphrey
#'
#' @export
permutation_bipartite_matches <- function(matches_by_sinks,
                                          briers_by_sinks,
                                          x_mat,
                                          n_sinks,
                                          silent = !interactive()) {
    best_brier_inds <- lapply(briers_by_sinks, function(x) {
        which(rank(x, ties.method = "first") == length(x))
    })

    if (!silent) {
        message("running permutations, will be a little slow")
    }

    permutation_briers <- lapply(1L:length(n_sinks), function(j) {
        if (!silent) {
            n_sinks[j]
        }
        best_brier_ind <- best_brier_inds[[j]]
        permutation_brier(x_mat,
                          match_list = matches_by_sinks[[j]][[best_brier_ind]])
    })

    ## compute the permutation score for each match
    permutation_brier_scores <- lapply(1L:length(n_sinks), function(j) {
        permutation_vec <- permutation_briers[[j]]
        unlist(lapply(briers_by_sinks[[j]], function(x) {
            mean(x <= permutation_vec)
        }))
    })

    ## now that we're doing one-sided brier,
    ## the lowest value will just be the best
    ## so will highest brier

    best_matches <- lapply(1L:length(n_sinks), function(j) {
        best_brier_ind <- best_brier_inds[[j]]

        stopifnot(permutation_brier_scores[[j]][best_brier_ind] ==
                  min(permutation_brier_scores[[j]]))

        list(
            n_sinks = n_sinks[j],
            raw_brier = briers_by_sinks[[j]][best_brier_ind],
            permutation_brier = permutation_brier_scores[[j]][best_brier_ind],
            match_list = matches_by_sinks[[j]][[best_brier_ind]]
        )
    })

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
#' @export
all_nonbipartite_matches <- function(x_mat,
                                     cov_x,
                                     weight_list,
                                     tolerance_list = gen_tolerance_list(),
                                     match_method = c("with_replacement",
                                                      "optimal",
                                                      "greedy"),
                                     n_sinks = 0,
                                     caliper_list = gen_caliper_list(),
                                     propensity_list =
                                         match_propensity_list(NULL),
                                     sqrt_mahal = TRUE,
                                     keep_all_with_replacement = FALSE){
    if (!is.null(propensity_list)) {
        if (!is.null(caliper_list)) {
            stop("don't use both `caliper_list` and `propensity_list`: ",
                 " If you do want both, create the combined caliper separately")
        }

        ## generate propensity score
        prop_list_names <- c("propensity_function",
                             "oos_propensity",
                             "n_folds")
        prop_score <- propensity_score(
            x_mat = x_mat,
            treat_vec = tolerance_list[["tolerance_vec"]],
            propensity_list = propensity_list[prop_list_names])
        caliper_list = gen_caliper_list(
            caliper_vec = prop_score,
            caliper_max = sd(prop_score) * propensity_list[["caliper_sd_mult"]],
            continuous_mult = propensity_list[["continuous_mult"]])
    }

    if (!is.null(caliper_list)) {
        caliper_dist_mat <- create_caliper(caliper_list)
    }

    lapply(weight_list, function(weight_vec) {
        w_dist_mat <- weighted_mahal(x_mat,
                                     cov_x = cov_x,
                                     weight_vec = weight_vec,
                                     sqrt_mahal = sqrt_mahal)

        if(!is.null(caliper_list)){
            w_dist_mat <- w_dist_mat + caliper_dist_mat
        }

        nonbipartite_matches(dist_mat = w_dist_mat,
                             tolerance_list = tolerance_list,
                             match_method = match_method,
                             n_sinks = n_sinks,
                             keep_all_with_replacement = keep_all_with_replacement)
    })
}
