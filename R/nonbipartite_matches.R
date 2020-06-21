#' Simple wrapper to unify tolerance input
#'
#' @param tolerance_vec Default NULL; numeric "continuous treatment"
#'   vector that we use to form
#'   non-bipartite matches: units i and j can be matched if
#'   \eqn{\mid \texttt{tolerance_vec}[i] - \texttt{tolerance_vec}[j]\mid >
#'   \texttt{tolerance_min}}.
#' @param tolerance_min See above for what this does - blocks
#'   matches that are too close on \code{tolerance_vec}. Something like minimum
#'   relevant difference to be a "treatment" effect. Default \code{NULL} gives
#'   zero, i.e. only blocks equality of \code{tolerance_vec}.
#' @param tolerance_max Optionally we may want to also ensure
#'   our "treatment" values aren't too far apart. E.g. we may think our
#'   assumptions are reasonable valid for small differences
#'   in the tolerance vector, but not for large. Or another way:
#'   we're asking for say marginal effects: how bad is one extra beer a
#'   day?
#' @return Either \code{NULL}, or a list with the same names as the input,
#'   with validated values.
#' @author Colman Humphrey
#'
#' @export
gen_tolerance_list <- function(tolerance_vec = NULL,
                               tolerance_min = NULL,
                               tolerance_max = NULL) {
    if (is.null(tolerance_vec)) {
        if (!is.null(tolerance_min) || !is.null(tolerance_max)) {
            stop(
                "can't set `tolerance_{min, max}` without setting ",
                "`tolerance_vec`"
            )
        }
        return(NULL)
    }

    if (is.null(tolerance_min)) {
        tolerance_min <- 0
    } else {
        if (length(tolerance_min) > 1) {
            stop("`tolerance_min` should be length one (for now) or NULL")
        }
    }

    if (!is.null(tolerance_max) && length(tolerance_max) > 1) {
        stop("`tolerance_max` should be length one (for now) or NULL")
    }

    return(list(
        tolerance_vec = tolerance_vec,
        tolerance_min = tolerance_min,
        tolerance_max = tolerance_max
    ))
}
#' Converts tolerance list to caliper list
#' @inheritParams gen_caliper_list
#' @param tolerance_list Result from \code{gen_tolerance_list}
#' @param use_min Logical, should the caliper max be the tolerance min?
#'   Use tolerance max if not. Default \code{TRUE}.
#' @return List from \code{gen_caliper_list}
#' @author Colman Humphrey
#'
#' @keywords internal
tolerance_to_caliper_list <- function(tolerance_list,
                                      use_min = TRUE,
                                      continuous_mult = 1) {
    stopifnot(is_tf(use_min))

    cal_max <- ifelse(use_min,
        tolerance_list[["tolerance_min"]],
        tolerance_list[["tolerance_max"]]
    )
    gen_caliper_list(
        caliper_vec = tolerance_list[["tolerance_vec"]],
        caliper_max = cal_max,
        continuous_mult = continuous_mult
    )
}

#' Generating nonbipartite matched pairs
#'
#' Generates matched pairs either:
#' \describe{
#'   \item{With Replacement}{Finds smallest control for each treatment}
#'   \item{Without Replacement, Greedy}{Greedily generates pairs. Note that
#'   the order for choosing the greedy pairs is random, which is not the only
#'   possible solution.}
#'   \item{Without Replacement, Optimally}{Minimum total distance}
#' }
#' If you're happy to allow units to be the control multiple times,
#' then the first way is fast and optimal.
#'
#' If not, you have to trade off speed vs optimality. Greedy runs
#' over all units in a random order, so if you want to run greedy a bunch of
#' times and take the best, it would still be (likely) much faster than
#' running optimal matching.
#'
#' NOTE:
#' with sinks -
#' the optimal match gives "half" the number of sinks you might expect.
#' If you have a pool of 400 matchable things,
#' NBP gives you 200 matches if no sinks.
#' If you say sinks = 20 FOR OPTIMAL,
#' you'll get 190 matches, not 180.
#' We leave this behaviour in so that this algorithm
#' can generate expected results if someone is used to nbpmatch
#'
#' @inheritParams all_nonbipartite_matches
#' @inheritParams bipartite_matches
#' @param n_sinks Vector of sinks per match. Note that this is NOT the same
#'   for each \code{match_method}: for both \code{"with_replacement"}
#'   and \code{"greedy"}, it subtracts one pair. But for optimal matching,
#'   it removes one full unit. \code{"greedy"} is the weird one here, but it
#'   wouldn't truly make sense given the code to try and replicate the optimal
#'   implementation. Be safe and don't use greedy...
#' @param keep_all_with_replacement logical, default FALSE.
#'   When nbp matching with replacement, you can in some cases form nearly
#'   as many pairs as there are units (if using a tolerance vec, then maybe not
#'   the lowest value pairs, since they'll have no control unit).
#'   But this would form a contrast with all other methods, so we will cut
#'   down to using half by default.
#' @return basic return value is a list with three elements:
#'   \describe{
#'     \item{\code{treat_index}}{index of treated units}
#'     \item{\code{control_index}}{index of control units}
#'     \item{\code{distance}}{distances between the pairs}
#'   }
#'   If \code{n_sinks} is not NULL, you'll get a list of such objects, each
#'   with an extra element: the number of sinks used.
#' @author Colman Humphrey
#' @export
nonbipartite_matches <- function(dist_mat,
                                 tolerance_list = gen_tolerance_list(),
                                 match_method = c(
                                     "with_replacement",
                                     "optimal",
                                     "greedy"
                                 ),
                                 n_sinks = NULL,
                                 keep_all_with_replacement = FALSE) {
    stopifnot(is.matrix(dist_mat))
    stopifnot(min(dist_mat) >= 0)

    match_method <- match.arg(match_method)

    if (nrow(dist_mat) != ncol(dist_mat) ||
        sum(abs(dist_mat - t(dist_mat)), na.rm = TRUE) >
            (0.001 * max(ifelse(dist_mat == Inf, 0, dist_mat)))) {
        stop("dist_mat should be square and symmetric", call. = FALSE)
    }

    ## dealing with tolerance:
    if (!is.null(tolerance_list)) {
        if (length(tolerance_list[["tolerance_vec"]]) != nrow(dist_mat)) {
            stop("`tolerance_vec` must be the same length as the ",
                "dimensions of `dist_mat`",
                call. = FALSE
            )
        }

        temp_dist <- create_caliper(
            caliper_list = tolerance_to_caliper_list(tolerance_list,
                use_min = TRUE
            ),
            treat_vec = NULL
        )
        dist_mat[temp_dist == 0] <- Inf

        if (!is.null(tolerance_list[["tolerance_max"]])) {
            temp_dist <- create_caliper(
                caliper_list = tolerance_to_caliper_list(tolerance_list,
                    use_min = FALSE
                ),
                treat_vec = NULL
            )
            dist_mat[temp_dist > 0] <- Inf
        }
    } else {
        if (missing(tolerance_list)) {
            warning("assuming that all pairs are matchable, ",
                "not neccesarily correct. If this is what you want, ",
                "you can silence this warning by explicitly supplying ",
                "`tolerance_vec = NULL`",
                call. = FALSE
            )
        }

        ## no self-matching
        diag(dist_mat) <- Inf
    }

    stopifnot(is_tf(keep_all_with_replacement))

    if (!is.null(n_sinks)) {
        stopifnot(is.numeric(n_sinks) &&
            min(n_sinks) >= 0L &&
            !any(is.na(n_sinks)) &&
            length(unique(n_sinks)) == length(n_sinks))
    }

    ## ------------------------------------

    if (match_method == "with_replacement") {
        return(simple_sink_wrap(
            with_replacement_nbp_match(
                dist_mat,
                tolerance_list[["tolerance_vec"]],
                keep_all = keep_all_with_replacement
            ),
            n_sinks
        ))
    }

    if (match_method == "greedy") {
        return(simple_sink_wrap(
            greedy_nbp_match(
                dist_mat,
                tolerance_list[["tolerance_vec"]]
            ),
            n_sinks
        ))
    }

    optimal_nbp_sink_wrap(dist_mat,
        tolerance_vec = NULL,
        n_sinks = n_sinks
    )
}


#' @inheritParams bipartite_matches
#' @keywords internal
with_replacement_nbp_match <- function(dist_mat,
                                       tolerance_vec,
                                       keep_all = FALSE) {
    tol_small <- which(outer(tolerance_vec, tolerance_vec, "-") <= 0,
        arr.ind = TRUE
    )
    dist_mat[tol_small] <- Inf

    min_index <- min_different_rank(dist_mat)

    match_list <- list(
        treat_index = 1:nrow(dist_mat),
        control_index = min_index,
        distance = dist_mat[cbind(
            1:nrow(dist_mat),
            min_index
        )]
    )

    if (!keep_all) {
        ## while this seems arbitrary, not limiting to half would mean this
        ## would give twice as many pairs as all other results
        keep_ind <- rank(match_list[["distance"]],
            ties.method = "random"
        ) <= (nrow(dist_mat) / 2)
        match_list <- lapply(match_list, function(x) {
            x[keep_ind]
        })
    }

    ## might also be pointless now
    reorder_nbp(match_list, tolerance_vec)
}
#' @inheritParams with_replacement_nbp_match
#' @keywords internal
remove_duplicates <- function(match_list) {
    min_ind <- pmin(
        match_list[["treat_index"]],
        match_list[["control_index"]]
    )
    max_ind <- pmax(
        match_list[["treat_index"]],
        match_list[["control_index"]]
    )

    duplicated_ind <- duplicated(cbind(min_ind, max_ind))

    lapply(match_list, function(x) {
        x[!duplicated_ind]
    })
}

#' @inheritParams bipartite_matches
#' @keywords internal
greedy_nbp_match <- function(dist_mat,
                             tolerance_vec) {
    min_vals <- apply(dist_mat, 1, min)

    result_mat <- matrix(NA, nrow = nrow(dist_mat), ncol = 3)

    while (min(min_vals) < Inf) {
        random_value <- sample(1:length(min_vals),
            size = 1,
            prob = 1 / (min_vals + 1)
        )
        match_ind <- which(rank(dist_mat[random_value, ],
            ties.method = "random"
        ) == 1L)

        if (tolerance_vec[random_value] < tolerance_vec[match_ind]) {
            result_mat[random_value, 1] <- match_ind
            result_mat[random_value, 2] <- random_value
        } else {
            result_mat[random_value, 1] <- random_value
            result_mat[random_value, 2] <- match_ind
        }
        result_mat[random_value, 3] <- dist_mat[random_value, match_ind]

        ## blocked:
        dist_mat[, match_ind] <- Inf
        dist_mat[match_ind, ] <- Inf
        dist_mat[, random_value] <- Inf
        dist_mat[random_value, ] <- Inf
        min_vals <- apply(dist_mat, 1, min)
    }
    result_mat <- result_mat[!is.na(result_mat[, 1]), ]

    list(
        treat_index = result_mat[, 1],
        control_index = result_mat[, 2],
        distance = result_mat[, 3]
    )
}
#' @inheritParams nonbipartite_matches
#' @keywords internal
optimal_nbp_match <- function(dist_mat,
                              tolerance_vec = NULL,
                              n_sinks = 0) {
    ## create the match (timely usually)
    nbp_match <- nbpMatching::nonbimatch(
        add_nbp_sinks(
            dist_mat = dist_mat,
            n_sinks = n_sinks
        )
    )

    ## fix it; remove sinks etc etc
    fix_nbp_match(nbp_match,
        nrow_match = nrow(dist_mat),
        tolerance_vec
    )
}


#' This function gives a distance matrix that nbpmatching likes
#'
#' Note: if your \code{dist_mat} is even (e.g. number of rows/cols),
#' you should supply an even number of sinks, and odd if odd -
#' else you'll get a ghost value added and a warning. This is true
#' even for adding zero sinks to an odd.
#' @param dist_mat A symmetric distance matrix, e.g. result of
#'   \code{weighted_mahal}
#' @param n_sinks How many potential matches to throw away?
#' @keywords internal
add_nbp_sinks <- function(dist_mat,
                          n_sinks = 0L) {
    if ((nrow(dist_mat) + n_sinks) %% 2L != 0L) {
        warning(
            "There must be an even number of elements (including sinks); ",
            "adding an extra sink"
        )
        n_sinks <- n_sinks + 1L
    }

    if (n_sinks > 0) {
        distmat_add_zeros <- matrix(
            0,
            nrow(dist_mat) + n_sinks,
            ncol(dist_mat) + n_sinks
        )
        distmat_add_zeros[
            1:nrow(dist_mat),
            1:ncol(dist_mat)
        ] <- dist_mat
        ## they can't match each other:
        distmat_add_zeros[
            1:n_sinks + nrow(dist_mat),
            1:n_sinks + ncol(dist_mat)
        ] <- Inf
    } else {
        distmat_add_zeros <- dist_mat
    }

    ## get in format
    nbpMatching::distancematrix(distmat_add_zeros)
}


#' Converting nbp output to common output for this package.
#'
#' Takes result of NBP matching and deletes sinks,
#' and I guess infinite distances; finally also returns only some columns.
#' @param nbp_match Result from nbpMatching::nonbimatch.
#' @param nrow_match How many real rows were there?
#'   You feed sinks to nbpmatching sort of manually, so
#'   this is one way to remove them.
#' @param tolerance_vec If this is given, we sort the matches by high / low.
#' @keywords internal
fix_nbp_match <- function(nbp_match,
                          nrow_match,
                          tolerance_vec = NULL) {
    halves <- nbp_match[["halves"]][, c("Group1.Row", "Group2.Row", "Distance")]
    names(halves) <- c("treat_index", "control_index", "distance")

    nonphantom_ind <- pmax(
        halves[["treat_index"]],
        halves[["control_index"]]
    ) <= nrow_match
    noninf_ind <- halves[["distance"]] != Inf

    halves <- halves[nonphantom_ind & noninf_ind, ]

    rownames(halves) <- NULL

    reorder_nbp(
        as.list(halves),
        tolerance_vec
    )
}
#' Reorders list first by treat index, then within pairs by
#' highest tol val if given
#'
#' @inheritParams nonbipartite_matches
#'
#' @keywords internal
reorder_nbp <- function(match_list,
                        tolerance_vec = NULL) {
    treat_order <- order(match_list[["treat_index"]])
    match_list <- lapply(match_list, function(x) {
        x[treat_order]
    })

    if (!is.null(tolerance_vec)) {
        larger_ind <- tolerance_vec[match_list[["treat_index"]]] >
            tolerance_vec[match_list[["control_index"]]]

        return(list(
            treat_index = ifelse(larger_ind,
                match_list[["treat_index"]],
                match_list[["control_index"]]
            ),
            control_index = ifelse(larger_ind,
                match_list[["control_index"]],
                match_list[["treat_index"]]
            ),
            distance = match_list[["distance"]]
        ))
    }

    match_list
}
