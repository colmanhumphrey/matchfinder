#' Wraps a match_list from a simple match result and slices off
#' various numbers of sinks
#'
#' Note that for \code{with_replacement} matching, this is the
#' right thing to do, but for greedy matching with a random ordering
#' it isn't obviously correct, but does give the nice property
#' of more sinks = lower distances.
#' @param simple_match_list match result from one of \code{with_replacement_match},
#'   \code{greedy_match}, \code{with_replacement_nbp_match} or
#'   \code{greedy_nbp_match}.
#' @param n_sinks default NULL, vector of sink values to use.
#' @return list of lists; see parent function
#' @author Colman Humphrey
#'
#' @keywords internal
simple_sink_wrap <- function(simple_match_list,
                             n_sinks = NULL) {
    if (is.null(n_sinks)) {
        n_sinks <- 0L
    }

    ## ------------------------------------

    dist_ranks <- rank(simple_match_list[["distance"]],
        ties.method = "random"
    )
    setNames(lapply(n_sinks, function(sink_val) {
        keep_ind <- dist_ranks <= (length(dist_ranks) - sink_val)
        match_list <- lapply(simple_match_list, function(x) {
            x[keep_ind]
        })
        match_list[["num_sinks"]] <- sink_val
        match_list
    }), n_sinks)
}


#' Given a vector of sink values, generates an optimal match
#' for each.
#'
#' Will be slow; you can't just generate one match and subset from it.
#' @inheritParams bipartite_matches
#' @param n_sinks default NULL, vector of sink values to use.
#' @return list of lists; see parent function
#' @author Colman Humphrey
#'
#' @keywords internal
optimal_sink_wrap <- function(dist_mat,
                              treat_vec,
                              n_sinks,
                              tol_val) {
    if (is.null(n_sinks)) {
        n_sinks <- 0L
    }

    ## ------------------------------------

    setNames(lapply(n_sinks, function(sink_val) {
        match_list <- optimal_match(
            dist_mat,
            treat_vec,
            sink_val,
            tol_val
        )
        match_list$num_sinks <- sink_val
        match_list
    }), n_sinks)
}


#' Given a vector of sink values, produces optimal NBP matches for each
#'
#' @inheritParams nonbipartite_matches
#' @author Colman Humphrey
#'
#' @keywords internal
optimal_nbp_sink_wrap <- function(dist_mat,
                                  tolerance_vec,
                                  n_sinks = NULL) {
    if (is.null(n_sinks)) {
        n_sinks <- 0L
    }

    ## ------------------------------------

    setNames(lapply(n_sinks, function(sink_val) {
        match_list <- optimal_nbp_match(dist_mat,
            tolerance_vec = tolerance_vec,
            n_sinks = sink_val
        )
        match_list[["num_sinks"]] <- sink_val
        match_list
    }), n_sinks)
}
