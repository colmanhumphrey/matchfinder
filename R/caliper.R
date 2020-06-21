#' Simple wrapper to unify caliper input
#'
#' @param caliper_vec Default NULL; numeric vector that "blocks"
#'   matches if they're too (further than \code{caliper_max}) on this value
#' @param caliper_max The maximum allowed difference (exactly this difference
#'   is allowed).
#' @param continuous_mult The value to multiply differences above caliper max.
#'   Set as \code{Inf} to have infinite penalties, i.e. block matches above
#'   the max.
#' @return Either \code{NULL}, or a list with the same names as the input, after checking
#'   values.
#' @author Colman Humphrey
#'
#' @export
gen_caliper_list <- function(caliper_vec = NULL,
                             caliper_max = NULL,
                             continuous_mult = 100) {
    if (is.null(caliper_vec)) {
        if (!is.null(caliper_max)) {
            stop("can't give `caliper_max` without `caliper_vec`")
        }
        return(NULL)
    }

    if (is.null(caliper_max)) {
        stop("supply `caliper_max` if using calipers")
    }

    if (length(caliper_max) > 1L) {
        stop("`caliper_max` should be length one")
    }

    if (length(continuous_mult) > 1L) {
        stop("`continuous_mult` should be length one")
    }

    return(list(
        caliper_vec = caliper_vec,
        caliper_max = caliper_max,
        continuous_mult = continuous_mult
    ))
}


#' Creates caliper penalties
#'
#' Given a vector, \code{caliper_vec}, this function sets penalties
#' for any pairwise differences above \code{caliper_max}; either
#' an \code{Inf} penalty (the default), or a continuous penalty. Typical usage will
#' then add the resulting matrix from this function onto a distance
#' matrix, say from pairwise Mahalanobis.
#'
#' @param caliper_list Result of \code{gen_caliper_list}
#' @param treat_vec Optional; if you only want pairs between treat and control.
#' @return A matrix, either square (\code{|caliper_vec| x |caliper_vec|})
#'   or else \code{sum(treat_vec == 1) x sum(treat_vec == 0)}.
#' @author Colman Humphrey
#'
#' @export
create_caliper <- function(caliper_list,
                           treat_vec = NULL) {
    if (is.null(treat_vec)) {
        caliper_treat <- caliper_list[["caliper_vec"]]
        caliper_control <- caliper_list[["caliper_vec"]]
    } else {
        if (!(length(treat_vec) == length(caliper_list[["caliper_vec"]]))) {
            stop("treat_vec not the same length as `caliper_vec`")
        }
        ## in case it's logical...:
        treat_vec <- treat_vec * 1L

        caliper_treat <- caliper_list[["caliper_vec"]][treat_vec == 1]
        caliper_control <- caliper_list[["caliper_vec"]][treat_vec == 0]
    }

    abs_caliper_diff <- abs(outer(caliper_treat,
        caliper_control,
        FUN = "-"
    )) - caliper_list[["caliper_max"]]
    abs_caliper_diff[abs_caliper_diff < 0] <- 0

    continuous_mult <- caliper_list[["continuous_mult"]]

    if (!is.null(continuous_mult) && continuous_mult < Inf) {
        penalty_mat <- continuous_mult * abs_caliper_diff
    } else {
        penalty_mat <- abs_caliper_diff
        penalty_mat[penalty_mat > 0] <- Inf
    }

    penalty_mat
}
