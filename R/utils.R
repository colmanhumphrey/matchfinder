#' Minimum ranked element, not equal to the row index
#'
#' @param dist_mat (Probably symmetric) matrix.
#' @return vector of indices
#' @author Colman Humphrey
#'
#' @keywords internal
min_different_rank <- function(dist_mat) {
    min_blocked_rank(
        dist_mat,
        seq_len(nrow(dist_mat))
    )
}


#' Minimum ranked element, not equal to the blocked index
#'
#' @param dist_mat Some distance matrix.
#' @param blocked_ind For each row, blocked element.
#' @return Vector of indices
#' @author Colman Humphrey
#'
#' @keywords internal
min_blocked_rank <- function(dist_mat,
                             blocked_ind = NULL) {
    if (is.null(blocked_ind)) {
        return(unlist(lapply(seq_len(nrow(dist_mat)), function(j) {
            ranks <- rank(dist_mat[j, ], ties.method = "random")
            which(ranks == min(ranks))
        })))
    }

    stopifnot(length(blocked_ind) == nrow(dist_mat))

    unlist(lapply(seq_len(nrow(dist_mat)), function(j) {
        ranks <- rank(dist_mat[j, ], ties.method = "random")
        ranks[blocked_ind[j]] <- max(ranks) + 1L
        which(ranks == min(ranks))
    }))
}


#' example from ?sample, when length could be one
#' @export
fixed_sample <- function(x, ...) x[sample.int(length(x), ...)]


#' Takes a subsample with a given number of unique elements
#' @param primary vector to subsample from
#' @param unique how many unique elements you want
#' @return vector equal in length to primary, with `unique` unique elements
#' @author Colman Humphrey
#' @keywords internal
unique_size_sub <- function(primary, unique) {
    sub_sample <- fixed_sample(primary, unique, replace = FALSE)

    c(
        sub_sample,
        fixed_sample(sub_sample, length(primary) - unique, replace = TRUE)
    )
}


#' Test if \code{x} is a length-one logical
#'
#' @param x ideally a logical!
#' @return TRUE if x is a length one logical, else FALSE
#' @author Colman Humphrey
#'
#' @keywords internal
is_tf <- function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
}


#' Produces a list of indices split into k-folds
#'
#' @param index_length length of vector to index
#' @param num_folds how many folds you need
#' @return list of vectors of indices
#' @author Colman Humphrey
#'
#' @export
fold_indexing <- function(index_length,
                          num_folds) {
    random_index <- sample(1L:index_length, size = index_length)

    gap_len <- index_length %% num_folds
    sub_length <- index_length - gap_len

    cut_indices <- seq(1L, sub_length + 1L,
        length.out = num_folds + 1L
    ) +
        cumsum(c(0L, rep(1L, gap_len), rep(0L, num_folds - gap_len)))

    lapply(1L:num_folds, function(j) {
        random_index[cut_indices[j]:(cut_indices[j + 1L] - 1L)]
    })
}


#' Simple brier score calc
#'
#' @param predict \eqn{P(Y = 1)} for each value
#' @param outcome result (in \eqn{\{0, 1\}})
#' @param avg logica, default TRUE, do you want the mean?
#' @return length one double: the total brier sum, or the avg (default)
#' @author Colman Humphrey
#' @keywords internal
calc_brier <- function(predict,
                       outcome,
                       avg = TRUE) {
    if (!avg) {
        return(sum((outcome - predict)^2))
    } else {
        return(mean((outcome - predict)^2))
    }
}


#' Given a logical swap vector, switches treated and control units around
#'
#' @param match_list typical \code{match_list} entity
#' @param swap logical vector indicating which pairs should swap
#' @return another match_list, but now with swapped pairs
#' @author Colman Humphrey
swap_pairs <- function(match_list,
                       swap) {
    stopifnot(length(swap) == length(match_list[["treat_index"]]))

    primary_list <- list(
        treat_index = ifelse(!swap,
            match_list[["treat_index"]],
            match_list[["control_index"]]
        ),
        control_index = ifelse(swap,
            match_list[["treat_index"]],
            match_list[["control_index"]]
        )
    )

    if (all(c("treat_index_within", "control_index_within") %in%
        names(match_list))) {
        primary_list[["treat_index_within"]] <-
            ifelse(!swap,
                match_list[["treat_index_within"]],
                match_list[["control_index_within"]]
            )

        primary_list[["control_index_within"]] <-
            ifelse(swap,
                match_list[["treat_index_within"]],
                match_list[["control_index_within"]]
            )

        primary_list <- primary_list[c(1L, 3L, 2L, 4L)]
    }

    primary_list[["distance"]] <- match_list[["distance"]]

    primary_list
}

#' Converts \code{rank_cols} in all allowed forms to an integer index
#'
#' Takes in ranked cols, either given by integer (numeric is fine) index;
#' logical index; or named. Converts to integer index
#' @param rank_cols integer/number, or logical, or names within
#'   \code{colnames(x_mat)}
#' @param x_mat the x matrix of interest
#' @author Colman Humphrey
#' @keywords internal
rank_integer_index <- function(rank_cols, x_mat) {
    if (is.null(rank_cols)) {
        return(vector("integer", 0L))
    }

    if (all(rank_cols %in% 1L:ncol(x_mat))) {
        return(rank_cols)
    }

    if (is.logical(rank_cols)) {
        if (length(rank_cols) == ncol(x_mat)) {
            return(which(rank_cols))
        } else {
            stop("logical `rank_cols` must be same length as `ncol(x_mat)`")
        }
    } else {
        if (is.null(colnames(x_mat))) {
            stop("x_mat must have colnames to use named rank_cols")
        }

        if (any(!(rank_cols %in% colnames(x_mat)))) {
            stop("not all rank_cols are present in x_mat colnames")
        }

        return(which(colnames(x_mat) %in% rank_cols))
    }
}


#' Converts indicated columns to ranked versions of themselves
#'
#' This function takes a numeric matrix and converts any columns
#' indicated by the \code{rank_cols} input to their ranks
#' (break ties however you want), scaled down by \code{nrow(x_mat)}.
#' @param x_mat numeric matrix (adjust non-numeric columns prior)
#' @param rank_cols names or index of columns to be converted to ranks
#'   before analysis
#' @param ties_method how to break ties in ranks, by default uses
#'   "average" (same as the rank function's default). See \code{?rank}
#'   for further options
#' @return x_mat again, potentially adjusted for ranks
#' @author Colman Humphrey
#' @keywords internal
ranked_x <- function(x_mat,
                     rank_cols = NULL,
                     ties_method = "average") {
    ## convert to ranks
    for (rank_change in rank_integer_index(rank_cols, x_mat)) {
        ## divides by nrow, better scaling
        x_mat[, rank_change] <- rank(x_mat[, rank_change],
            ties.method = ties_method
        ) / nrow(x_mat)
    }

    x_mat
}


#' Find nearest index to vector given set of indices
#' @param find_vec vector to find close matches in
#' @param given_index index for which we need to find the friends
#' @return set of indices of values closest to each given index
#' @author Colman Humphrey
near_given_match <- function(find_vec,
                             given_index) {
    ## runif breaks ties randomly between orders for same element
    order_tol <- order(
        find_vec,
        runif(length(find_vec))
    )

    given_order_index <- match(given_index, order_tol)

    ## ----------------
    ## find the nearest elements above our given_index

    ## the pmin is to avoid indexing above the vector
    ## won't happen in a proper match, since all controls need a higher
    ## tol treat
    ## but just in case / in cases with many ties or whatever
    above_index <- order_tol[pmin(
        given_order_index + 1L,
        length(order_tol)
    )]
    tol_above <- find_vec[above_index]
    ## if our control is the max, we won't use above... of course
    tol_above[given_order_index == length(order_tol)] <- Inf

    ## ----------------
    ## find the nearest elements below our control index

    ## here, we do need this pmax even in proper matches:
    ## it's very possible to use the smallest tolerance value as a control
    below_index <- order_tol[pmax(given_order_index - 1L, 1L)]
    tol_below <- find_vec[below_index]
    ## if our control is the min, we won't be using below... of course
    tol_below[given_order_index == 1L] <- Inf

    ## ----------------
    ## note that the above Infs cannot co-occur, so we're good

    given_tolerance <- find_vec[given_index]
    ifelse(abs(tol_above - given_tolerance) > abs(given_tolerance - tol_below),
        below_index, above_index
    )
}

#' exp / 1 + exp
#' @param x numeric vector
#' @return numeric vector, now in (0, 1)
#' @author Colman Humphrey
expit <- function(x) {
    exp(x) / (1 + exp(x))
}


#' Symmetrises a matrix in a boring way, and zeros the diag
#'
#' This is useful for generating pairwise distance matrices
#' @param mat a matrix
#' @keywords internal
sym_mat <- function(mat) {
    sym <- mat + t(mat)
    diag(sym) <- 0
    sym
}


#' Generic binary search algo: finding the input.
#'
#' @param target_value The value the function should achieve.
#' @param monotone_function The function we want to find the
#'   relevant input for.
#' @param init_bounds Default NULL; supply a length-two vector if
#'   you wish to supply your own initial bounds. Note that if your function
#'   only works in a certain range, or is only monotonic in a certain range,
#'   then you'll very likely need to supply your own bounds. If not given,
#'   we'll search in both directions from zero.
#' @param error_gap Binary search until the gap between the
#'   function on the input and the target value is smaller than this number.
#' @param max_iters How many iterations to try at most. You can supply
#'   \code{Inf} if you really want.
#' @return The input that gives \code{target_value} as output.
#' @author Colman Humphrey
#'
#' @export
binary_search <- function(target_value,
                          monotone_function,
                          init_bounds = NULL,
                          error_gap = 1e-6,
                          max_iters = 100L) {
    stopifnot(error_gap > 0)

    test_bounds <- init_bounds
    if (is.null(init_bounds)) {
        test_bounds <- c(0, 1)
    }
    if (monotone_function(test_bounds[1L]) >
        monotone_function(test_bounds[2L])) {
        return(binary_search(
            target_value = -target_value,
            monotone_function = function(x) {
                -monotone_function(x)
            },
            init_bounds = init_bounds,
            max_iters = max_iters
        ))
    }

    if (is.null(init_bounds)) {
        if (monotone_function(0) > target_value) {
            upper_bound <- 0
            lower_bound <- -1
            while (monotone_function(lower_bound) > target_value) {
                lower_bound <- lower_bound * 1.1 - 1
            }
        } else {
            lower_bound <- 0
            upper_bound <- 1
            while (monotone_function(upper_bound) < target_value) {
                upper_bound <- upper_bound * 1.1 + 1
            }
        }
    } else {
        stopifnot(length(init_bounds) == 2L)
        stopifnot(init_bounds[1L] < init_bounds[2L])

        lower_bound <- init_bounds[1]
        upper_bound <- init_bounds[2]
    }

    stopifnot(sign(monotone_function(lower_bound) - target_value) *
        sign(monotone_function(upper_bound) - target_value) == -1L)

    input_val <- (lower_bound + upper_bound) / 2
    mono_val <- monotone_function(input_val)

    iters <- 0L

    while (abs(mono_val - target_value) > error_gap && iters < max_iters) {
        if (mono_val > target_value) {
            upper_bound <- input_val
        } else {
            lower_bound <- input_val
        }

        input_val <- (lower_bound + upper_bound) / 2
        mono_val <- monotone_function(input_val)

        iters <- iters + 1L
    }

    if (iters == max_iters) {
        stop("reached max iterations")
    }

    input_val
}


##' Validates a match on a tolerance list
##'
##' Checks that the differences are all at least
##' at the min, and if a max is given, checks that too.
##' Gives reasonable error messages on errors
##' @inheritParams match_estimate_tolerance
##' @return List with two elements:
##' \itemize{
##'  \item{\code{error}}{Boolean - is all good?}
##'  \item{\code{message}}{If error, then gives error string}
##' }
##' @author Colman Humphrey
##'
##' @keywords internal
tolerance_check <- function(match_list,
                            tolerance_list) {
    tol_vec <- tolerance_list[["tolerance_vec"]]

    tol_diffs <- tol_vec[match_list[["treat_index"]]] -
        tol_vec[match_list[["control_index"]]]

    ## equality isn't allowed
    min_violations <- sum(tol_diffs <= tolerance_list[["tolerance_min"]])
    if (min_violations > 0L) {
        neg_violations <- sum(tol_diffs < 0)
        if (neg_violations > 0L) {
            error_message <-
                if (neg_violations == length(tol_diffs)) {
                    ## all wrong direction
                    paste0(
                        "tolerance direction violated by match: ",
                        "all pairs have treatment tolerance ",
                        "with lower value than control tolerance, ",
                        "did the treatment and control get switched?"
                    )
                } else {
                    paste0(
                        "tolerance direction violated by match: ",
                        neg_violations, " pairs have treatment tolerance ",
                        "with lower value than control tolerance"
                    )
                }
            if (neg_violations < min_violations) {
                error_message <- paste0(
                    error_message,
                    "; a further ", min_violations - neg_violations,
                    " pairs have min constraint violated by match ",
                    "(differences below or at `tolerance_min`)"
                )
            }

            return(list(
                error = TRUE,
                message = error_message
            ))
        }

        return(list(
            error = TRUE,
            message = paste0(
                "",
                min_violations, " pairs have difference below or at ",
                "`tolerance_min`"
            )
        ))
    }


    if (!is.null(tolerance_list[["tolerance_max"]])) {
        ## equality is allowed
        max_violations <- sum(tol_diffs > tolerance_list[["tolerance_max"]])

        if (max_violations > 0L) {
            return(list(
                error = TRUE,
                message = paste0(
                    "max constraint violated by match: ",
                    max_violations, " pairs have difference above ",
                    "`tolerance_max`"
                )
            ))
        }
    }

    return(list(
        error = FALSE,
        message = ""
    ))
}
