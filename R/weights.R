#' Generating random weights to use in matching procedures,
#' each summing to one.
#'
#' This function generates a vector of weights for each
#' column in your input data. You must set prior weights
#' (which can all be the same), and you must set a number of vectors.
#' After that you have some options: set minimum weights to make sure
#' that we don't have too little of any given column,
#' and set up hierarchical info if you have categories you
#' want to lump together (see below)
#' @param prior_weights Must be equal to the length of your columns -
#'   i.e. the length of the weight vectors this function will produce.
#'   We'll use this to generate uniform random variables from 0 to this
#'   number. Can be the same value repeated if desired.
#' @param number_vectors How many weight vectors you want.
#' @param minimum_weights If you want to set minimums, either globally or
#'   per column. Note that this will give this minimum, and then add random
#'   weights on top of that.
#' @param hierarchical_list List per group / category of variable:
#'   \describe{
#'     \item{\code{"index"}}{Vector of indices that this group corresponds with}
#'     \item{\code{"weight"}}{Weight for this group}
#'     \item{\code{"variance"}}{(Optional) how much variance you want this group
#'   to have}
#'   }
#'   If you use this, we'll still combine it with
#'   `prior_weights` and `minimum_weights`
#' @return list of weight vectors
#' @author Colman Humphrey
#'
#' @export
generate_random_weights <- function(prior_weights,
                                    number_vectors,
                                    minimum_weights = NULL,
                                    hierarchical_list = NULL) {
    if (is.null(minimum_weights)) {
        minimum_weights <- rep(0, length(prior_weights))
    } else {
        if (length(minimum_weights) == 1) {
            minimum_weights <- rep(minimum_weights, length(prior_weights))
        }

        if (length(minimum_weights) != length(prior_weights)) {
            stop(
                "`minimum_weights` must be length 1 ",
                "or same length as `prior_weights`"
            )
        }
    }

    if (sum(minimum_weights) > 1) {
        stop(
            "`sum(minimum_weights)` must be less than 1, ",
            "since we're returning scaled weights"
        )
    }

    stopifnot(length(number_vectors) == 1L &&
        is.numeric(number_vectors) &&
        number_vectors >= 1L)

    ## ------------------------------------

    if (!is.null(hierarchical_list)) {
        all_index <- unlist(lapply(hierarchical_list, function(x) {
            x[["index"]]
        }))

        if (!isTRUE(all.equal(
            sort(all_index),
            1L:length(prior_weights)
        ))) {
            stop(
                "the hierarchical list must exactly cover ",
                "`1L:length(prior_weights)`"
            )
        }

        hierarchical_list <- lapply(hierarchical_list, function(x) {
            if (is.null(x[["variance"]])) {
                ## scale invariance
                x[["gamma_alpha"]] <- 1
                x[["gamma_beta"]] <- 1 / x[["weight"]]
            } else {
                x[["gamma_beta"]] <- x[["weight"]] / x[["variance"]]
                x[["gamma_alpha"]] <- x[["weight"]] * x[["gamma_beta"]]
            }
            return(x)
        })

        return(hierarchical_random_weights(
            prior_weights = prior_weights,
            number_vectors = number_vectors,
            minimum_weights = minimum_weights,
            hierarchical_list
        ))
    }

    rescale_to_min <- (function(sum_weights, minimum_weights) {
        function(weight_vec) {
            (weight_vec / sum(weight_vec)) *
                (1 - sum_weights) + minimum_weights
        }
    })(sum(minimum_weights), minimum_weights)

    lapply(1L:number_vectors, function(j) {
        temp_weight <- runif(length(prior_weights), 0, prior_weights)
        rescale_to_min(temp_weight)
    })
}


#' Generates weights in a hierarchical setting
#' @inheritParams generate_random_weights
#' @keywords internal
hierarchical_random_weights <- function(prior_weights,
                                        number_vectors,
                                        minimum_weights,
                                        hierarchical_list) {
    index_list <- lapply(hierarchical_list, function(cat) {
        cat[["index"]]
    })

    order_vec <- order(unlist(index_list))

    rescale_to_min <- (function(sum_weights, minimum_weights) {
        function(weight_vec) {
            (weight_vec / sum(weight_vec)) *
                (1 - sum_weights) + minimum_weights
        }
    })(sum(minimum_weights), minimum_weights)

    lapply(1L:number_vectors, function(j) {
        weight_list <- lapply(hierarchical_list, function(cat) {
            gamma_var <- rgamma(1,
                shape = cat[["gamma_alpha"]],
                rate = cat[["gamma_beta"]]
            )
            pre_weights <- runif(
                length(cat[["index"]]),
                0,
                prior_weights[cat[["index"]]]
            )
            (pre_weights / sum(pre_weights)) * gamma_var
        })
        weight_vec <- unlist(weight_list)[order_vec]

        rescale_to_min(weight_vec)
    })
}
