#' Calculates the standard deviation for a match, on the outcome
#'
#' @inheritParams all_bipartite_matches
#' @param y_vector Outcome vector (not used in match generation).
#' @param match_list A particular match result.
#' @param use_all_controls logical; if using bipartite matches, should we
#'   estimate the variance using all controls possible? Default \code{TRUE}.
#'   If \code{FALSE}, will only use controls from the match given
#' @return a single float, the standard deviation (not standard error)
#'   for the match
#' @author Colman Humphrey
#'
#' @export
bipartite_match_sd <- function(x_mat,
                               cov_x,
                               y_vector,
                               match_list,
                               treat_vec,
                               caliper_list = gen_caliper_list(),
                               propensity_list = match_propensity_list(NULL),
                               weight_vec = NULL,
                               use_all_controls = TRUE,
                               sqrt_mahal = TRUE) {
    if (!is.null(caliper_list) && !is.null(propensity_list)) {
        stop(
            "don't use both `caliper_list` and `propensity_list`: ",
            " If you do want both, create the combined caliper separately"
        )
    }
    ## ------------------------------------
    ## two variance components

    ## ------------------------------------
    ## V1 easy: usual variance
    ## exists if matches are made without replacement

    var_difference <- var(y_vector[match_list[["treat_index"]]] -
        y_vector[match_list[["control_index"]]])

    ## ----------------
    ## V2 harder: the added variance from using repeats

    all_unique_controls <- !any(duplicated(match_list[["control_index"]]))

    if (all_unique_controls) {
        var_repeated <- 0
    } else {
        if (!is.null(propensity_list)) {
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
                caliper_max = sd(prop_score) *
                    propensity_list[["caliper_sd_mult"]],
                continuous_mult = propensity_list[["continuous_mult"]]
            )
        } # else caliper list is as it was

        var_repeated <- gen_bipartite_repeated_variance(
            x_mat,
            cov_x = cov_x,
            y_vector = y_vector,
            control_index = match_list[["control_index"]],
            treat_vec = treat_vec,
            caliper_list = caliper_list,
            weight_vec = weight_vec,
            use_all_controls = use_all_controls,
            sqrt_mahal = sqrt_mahal
        )
    }

    ## ------------------------------------

    var_total <- var_difference + var_repeated

    sqrt(var_total)
}


#' Abadie, Alberto and Guido W Imbens (2006).
#' “Large sample properties of matching estimators
#' for average treatment effects”.
#' Econometrica 74.1, pp. 235–267. - see page 251.
#' More info at “Matching on the estimated propensity score”.
#' Econometrica 84.2, pp. 781–807.
#' @inheritParams bipartite_match_sd
#' @keywords internal
gen_bipartite_repeated_variance <- function(x_mat,
                                            cov_x,
                                            y_vector,
                                            control_index,
                                            treat_vec,
                                            caliper_list = gen_caliper_list(),
                                            weight_vec = NULL,
                                            use_all_controls = TRUE,
                                            sqrt_mahal = TRUE) {
    all_control_index <- which(treat_vec == 0L)

    count_frame <- data.frame(table(control_index))
    names(count_frame) <- c("control_index", "count")

    ## the second term should actually be K_{sq,i},
    ## but when using just one control, we have K_i = K_{sq,i}
    ## 2020-10-05: I don't really understand the above comment
    k_sq_minus_k <- count_frame[["count"]]^2 - count_frame[["count"]]

    unique_control_index <- as.numeric(as.character(
        count_frame[["control_index"]]
    ))

    use_all_controls <- ifelse(is.null(use_all_controls),
        TRUE, use_all_controls
    )

    ## ------------------------------------

    control_within_all_index <- match(
        unique_control_index,
        all_control_index
    )

    if (use_all_controls) {
        control_controls <- all_control_index
    } else {
        control_controls <- unique_control_index
    }

    blocked_ind <- match(
        unique_control_index,
        control_controls
    )
    partial_index <- list(
        unique_control_index,
        control_controls
    )

    control_dist_mat <- weighted_mahal(
        x_mat = x_mat,
        cov_x = cov_x,
        weight_vec = weight_vec,
        sqrt_mahal = sqrt_mahal,
        partial_index = partial_index
    )

    if (!is.null(caliper_list)) {
        caliper_full <- create_caliper(caliper_list)
        control_dist_mat <- control_dist_mat +
            caliper_full[cbind(
                partial_index[[1]],
                partial_index[[2]]
            )]
    }

    min_control_match <- min_blocked_rank(
        control_dist_mat,
        blocked_ind
    )

    control_match <- control_controls[min_control_match]

    ## ------------------------------------

    sigma_xw_sq <- 0.5 * (y_vector[unique_control_index] -
        y_vector[control_match])^2

    ## total variance due to rep
    sum(k_sq_minus_k * sigma_xw_sq) / (length(control_index) - 1)
}

#' @inheritParams match_standard_scaled_dev
#' @keywords internal
gen_nonbipartite_repeated_variance <- function(x_mat,
                                               cov_x,
                                               y_vector,
                                               control_index,
                                               tolerance_list =
                                                   gen_tolerance_list(),
                                               caliper_list =
                                                   gen_caliper_list(),
                                               weight_vec = NULL,
                                               use_all_controls = TRUE,
                                               sqrt_mahal = TRUE) {
    count_frame <- data.frame(table(control_index))
    names(count_frame) <- c("control_index", "count")

    ## the second term should actually be K_{sq,i},
    ## but when using just one control, we have K_i = K_{sq,i}
    k_sq_minus_k <- count_frame[["count"]]^2 - count_frame[["count"]]

    unique_control_index <- as.numeric(as.character(
        count_frame[["control_index"]]
    ))

    ## ------------------------------------

    partial_index <- list(
        unique_control_index,
        seq_len(nrow(x_mat))
    )

    control_dist_mat <- weighted_mahal(
        x_mat = x_mat,
        cov_x = cov_x,
        weight_vec = weight_vec,
        sqrt_mahal = sqrt_mahal,
        partial_index = partial_index
    )

    if (!is.null(caliper_list)) {
        caliper_full <- create_caliper(caliper_list)
        control_dist_mat <- control_dist_mat +
            caliper_full[cbind(
                partial_index[[1]],
                partial_index[[2]]
            )]
    }

    min_control_match <- min_blocked_rank(
        control_dist_mat,
        blocked_ind = unique_control_index
    )

    ## ----------------

    tolerance_match <- near_given_match(
        tolerance_list[["tolerance_vec"]],
        given_index = unique_control_index
    )

    ## ----------------

    avg_y_control_controls <- (y_vector[min_control_match] +
        y_vector[tolerance_match]) / 2
    avg_y_control_tols <- (
        tolerance_list[["tolerance_vec"]][min_control_match] +
            tolerance_list[["tolerance_vec"]][tolerance_match]) / 2

    y_diffs <- y_vector[unique_control_index] -
        avg_y_control_controls
    tol_diffs <- tolerance_list[["tolerance_vec"]][unique_control_index] -
        avg_y_control_tols

    ## somewhat arbitrary work here:
    if (!is.null(tolerance_list[["tolerance_min"]])) {
        tol_diffs <- pmax(tol_diffs, tolerance_list[["tolerance_min"]])
    } else {
        sorted_tol <- sort(unique(tolerance_list[["tolerance_vec"]]))
        min_tol_diff <- min(sorted_tol[2L:length(sorted_tol)] -
            sorted_tol[1L:(length(sorted_tol) - 1L)])
        tol_diffs <- pmax(tol_diffs, min_tol_diff)
    }

    if (!is.null(tolerance_list[["tolerance_max"]])) {
        tol_diffs <- pmin(tol_diffs, tolerance_list[["tolerance_max"]])
    }

    ## ------------------------------------

    sigma_xw_sq_both <- (2 / 3) * (y_diffs / tol_diffs)^2

    sum(k_sq_minus_k * sigma_xw_sq_both) / (length(control_index) - 1)
}


#' Simlar to bipartite_match_sd, except for use with
#' non-bipartite matches.
#'
#' For when your result isn't just based
#' on the differences in the y-vector, but those differences are
#' weighed inversely proportional to the difference between the
#' tolerance values.
#' To get on the standard error level, divide the result
#' by \code{sqrt(length(match_list[["treat_index"]]))}
#' This will be conservative in general with
#' repeats.
#' @inheritParams bipartite_match_sd
#' @inheritParams nonbipartite_matches
#' @return a single float, the standard deviation (not standard error)
#'   for the match
#' @author Colman Humphrey
#'
#' @export
nonbipartite_match_sd_scaled <- function(x_mat,
                                         cov_x,
                                         y_vector,
                                         match_list,
                                         tolerance_list = gen_tolerance_list(),
                                         caliper_list = gen_caliper_list(),
                                         weight_vec = NULL,
                                         use_regression = TRUE,
                                         sqrt_mahal = TRUE) {
    ## two variance components

    ## ------------------------------------
    ## V1 easy: usual variance
    ## exists if matches are made without replacement

    y_diffs <- y_vector[match_list[["treat_index"]]] -
        y_vector[match_list[["control_index"]]]
    tol_diffs <-
        tolerance_list[["tolerance_vec"]][match_list[["treat_index"]]] -
        tolerance_list[["tolerance_vec"]][match_list[["control_index"]]]


    var_difference <- if (use_regression) {
                          reg_res <- lm(y_diffs ~ tol_diffs + 0)
                          vcov(reg_res)[1, 1] * length(y_diffs)
                      } else {
                          var(y_diffs / tol_diffs)
                      }

    ## ----------------
    ## V2 harder: the added variance from using repeats

    all_unique_controls <- !any(duplicated(match_list[["control_index"]]))

    if (all_unique_controls) {
        var_repeated <- 0
    } else {
        var_repeated <- gen_nonbipartite_repeated_variance(
            x_mat = x_mat,
            cov_x = cov_x,
            y_vector = y_vector,
            control_index = match_list[["control_index"]],
            tolerance_list = tolerance_list,
            caliper_list = caliper_list,
            weight_vec = weight_vec,
            use_all_controls = use_all_controls,
            sqrt_mahal = sqrt_mahal
        )
    }

    ## ------------------------------------

    var_total <- var_difference + var_repeated

    sqrt(var_total)
}


##' Approximate the standard deviation of the ratio of two vectors
##'
##' The delta method gives us:
##' \deqn{
##'    Var(N/D) \approx \frac{\mu_N^2}{\mu_D^2}
##'                     \bigg[ \frac{\sigma^2_N}{\mu_N^2} +
##'                            \frac{\sigma^2_D}{\mu_D^2} -
##'                            2 \frac{Cov(N, D)}{\mu_N \mu_D}]
##' }
##' Of course it's not necessarily clear we'd have the covariance
##' and not the vectors, but there are indeed some interesting
##' cases.
##'
##' From e.g. Kendall’s Advanced Theory of Statistics or other
##' places with delta method / Taylor expansion
##' @param numerator_mean mean of the numerator
##' @param numerator_sd standard dev of the numerator
##' @param denominator_mean mean of the denominator
##' @param denominator_sd standard dev of the denominator
##' @param covariance covariance of the numerator and denominator
##' @return Returns a single float: the sd estimate
##' @author Colman Humphrey
##'
##' @export
approx_ratio_sd <- function(numerator_mean,
                            numerator_sd,
                            denominator_mean,
                            denominator_sd,
                            covariance) {
    if (numerator_sd <= 0 || denominator_sd <= 0) {
        stop("sd inputs must be positive")
    }

    if (denominator_mean <= 0) {
        stop("denominator should have positive mean")
    }

    if (denominator_mean * 5 < denominator_sd) {
        warning(
            "denominator sd very large compared to mean;",
            "could have negative values?"
        )
    }

    if (abs(covariance) > numerator_sd * denominator_sd) {
        stop(
            "covariance is bigger in absolute value than is possible ",
            "given the two standard deviations"
        )
    }

    mean_ratio <- (numerator_mean / denominator_mean)^2
    n_ratio <- (numerator_sd / numerator_mean)^2
    d_ratio <- (denominator_sd / denominator_mean)^2
    cov_scaled <- covariance / (numerator_mean * denominator_mean)

    sqrt(mean_ratio * (n_ratio + d_ratio - 2 * cov_scaled))
}

##' Quick and dirty method to approximate the sd and mean
##' of the ratio of the differences
##'
##' For non-bipartite work, we work with something like:
##' \deqn{
##'     \frac{y_{\text{treat}} - y_{\text{control}}}
##'          {tol_{\text{treat}} - tol_{\text{control}}}
##' }
##' And not just \eqn{y_{\text{treat}} - y_{\text{control}}} as
##' our "difference": here we say we're caring about the change
##' in \eqn{y} per unit tol.
##' This function estimates the standard deviation of this
##' if you just did random pairings **that obeyed the tolerance
##' rules** with min and max.
##' This function computes approx sd for the **optimal** case.
##' Also yes this is a terrible function name.
##' @param y_vector Vector of outcomes we care about
##' @param tolerance_list Usual tolerance list (see \code{gen_tolerance_list})
##' @param use_regression \code{TRUE} to use regression, \code{FALSE} to use
##'   the ratios directly. Default \code{TRUE}.
##' @param samples How many samples to use. Defaults to a number between 50
##'   and 200, depending on the length of the various vectors
##' @return Returns a single number: the average of the sds of \code{samples}
##'   random runs
##' @author Colman Humphrey
##'
##' @export
y_tolerance_diff_ratio <- function(y_vector,
                                   tolerance_list,
                                   use_regression = TRUE,
                                   samples = NULL) {
    tol_vector <- tolerance_list[["tolerance_vec"]]
    if (length(y_vector) != length(tol_vector)) {
        stop("y_vector and the tolerance vector should have same length")
    }

    if (is.null(samples)) {
        samples <- min(max(50L, ceiling(20000L / length(y_vector))), 200L)
    }

    tol_matches <- lapply(1L:(samples * 2L), function(j) {
        tol_random_sample(tolerance_list)
    })
    pair_count <- unlist(lapply(tol_matches, function(x) {
        length(x[["treat_index"]])
    }))
    keep_pairs <- rank(pair_count, ties.method = "random") > samples
    tol_matches <- tol_matches[keep_pairs]

    ##------------------------------------

    random_ratio_res <- lapply(tol_matches, function(match_list) {
        return(list(
            match_mean = match_estimate_tolerance(
                match_list = match_list,
                y_vector = y_vector,
                tolerance_list = tolerance_list,
                use_regression = use_regression
            ),
            match_sd = nonbipartite_match_sd_scaled(
                x_mat = matrix(NA, 1L, 1L),  # won't use it
                cov_x = matrix(NA, 1L, 1L),  # won't use it
                y_vector = y_vector,
                match_list = match_list,
                tolerance_list = tolerance_list,
                use_regression = use_regression
            )
        ))
    })

    list(
        mean = mean(unlist(
            lapply(random_ratio_res, `[[`, "match_mean")
        )),
        sd = mean(unlist(
            lapply(random_ratio_res, `[[`, "match_sd")
        ))
    )
}


##' Generates a random pairing from a tolerance list in a naive way
##'
##' This function takes in a tolerance list and attempts to generate
##' a "simple" match list by naive random sampling.
##' @param tolerance_list Usual tol, see e.g. \code{gen_tolerance_list}
##' @param prior_pairs For use in recursion - what pairs are already formed
##'   (indeed we have none for the first iteration)
##' @param pairable_units For use in recursion - what units are available
##'   to form pairs (\code{NULL} is all, hence for the first iteration)
##' @param iteration What iteration number are we on
##' @param max_iterations Max iterations
##' @param verbose Boolean - spit out extra messages / warnings?
##' @return "Simple" match list, list with two elements:
##' \itemize{
##'    \item{\code{treat_index}} Treatment units
##'    \item{\code{control_index}} Control units
##' }
##' @author Colman Humphrey
##'
##' @keywords internal
tol_random_sample <- function(tolerance_list,
                              prior_pairs = matrix(NA, nrow = 0L, ncol = 2L),
                              pairable_units = NULL,
                              iteration = 0L,
                              max_iterations = 10L,
                              verbose = FALSE) {
    len_tol <- length(tolerance_list[["tolerance_vec"]])

    pairable_units <- if (is.null(pairable_units)) {
                          seq_len(len_tol)
                      } else {
                          pairable_units
                      }

    if (length(pairable_units) == 0L) {
        return(matrix_to_simple_match(matrix(NA, 0L, 2L)))
    }

    rand_sample <- fixed_sample(pairable_units)

    mid <- length(pairable_units) %/% 2
    pairs <- cbind(
        rand_sample[1L:mid],
        rand_sample[(mid + 1L):(2L * mid)]
    )
    extra_el <- if (length(pairable_units) %% 2 == 1L) {
                    pairable_units[length(pairable_units)]
                } else {
                    NULL
                }
    right_order <- tolerance_list[["tolerance_vec"]][pairs[, 1]] -
        tolerance_list[["tolerance_vec"]][pairs[, 2]] > 0

    if (!all(right_order)) {
        temp <- pairs[!right_order, 1]
        pairs[!right_order, 1] <- pairs[!right_order, 2]
        pairs[!right_order, 2] <- temp
    }

    diffs <- tolerance_list[["tolerance_vec"]][pairs[, 1]] -
        tolerance_list[["tolerance_vec"]][pairs[, 2]]

    ##-------------------------------------

    tol_min <- tolerance_list[["tolerance_min"]]
    tol_max <- if (is.null(tolerance_list[["tolerance_max"]])) {
                   max(diffs) + 1
               } else {
                   tolerance_list[["tolerance_max"]]
               }

    valid_diff <- tol_min < diffs & diffs < tol_max

    if (all(valid_diff)) {
        return(matrix_to_simple_match(rbind(prior_pairs, pairs)))
    }

    good_so_far <- rbind(prior_pairs, pairs[valid_diff, ])

    if (sum(!valid_diff) == 1L) {
        ## can't switch on one
        if (verbose) {
            warning("one unmatchable pair, will be dropped")
        }
        return(matrix_to_simple_match(good_so_far))
    }

    if (iteration == max_iterations) {
        ## best is 0.5, 0.25 is "half" available pairs
        if (nrow(good_so_far) < 0.3 * len_tol || verbose) {
            warning(
                "reached max iterations, but only ",
                nrow(good_so_far),
                " out of a possible ",
                len_tol %/% 2L
            )
        }

        return(matrix_to_simple_match(good_so_far))
    }

    ## -------------------------------------

    if (verbose) {
        message(
            "starting iteration ", iteration + 1L,
            " with ", nrow(good_so_far), " pairs so far",
            " (max: ", len_tol %/% 2, ")"
            )
    }

    return(tol_random_sample(
        tolerance_list,
        prior_pairs = good_so_far,
        pairable_units = c(pairs[!valid_diff, , drop = FALSE]),
        iteration = iteration + 1L,
        max_iterations = max_iterations,
        verbose = verbose
    ))
}

##' Simple function for \code{tol_random_sample} to convert matrices
##' of pairs into list
##'
##' @param pairs_matrix Matrix with two columns, one for the treatment
##'   indices, one for the control
##' @param treat_first Boolean; is the first column the treatment column?
##' @return See \code{tol_random_sample}
##' @author Colman Humphrey
##'
##' @keywords internal
matrix_to_simple_match <- function(pairs_matrix,
                                   treat_first = TRUE) {
    if (ncol(pairs_matrix) != 2L) {
        stop("`pairs_matrix` must have exactly two columns")
    }
    treat_ind <- ifelse(treat_first, 1L, 2L)
    control_ind <- 3L - treat_ind
    return(list(
        treat_index = pairs_matrix[, treat_ind],
        control_index = pairs_matrix[, control_ind]
    ))
}
