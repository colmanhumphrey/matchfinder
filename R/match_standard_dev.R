#' Calculates the standard deviation for a match, on the outcome
#'
#' @inheritParams all_bipartite_matches
#' @param y_vector Outcome vector (not used in match generation).x
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
        stop("don't use both `caliper_list` and `propensity_list`: ",
             " If you do want both, create the combined caliper separately")
    }
    ##------------------------------------
    ## two variance components

    ##------------------------------------
    ## V1 easy: usual variance
    ## exists if matches are made without replacement

    var_difference <- var(y_vector[match_list[["treat_index"]]] -
                          y_vector[match_list[["control_index"]]])

    ##----------------
    ## V2 harder: the added variance from using repeats

    all_unique_controls <- !any(duplicated(match_list[["control_index"]]))

    if (all_unique_controls) {
        var_repeated <- 0
    } else {
        if (!is.null(propensity_list)) {
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
            caliper_list = gen_caliper_list(
                caliper_vec = prop_score,
                caliper_max = sd(prop_score) * propensity_list[["caliper_sd_mult"]],
                continuous_mult = propensity_list[["continuous_mult"]])
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
            sqrt_mahal = sqrt_mahal)
    }

    ##------------------------------------

    var_total <- var_difference + var_repeated

    sqrt(var_total)
}


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
    k_sq_minus_k = count_frame[["count"]]^2 - count_frame[["count"]]

    unique_control_index <- as.numeric(as.character(
        count_frame[["control_index"]]))

    use_all_controls <- ifelse(is.null(use_all_controls),
                               TRUE, use_all_controls)

    ##------------------------------------

    control_within_all_index <- match(unique_control_index,
                                      all_control_index)

    if (use_all_controls) {
        control_controls <- all_control_index
    } else {
        control_controls <- unique_control_index
    }

    blocked_ind <- match(unique_control_index,
                         control_controls)
    partial_index <- list(unique_control_index,
                          control_controls)

    control_dist_mat <- weighted_mahal(x_mat = x_mat,
                                       cov_x = cov_x,
                                       weight_vec = weight_vec,
                                       sqrt_mahal = sqrt_mahal,
                                       partial_index = partial_index)

    if (!is.null(caliper_list)) {
        caliper_full <- create_caliper(caliper_list)
        control_dist_mat <- control_dist_mat +
            caliper_full[cbind(partial_index[[1]],
                               partial_index[[2]])]
    }

    min_control_match <- min_blocked_rank(control_dist_mat,
                                          blocked_ind)

    control_match <- control_controls[min_control_match]

    ##------------------------------------

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
                                               tolerance_list = gen_tolerance_list(),
                                               caliper_list = gen_caliper_list(),
                                               weight_vec = NULL,
                                               use_all_controls = TRUE,
                                               sqrt_mahal = TRUE) {
    count_frame <- data.frame(table(control_index))
    names(count_frame) <- c("control_index", "count")

    ## the second term should actually be K_{sq,i},
    ## but when using just one control, we have K_i = K_{sq,i}
    k_sq_minus_k = count_frame[["count"]]^2 - count_frame[["count"]]

    unique_control_index <- as.numeric(as.character(
        count_frame[["control_index"]]))

    ##------------------------------------

    partial_index <- list(unique_control_index,
                          1L:nrow(x_mat))

    control_dist_mat <- weighted_mahal(x_mat = x_mat,
                                       cov_x = cov_x,
                                       weight_vec = weight_vec,
                                       sqrt_mahal = sqrt_mahal,
                                       partial_index = partial_index)

    if (!is.null(caliper_list)) {
        caliper_full <- create_caliper(caliper_list)
        control_dist_mat <- control_dist_mat +
            caliper_full[cbind(partial_index[[1]],
                               partial_index[[2]])]
    }

    min_control_match <- min_blocked_rank(control_dist_mat,
                                          blocked_ind = unique_control_index)

    ##----------------

    tolerance_match <- near_given_match(
        tolerance_list[["tolerance_vec"]],
        given_index = unique_control_index)

    ##----------------

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

    ##------------------------------------

    sigma_xw_sq_both <- (2/3) * (y_diffs / tol_diffs)^2

    sum(k_sq_minus_k * sigma_xw_sq_both) / (length(control_index) - 1)
}


#' Simlar to bipartite_match_sd, except for use with
#' non-bipartite matches.
#'
#' For when your result isn't just based
#' on the differences in the y-vector, but those differences are
#' weighed inversely proportional to the difference between the
#' tolerance values
#' @inheritParams bipartite_match_sd
#' @param tolerance_vec the tolerance used to form the
#'   non-bipartite matches
#' @param tolerance_min min value we must obey in tol diffs
#' @param tolerance_max max value we must obey in tol diffs
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
                                         sqrt_mahal = TRUE) {
    ## two variance components

    ##------------------------------------
    ## V1 easy: usual variance
    ## exists if matches are made without replacement

    y_diffs <- y_vector[match_list[["treat_index"]]] -
        y_vector[match_list[["control_index"]]]
    tol_diffs <- tolerance_list[["tolerance_vec"]][match_list[["treat_index"]]] -
        tolerance_list[["tolerance_vec"]][match_list[["control_index"]]]

    var_difference <- var(y_diffs / tol_diffs)

    ##----------------
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
            sqrt_mahal = sqrt_mahal)
    }

    ##------------------------------------

    var_total <- var_difference + var_repeated

    sqrt(var_total)
}
