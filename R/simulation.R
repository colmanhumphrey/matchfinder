#' Generates input needed for a simulation run
#'
#' @param n_rows How many rows to generate.
#' @param n_cols How many columns to use.
#' @param x_generator Function that takes number of rows and columns
#'   and produces a data matrix of that dimension.
#' @param treat_prob_generator Function that takes a matrix and produces
#'   treatment probabilities.
#' @param mean_generator Function that takes a matrix and produces
#'   an expected value for each row.
#' @param error_generator Function that accepts a number of rows and
#'   generates an error vector (e.g. generates normal noise). Should
#'   generate zero-mean data.
#' @return List:
#' \describe{
#'   \item{\code{x_mat}}{Data Matrix}
#'   \item{\code{treat_vec}}{Treatment Vector}
#'   \item{\code{y_vec}}{Output Vector}
#' }
#' @author Colman Humphrey
#'
#' @export
generate_simulation_input <- function(n_rows = 500L,
                                      n_cols = 5L,
                                      x_generator =
                                          default_x_generator,
                                      treat_prob_generator =
                                          example_treat_prob_generator,
                                      mean_generator =
                                          example_mean_generator,
                                      error_generator =
                                          default_error_generator) {
    x_mat <- x_generator(n_rows = n_rows,
                         n_cols = n_cols)

    treat_prob <- treat_prob_generator(x_mat)
    treat_vec <- rbinom(length(treat_prob), 1, treat_prob)

    mu_pre_treat <- mean_generator(x_mat)

    list(
        x_mat = x_mat,
        treat_vec = treat_vec,
        y_vec = mu_pre_treat + treat_vec + error_generator(n_rows = n_rows)
    )
}


#' Generates the four different treatment functions we use
#' for simulations in the paper.
#'
#' @param target_mean Desired mean for each of the treatment probabilities,
#'   should be in (0, 1).
#' @return List of four functions that take in a matrix and
#'   output a vector of probabilities:
#' \describe{
#'   \item{\code{constant_treat_prob}}{Just repeats \code{target_mean}}
#'   \item{\code{logistic_treat_prob}}{Computes a logistic style relationship
#'   between the data matrix and the output vector}
#'   \item{\code{sparse_treat_prob}}{Based on the sign of the first data matrix
#'   column; gives a number a bit above the mean and a bit below.}
#'   \item{\code{sparse_nonlin_treat_prob}}{Cubic logit function of the
#'   first vector}
#' }
#' @author Colman Humphrey
#'
#' @export
paper_treatment_functions <- function(target_mean = 0.425) {
    stopifnot(length(target_mean) == 1 &&
              0 < target_mean && target_mean < 1)

    constant_treat_prob <- function(x_mat) {
        rep(target_mean, nrow(x_mat))
    }

    logistic_treat_prob <- function(x_mat) {
        coef_vec <- rnorm(ncol(x_mat), 0, 0.2)
        lin_vec <- c(x_mat %*% coef_vec)

        target_mean_expit(target_mean = target_mean,
                          linear_vector = lin_vec)
    }

    sparse_treat_prob <- function(x_mat) {
        median_adjusted <- x_mat[, 1] - median(x_mat[, 1])

        min_gap <- min(target_mean, 1 - target_mean) / 2

        ifelse(sign(median_adjusted) > 0,
               target_mean + min_gap,
               target_mean - min_gap)
    }

    sparse_nonlin_treat_prob <- function(x_mat) {
        mean_adj <- (x_mat[, 1] - mean(x_mat[, 1])) / (sd(x_mat[, 1]) * 2) + 1
        numer <- mean_adj^3 - mean_adj^2

        target_mean_expit(target_mean = target_mean,
                          linear_vector = numer)
    }

    ##------------------------------------

    list(
        constant_treat_prob = constant_treat_prob,
        logistic_treat_prob = logistic_treat_prob,
        sparse_treat_prob = sparse_treat_prob,
        sparse_nonlin_treat_prob = sparse_nonlin_treat_prob
    )
}


#' Generates the four different mean functions we use
#' for simulations in the paper.
#'
#' We make them all mean zero for no real reason, we don't really care
#' what the mean is.
#' @return List of four functions that take in a matrix and
#'   output a vector of means:
#' \describe{
#'   \item{\code{constant_mu}}{Vector of zeroes}
#'   \item{\code{linear_mu}}{The matrix multiplied by a random vector.}
#'   \item{\code{sign_mu}}{Literally \code{sign(x_mat[, 1])}}
#'   \item{\code{non_linear_mu}}{Messy non-linear function of the data matrix}
#' }
#' @author Colman Humphrey
#'
#' @export
paper_mean_functions <- function() {
    constant_mu <- function(x_mat) {
        rep(0, nrow(x_mat))
    }

    linear_mu <- function(x_mat) {
        coef_vec <- rnorm(ncol(x_mat), 0, 1)
        lin_vec <- c(x_mat %*% coef_vec)

        lin_vec - mean(lin_vec)
    }

    sign_mu <- function(x_mat) {
        median_adjusted <- x_mat[, 1] - median(x_mat[, 1])

        sign(median_adjusted) - mean(sign(median_adjusted))
    }

    non_linear_mu <- function(x_mat) {
        coef_vec <- rnorm(ncol(x_mat), 0, 0.5)
        min_abs <- apply(abs(x_mat), 1, min)

        end_col <- ncol(x_mat)

        lin_vec <- cos(x_mat[, 1]) +
            sin(c(x_mat %*% coef_vec)) * sign(x_mat[, end_col]) -
            pmin(tan(pi / 3 - min_abs * pi), 10)

        lin_vec - mean(lin_vec)
    }

    ##------------------------------------

    list(
        constant_mu = constant_mu,
        linear_mu = linear_mu,
        sign_mu = sign_mu,
        non_linear_mu = non_linear_mu
    )
}


#' Generates a function that generates a vector of sink lengths
#'
#' This function returns a function that accepts a treatment vector
#' and generates a vector of numbers to be used as sink counts.
#' @param start_frac Smaller fraction of units to use as sink number, default 0.
#' @param end_frac Larger fraction of units to use as sink number, default 0.8.
#' @param length_out How many sink values we want, default 9.
#' @return Function that accepts a \code{treat_vec}
#'   and returns a vector of numbers.
#' @author Colman Humphrey
#'
#' @export
n_sink_generator <- function(start_frac = 0,
                             end_frac = 0.8,
                             length_out = 9) {
    stopifnot(length(start_frac) == 1L &&
              length(end_frac) == 1L &&
              length(length_out) == 1L)
    stopifnot(0 <= start_frac && start_frac <= 1)
    stopifnot(0 <= end_frac && end_frac <= 1)
    stopifnot(length_out >= 1L)

    if (length_out == 1) {
        stopifnot(start_frac == end_frac)
    } else {
        stopifnot(start_frac < end_frac)
    }


    function(treat_vec) {
        num_treat <- sum(treat_vec)

        floor(seq(from = floor(start_frac * num_treat),
                  to = floor(end_frac * num_treat),
                  length.out = length_out))
    }
}


#' Takes in functions to generate simulation data, and computes
#' simulation results for our method and some competitors
#'
#' @inheritParams generate_simulation_input
#' @inheritParams all_bipartite_matches
#' @param n_sink_gen Default \code{n_sink_generator}, and you'll
#'   probably want to use that: this argument should be a function that
#'   accepts a \code{treat_vec} and produces a vector of sink numbers.
#' @param num_weight_vectors How many weight vectors to generate.
#' @param silent Default \code{!interactive()}, if you want to
#'   suppress messages.
#' @return Returns a named list:
#' \describe{
#'     \item{\code{naive_est}}{Just a number: mean difference between all
#'   treated units and all control}
#'     \item{\code{propensity_results}}{List of lists: each with
#'   \code{n_sinks} and the \code{est}}
#'     \item{\code{mahal_results}}{Same as above}
#'     \item{\code{weighted_results}}{List of lists: each with
#'   \code{n_sinks}, the raw brier score, the permutation brier score,
#'   and the \code{est}}
#' }
#' @author Colman Humphrey
#'
#' @export
compute_sim_result <- function(x_generator = default_x_generator,
                               treat_prob_generator,
                               mean_generator,
                               error_generator = default_error_generator,
                               n_sink_gen = n_sink_generator(),
                               match_method = "with_replacement",
                               n_rows = 500L,
                               n_cols = 5L,
                               num_weight_vectors = 100L,
                               silent = !interactive()) {

    sim_data <- generate_simulation_input(n_rows = n_rows,
                                          n_cols = n_cols,
                                          x_generator = x_generator,
                                          treat_prob_generator = treat_prob_generator,
                                          mean_generator = mean_generator,
                                          error_generator = error_generator)

    x_mat <- sim_data[["x_mat"]]
    y_vector <- sim_data[["y_vec"]]
    treat_vec <- sim_data[["treat_vec"]]

    rm(sim_data)

    n_sinks <- n_sink_gen(treat_vec)
    match_list_est_func <- (function(y_vector, treat_vec) {
        function(match_list) {
            match_estimate(match_list = match_list,
                           y_vector = y_vector,
                           treat_vec = treat_vec)
        }
    })(y_vector, treat_vec)

    list_est_func <- (function(n_sinks) {
        function(match_lists) {
            Map(function(n_sink, match_list) {
                list(
                    n_sinks = n_sink,
                    est = match_list_est_func(match_list)
                )
            },
            n_sinks,
            match_lists)
        }
    })(n_sinks)

    naive_est <- mean(y_vector[treat_vec == 1]) - mean(y_vector[treat_vec == 0])

    ##------------------------------------

    if (!silent) {
        message("propensity matches")
    }

    propensity_matches <- propensity_bipartite_matches(
        x_mat = x_mat,
        treat_vec = treat_vec,
        match_method = match_method,
        propensity_list = gen_propensity_list(
            propensity_function = propensity_score_linear,
            oos_propensity = FALSE
        ),
        n_sinks = n_sinks)

    propensity_ests <- list_est_func(propensity_matches)

    ##------------------------------------

    if (!silent) {
        message("mahal matches")
    }

    mahal_matches <- all_bipartite_matches(
        x_mat = x_mat,
        cov_x = covariance_with_ranks(x_mat),
        weight_list = list(rep(1 / n_cols, times = n_cols)),
        treat_vec = treat_vec,
        match_method = match_method,
        n_sinks = n_sinks)[[1]]

    mahal_ests <- list_est_func(mahal_matches)

    ##------------------------------------

    if (!silent) {
        message("all weighted matches")
    }

    weight_list <- generate_random_weights(
        prior_weights = rep(1 / n_cols, times = n_cols),
        number_vectors = num_weight_vectors,
        minimum_weights = rep(1 / (3 * n_cols), times = n_cols))

    sink_brier_matches <- sink_brier_bipartite_matches(
        x_mat = x_mat,
        cov_x = covariance_with_ranks(x_mat),
        weight_list = weight_list,
        treat_vec = treat_vec,
        match_method = match_method,
        n_sinks = n_sinks,
        silent = silent)

    permutation_results <- permutation_bipartite_matches(
        matches_by_sinks = sink_brier_matches[["matches_by_sinks"]],
        briers_by_sinks = sink_brier_matches[["briers_by_sinks"]],
        x_mat = x_mat,
        n_sinks = n_sinks,
        silent = silent)

    weighted_results <- lapply(
        permutation_results[["best_matches"]],
        function(match_results) {
            list(
                n_sinks = match_results[["n_sinks"]],
                raw_brier = match_results[["raw_brier"]],
                permutation_brier = match_results[["permutation_brier"]],
                est = match_list_est_func(match_results[["match_list"]]))
        })

    ##------------------------------------

    list(
        naive_est = naive_est,
        propensity_results = propensity_ests,
        mahal_results = mahal_ests,
        weighted_results = weighted_results)
}


#' Reshapes a list of simulations to a nice dataframe
#'
#' @param list_of_sims List of results from \code{compute_sim_results}.
#' @param treat_model_name Name of the treatment model.
#' @param mu_model_name Name of the mean generation model.
#' @param n_rows How many rows were used.
#' @param n_cols How many columns were used.
#' @param num_weight_vectors How many weight vectors were used
#' @return Data frame spreading results...
#' @author Colman Humphrey
#'
#' @export
reshape_list_of_sims <- function(list_of_sims,
                                 treat_model_name,
                                 mu_model_name,
                                 n_rows,
                                 n_cols,
                                 num_weight_vectors) {
    list_of_sims <- lapply(1L:length(list_of_sims), function(j) {
        if (is.null(list_of_sims[[j]][["id"]])) {
            list_of_sims[[j]][["id"]] <- j
        }
        return(list_of_sims[[j]])
    })

    do.call(rbind, lapply(list_of_sims, function(sim_res) {
        do.call(rbind, lapply(1L:length(sim_res[["weighted_results"]]), function(j) {
            data.frame(
                id = sim_res[["id"]],
                treat_model = treat_model_name,
                mu_model = mu_model_name,
                p_brier =
                    sim_res[["weighted_results"]][[j]][["permutation_brier"]],
                raw_brier =
                    sim_res[["weighted_results"]][[j]][["raw_brier"]],
                n_rows = n_rows,
                n_cols = n_cols,
                num_weight_vectors = num_weight_vectors,
                n_sinks = sim_res[["weighted_results"]][[j]][["n_sinks"]],
                naive_est = sim_res[["naive_est"]],
                propensity_est =
                    sim_res[["propensity_results"]][[j]][["est"]],
                mahal_est = sim_res[["mahal_results"]][[j]][["est"]],
                weighted_est = sim_res[["weighted_results"]][[j]][["est"]]
            )
        }))
    }))
}


#' Reshapes a list of simulations to a nice dataframe, by p-cut
#'
#' @inheritParams reshape_list_of_sims
#' @param p_cut Vector of one-sided permution p-values to use.
#' @return Data frame spreading results...
#' @author Colman Humphrey
#'
#' @export
reshape_p_cut_list <- function(list_of_sims,
                               treat_model_name,
                               mu_model_name,
                               p_cut,
                               n_rows,
                               n_cols,
                               num_weight_vectors) {
    list_of_sims <- lapply(1L:length(list_of_sims), function(j) {
        if (is.null(list_of_sims[[j]][["id"]])) {
            list_of_sims[[j]][["id"]] <- j
        }
        return(list_of_sims[[j]])
    })

    do.call(rbind, lapply(p_cut, function(p_cut_val) {
        do.call(rbind, lapply(list_of_sims, function(par_res) {
            p_briers <- unlist(lapply(
                par_res[["weighted_results"]], function(x) {
                    x[["permutation_brier"]]
                }))
            if (all(p_briers > p_cut_val)) {
                given_cut_ind <- which.min(p_briers)
            } else {
                given_cut_ind <- which(p_briers <= p_cut_val)[1]
            }

            data.frame(
                id = par_res[["id"]],
                treat_model = treat_model_name,
                mu_model = mu_model_name,
                p_cut = p_cut_val,
                p_brier =
                    par_res[["weighted_results"]][[given_cut_ind]][["permutation_brier"]],
                raw_brier =
                    par_res[["weighted_results"]][[given_cut_ind]][["raw_brier"]],
                n_rows = n_rows,
                n_cols = n_cols,
                num_weight_vectors = num_weight_vectors,
                n_sinks = par_res[["weighted_results"]][[given_cut_ind]][["n_sinks"]],
                naive_est = par_res[["naive_est"]],
                propensity_est =
                    par_res[["propensity_results"]][[given_cut_ind]][["est"]],
                mahal_est = par_res[["mahal_results"]][[given_cut_ind]][["est"]],
                weighted_est = par_res[["weighted_results"]][[given_cut_ind]][["est"]]
            )
        }))
    }))
}


#' Run \code{compute_sim_result} in parallel using \code{parallel::mclapply}
#'
#' @inheritParams compute_sim_result
#' @param num_cores How many cores to use.
#' @param iterations How many iterations to do.
#' @author Colman Humphrey
#'
#' @export
parallel_sim <- function(x_generator = default_x_generator,
                         treat_prob_generator,
                         mean_generator,
                         error_generator = default_error_generator,
                         n_sink_gen = n_sink_generator(),
                         match_method = "with_replacement",
                         n_rows = 500L,
                         n_cols = 5L,
                         num_weight_vectors = 100L,
                         num_cores = parallel::detectCores() - 1,
                         iterations = 100L,
                         names_list = NULL,
                         silent = !interactive()) {
    sims <- parallel::mclapply(1L:iterations, function(j) {
        if (!silent) {
            print(paste0("iteration ", j, "/", iterations))
        }
        compute_sim_result(x_generator = x_generator,
                           treat_prob_generator = treat_prob_generator,
                           mean_generator = mean_generator,
                           error_generator = error_generator,
                           n_sink_gen = n_sink_gen,
                           match_method = match_method,
                           n_rows = n_rows,
                           n_cols = n_cols,
                           num_weight_vectors = num_weight_vectors,
                           silent = silent)
    }, mc.cores = num_cores)
}
