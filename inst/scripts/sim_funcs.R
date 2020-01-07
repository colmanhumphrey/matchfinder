library(MASS)
library(matchfinder)

constant_treat_prob <- function(x_mat) {
    rep(0.425, nrow(x_mat))
}
logistic_treat_prob <- function(x_mat,
                                coef_vec,
                                alpha) {
    matchfinder:::expit(x_mat %*% coef_vec + alpha)
}
sparse_treat_prob <- function(x_mat) {
    ifelse(sign(x_mat[, 1]) > 0, 0.7, 0.2)
}
sparse_nonlin_treat_prob <- function(x_mat) {
    numer <- (x_mat[, 1]^3 - 2 * x_mat[, 1]^2) / 10
    matchfinder:::expit(numer)
}

constant_mu <- function(x_mat) {
    rep(0, nrow(x_mat))
}
linear_mu <- function(x_mat, coef_vec) {
    x_mat %*% coef_vec
}
sign_mu <- function(x_mat) {
    sign(x_mat[, 1])
}
non_linear_mu <- function(x_mat, coef_vec) {
    min_abs <- apply(abs(x_mat), 1, min)
    cos(x_mat[, 1]) + sin(x_mat %*% coef_vec) * sign(x_mat[, 3]) -
        tan(pi / 3 - min_abs * pi)
}


#' Generating simulation data according to one of our setups
#'
#' @param treat_model
#' @param mu_model
#' @param n_rows
#' @param n_cols
#' @param x_norm
#' @return
#' @author Colman Humphrey
sim_input_list <- function(treat_model = c("a", "b", "c", "d"),
                           mu_model = c("a", "b", "c", "d"),
                           n_rows = 500,
                           n_cols = 5,
                           x_norm = TRUE) {
    treat_model <- match.arg(treat_model)
    mu_model <- match.arg(mu_model)

    if (x_norm) {
        sig_mat_pre = matrix(runif((n_cols - 1) * n_cols, 0, 0.1),
                             n_cols, n_cols - 1)
        sig_mat_cov = sig_mat_pre %*% t(sig_mat_pre) +
            diag(x = runif(n_cols, 0, 0.3))
        sig_mat_cor = cov2cor(sig_mat_cov)
        x_mat = MASS::mvrnorm(n = n_rows,
                              mu = rep(0, n_cols),
                              Sigma = sig_mat_cor)
    } else {
        x_mat <- matrix(runif(n_rows * n_cols, -1, 1),
                        nrow = n_rows)
    }

    ##------------------------------------

    if (treat_model == "a") {
        treat_prob <- constant_treat_prob(x_mat)
    }
    if (treat_model == "b") {
        coef_vec <- rnorm(n_cols, 0, 0.2)

        mean_logit <- 0
        while(abs(mean_logit - 0.425) > 0.02) {
            alp_val <- rnorm(1, -0.2, 0.3)
            logit_treat_prob <- logistic_treat_prob(x_mat,
                                                    coef_vec,
                                                    alp_val)
            mean_logit <- mean(logit_treat_prob)
        }
        treat_prob <- logit_treat_prob
    }
    if (treat_model == "c") {
        treat_prob <- sparse_treat_prob(x_mat)
    }
    if (treat_model == "d") {
        treat_prob <- sparse_nonlin_treat_prob(x_mat)
    }

    treat_vec <- rbinom(length(treat_prob), 1, treat_prob)

    if (mu_model == "a") {
        mu_pre_treat <- constant_mu(x_mat)
    }
    if (mu_model == "b") {
        mu_pre_treat <- linear_mu(x_mat, rnorm(n_cols, 0, 0.2))
    }
    if (mu_model == "c") {
        mu_pre_treat <- sign_mu(x_mat)
    }
    if (mu_model == "d") {
        mu_pre_treat <- non_linear_mu(x_mat, rnorm(n_cols, 0, 0.2))
    }

    list(
        x_mat = x_mat,
        treat_vec = treat_vec,
        y_vec = mu_pre_treat + treat_vec + rnorm(n_rows, 0, 1)
    )
}

ml_est <- function(ml, treat_vec, y_vector) {
    stopifnot(all(treat_vec[ml[["treat_index"]]] == 1L))
    stopifnot(all(treat_vec[ml[["control_index"]]] == 0L))

    mean(y_vector[ml[["treat_index"]]] -
         y_vector[ml[["control_index"]]])
}

n_sink_func_func <- function(start_frac = 0,
                             end_frac = 0.8,
                             length_out = 9) {
    function(treat_vec) {
        num_treat <- sum(treat_vec)
        floor(seq(from = floor(start_frac * num_treat),
                  to = floor(end_frac * num_treat),
                  length.out = length_out))
    }
}

compute_sim_result <- function(treat_model = c("a", "b", "c", "d"),
                               mu_model = c("a", "b", "c", "d"),
                               n_sink_gen = n_sink_func_func(),
                               n_rows = 500,
                               n_cols = 5,
                               num_weight_vectors = 20,
                               match_method = "with_replacement",
                               x_norm = TRUE,
                               silent = TRUE) {
    treat_model <- match.arg(treat_model)
    mu_model <- match.arg(mu_model)

    sim_data <- sim_input_list(treat_model = treat_model,
                               mu_model = mu_model,
                               n_rows = n_rows,
                               n_cols = n_cols,
                               x_norm = x_norm)

    x_mat <- sim_data[["x_mat"]]
    y_vector <- sim_data[["y_vec"]]
    treat_vec <- sim_data[["treat_vec"]]

    n_sinks <- n_sink_gen(treat_vec)

    naive_sd <- sqrt(var(y_vector[treat_vec == 1]) / sum(treat_vec == 1) +
                     var(y_vector[treat_vec == 0]) / sum(treat_vec == 0))
    naive_est <- mean(y_vector[treat_vec == 1]) -
        mean(y_vector[treat_vec == 0])

    ml_est_func <- (function(treat_vec, y_vector) {
        function(ml){
            ml_est(ml, treat_vec, y_vector)
        }
    }) (treat_vec, y_vector)

    ##----------------

    if (!silent) {
        message("propensity matches")
    }

    propensity_matches <- propensity_bipartite_matches(
        x_mat,
        treat_vec = treat_vec,
        match_method = match_method,
        propensity_list = gen_propensity_list(propensity_score_linear,
                                              oos_propensity = FALSE),
        n_sinks = n_sinks)

    propensity_ests <- lapply(1L:length(n_sinks), function(j) {
        list(
            n_sinks = n_sinks[j],
            est = ml_est_func(propensity_matches[[j]])
        )
    })

    ##----------------
    ## mahalanobis

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

    mahal_ests <- lapply(1L:length(n_sinks), function(j) {
        list(
            n_sinks = n_sinks[j],
            est = ml_est_func(mahal_matches[[j]])
        )
    })

    ##----------------

    if (!silent) {
        message("all weighted matches")
    }

    weight_list <- generate_random_weights(
        prior_weights = rep(1 / n_cols, times = n_cols),
        number_vectors = num_weight_vectors,
        minimum_weights = rep(1 / (3 * n_cols), times = n_cols))

    all_matches <- all_bipartite_matches(
        x_mat = x_mat,
        cov_x = covariance_with_ranks(x_mat),
        weight_list = weight_list,
        treat_vec = treat_vec,
        match_method = match_method,
        n_sinks = n_sinks)

    all_by_sinks <- lapply(n_sinks, function(x) {
        lapply(all_matches, function(y) {
            y[[as.character(x)]]
        })
    })

    if (!silent) {
        message("getting briers")
    }

    briers_by_sinks <- lapply(all_by_sinks, function(x){
        if (!silent) {
            print(x[[1]]["num_sinks"])
        }
        unlist(lapply(x, function(y) {
            brier_score_cv(x_mat,
                           y)
        }))
    })

    best_brier_inds <- lapply(briers_by_sinks, function(x) {
        which(rank(x, ties.method = "random") == length(x))
    })

    if (!silent) {
        message("getting permutation briers")
    }

    permutation_briers <- lapply(1L:length(n_sinks), function(j) {
        if (!silent) {
            n_sinks[j]
        }
        best_brier_ind <- best_brier_inds[[j]]
        permutation_brier(x_mat,
                          match_list = all_by_sinks[[j]][[best_brier_ind]])
    })

    permutation_brier_score <- lapply(1L:length(n_sinks), function(j) {
        permutation_vec <- permutation_briers[[j]]
        unlist(lapply(briers_by_sinks[[j]], function(x) {
            mean(x <= permutation_vec)
        }))
    })

    ## now that we're doing one-sided brier,
    ## the lowest value will just be the best
    ## so will highest brier

    weighted_results <- lapply(1L:length(n_sinks), function(j) {
        best_brier_ind <- best_brier_inds[[j]]

        stopifnot(permutation_brier_score[[j]][best_brier_ind] ==
                  min(permutation_brier_score[[j]]))

        list(
            n_sinks = n_sinks[j],
            raw_brier = briers_by_sinks[[j]][best_brier_ind],
            permutation_brier = permutation_brier_score[[j]][best_brier_ind],
            est = ml_est_func(all_by_sinks[[j]][[best_brier_ind]])
        )
    })

    ##------------------------------------

    list(
        naive_est = naive_est,
        propensity_results = propensity_ests,
        mahal_results = mahal_ests,
        weighted_results = weighted_results
    )
}

parallel_sim <- function(treat_model = c("a", "b", "c", "d"),
                         mu_model = c("a", "b", "c", "d"),
                         n_sink_gen = n_sink_func_func(),
                         n_rows = 1000,
                         n_cols = 5,
                         num_weight_vectors = 100,
                         match_method = "with_replacement",
                         num_cores = parallel::detectCores() - 1,
                         iterations = 100,
                         x_norm = TRUE,
                         silent = FALSE) {
    treat_model <- match.arg(treat_model)
    mu_model <- match.arg(mu_model)

    parallel::mclapply(1L:iterations, function(j) {
        if (!silent) {
            print(j)
        }
        compute_sim_result(treat_model = treat_model,
                           mu_model = mu_model,
                           n_sink_gen = n_sink_gen,
                           n_rows = n_rows,
                           n_cols = n_cols,
                           num_weight_vectors = num_weight_vectors,
                           match_method = match_method,
                           x_norm = x_norm,
                           silent = TRUE)
    }, mc.cores = num_cores)
}
