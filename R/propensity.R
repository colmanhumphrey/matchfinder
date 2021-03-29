#' Takes in elements needed for propensity work, checks input and
#' builds a named list.
#'
#' @param propensity_function A function that accepts a list
#'   with four elements: \code{x_train}, \code{x_test},
#'   \code{y_train}, \code{y_test}, and forms predictions
#'   using \code{x_test} (I guess \code{y_test} isn't used yet)
#' @param oos_propensity Logical, do you want to predict out of sample
#'   for the propensity score? Most people don't, and indeed \code{FALSE}
#'   is the default.
#' @param n_folds Default NULL; how many folds you want if using
#'   out of sample propensity.
#' @return Named list, same names as input params.
#' @author Colman Humphrey
#'
#' @export
gen_propensity_list <- function(propensity_function = propensity_score_xgb(),
                                oos_propensity = FALSE,
                                n_folds = NULL) {
    stopifnot(is_tf(oos_propensity))

    if (missing(oos_propensity) && !is.null(n_folds)) {
        oos_propensity <- TRUE
    }

    if (!oos_propensity) {
        if (!is.null(n_folds)) {
            stop("n_folds shouldn't be set if not using out of sample")
        }
    } else {
        n_folds <- ifelse(is.null(n_folds), 5L, n_folds)
        stopifnot(is.numeric(n_folds) && length(n_folds) == 1L &&
            !is.na(n_folds))
    }

    list(
        propensity_function = propensity_function,
        oos_propensity = oos_propensity,
        n_folds = n_folds
    )
}


#' Generates the propensity parameters used for using propensity-based calipers
#'
#' We use this for input for \code{all_propensity_caliper_matches} (and
#' likely a nonbipartite version soon).
#' @inheritParams gen_propensity_list
#' @param caliper_sd_mult We'll set the maximum gap between units
#'   as \code{sd(propensity_score) * k}, where this parameter is the
#'   value k. Default 0.6.
#' @param continuous_mult See e.g. \code{gen_caliper_list}: instead of
#'   blocking matches that are "too far apart" on the caliper, we'll
#'   add a penalty for going above.
#' @return list with names equal to all input params
#' @export
match_propensity_list <- function(propensity_function = propensity_score_xgb(),
                                  oos_propensity = FALSE,
                                  n_folds = NULL,
                                  caliper_sd_mult = 0.6,
                                  continuous_mult = 100) {
    if (is.null(propensity_function)) {
        return(NULL)
    }

    plain_prop_list <- gen_propensity_list(
        propensity_function = propensity_function,
        oos_propensity = oos_propensity,
        n_folds = n_folds
    )

    list(
        propensity_function = plain_prop_list[["propensity_function"]],
        oos_propensity = plain_prop_list[["oos_propensity"]],
        n_folds = plain_prop_list[["n_folds"]],
        caliper_sd_mult = caliper_sd_mult,
        continuous_mult = continuous_mult
    )
}


#' Propensity match for (default) bipartite
#' @inheritParams all_bipartite_matches
#' @param propensity_list See \code{gen_propensity_list}
#'
#' @export
propensity_bipartite_matches <- function(x_mat,
                                         treat_vec,
                                         match_method = c(
                                             "with_replacement",
                                             "optimal",
                                             "greedy"
                                         ),
                                         propensity_list =
                                             gen_propensity_list(),
                                         n_sinks = 0,
                                         caliper_list = gen_caliper_list(),
                                         sqrt_mahal = TRUE,
                                         tol_val = NULL) {
    ## in case of logical
    treat_vec <- treat_vec * 1L

    ## generate propensity score
    prop_score <- propensity_score(
        x_mat = x_mat,
        treat_vec = treat_vec,
        propensity_list = propensity_list
    )
    prop_dist_mat <- abs(outer(
        prop_score[treat_vec == 1],
        prop_score[treat_vec == 0],
        "-"
    ))

    if (!is.null(caliper_list)) {
        prop_dist_mat <- prop_dist_mat + create_caliper(caliper_list,
            treat_vec = treat_vec
        )
    }

    bipartite_matches(
        dist_mat = prop_dist_mat,
        treat_vec = treat_vec,
        match_method = match_method,
        n_sinks = n_sinks,
        tol_val = tol_val
    )
}

#' Propensity match for nbp
#' @inheritParams all_nonbipartite_matches
#' @param propensity_list See \code{gen_propensity_list}
propensity_nonbipartite_matches <- function(x_mat,
                                            tolerance_list =
                                                gen_tolerance_list(),
                                            propensity_list =
                                                gen_propensity_list(),
                                            match_method = c(
                                                "with_replacement",
                                                "optimal",
                                                "greedy"
                                            ),
                                            n_sinks = 0,
                                            caliper_list = gen_caliper_list(),
                                            sqrt_mahal = TRUE,
                                            keep_all_with_replacement = FALSE) {
    ## generate propensity score
    prop_score <- propensity_score(
        x_mat = x_mat,
        treat_vec = tolerance_list[["tolerance_vec"]],
        propensity_list = propensity_list
    )
    prop_dist_mat <- abs(outer(prop_score, prop_score, "-"))

    if (!is.null(caliper_list)) {
        prop_dist_mat <- prop_dist_mat + create_caliper(caliper_list,
            treat_vec = treat_vec
        )
    }

    nonbipartite_matches(
        dist_mat = prop_dist_mat,
        tolerance_list = tolerance_list,
        match_method = match_method,
        n_sinks = n_sinks,
        keep_all_with_replacement = keep_all_with_replacement
    )
}


#' Calculates propensity scores for a given matrix and treatment vector
#'
#' This function takes in an input matrix and a treatment vector,
#' along with a function that makes predictions (default xgboost method given)
#' and returns a predicted probability of treatment for each unit, either
#' using in-sample or out-of-sample fits.
#' @param x_mat Standard input matrix (already rank adjusted).
#' @param treat_vec Usual 0/1 treatment vector.
#' @param propensity_list See \code{gen_propensity_list}
#' @return Returns a vector equal in length to treat_vec of propensity.
#'   score.
#' @author Colman Humphrey
#'
#' @export
propensity_score <- function(x_mat,
                             treat_vec,
                             propensity_list = gen_propensity_list()) {
    propensity_function <- propensity_list[["propensity_function"]]

    if (!propensity_list[["oos_propensity"]]) {
        ## very simple...:
        train_treat_list <- list(
            x_train = x_mat,
            x_test = x_mat,
            y_train = treat_vec,
            y_test = treat_vec
        )
        return(propensity_function(train_treat_list))
    }

    ## ------------------------------------

    fold_res <- fold_indexing(nrow(x_mat), propensity_list[["n_folds"]])

    prediction_result <- lapply(fold_res, function(fold_inds) {
        fold_lgl <- (1L:nrow(x_mat)) %in% fold_inds
        train_treat_list <- list(
            x_train = x_mat[!fold_lgl, , drop = FALSE],
            x_test = x_mat[fold_inds, , drop = FALSE],
            y_train = treat_vec[!fold_lgl],
            y_test = treat_vec[fold_inds]
        )
        propensity_function(train_treat_list)
    })

    unlist(prediction_result)[order(unlist(fold_res))]
}


#' Function factory to predict treatment using xgboost
#'
#' This operates the exact same as \code{match_predict_xgb}, which
#' is itself quite a simple wrap of regular xgboost code.
#' Main difference here is the input data is even simpler.
#' The returned function accepts one parameter, \code{train_test_list},
#' a list with \code{x_train}, \code{x_test}, \code{y_train}, \code{y_test}
#' @inheritParams match_predict_xgb
#' @return returns a function that accepts \code{train_test_list}
#' and this returns a vector of predictions for the test data
#' @author Colman Humphrey
#'
#' @export
propensity_score_xgb <- function(nrounds = 50,
                                 nthread = 1,
                                 params = list(eta = 0.1, max.depth = 3),
                                 ...) {
    match_predict_xgb(nrounds = nrounds,
                      nthread = nthread,
                      params = params,
                      ...)
}


#' Function factory to predict treatment using \code{glm} (binomial)
#' or \code{lm}
#'
#' Does simple wrap around \code{glm} / \code{lm}
#' The returned function accepts one parameter, \code{train_test_list},
#' a list with \code{x_train}, \code{x_test}, \code{y_train}, \code{y_test}
#' @inheritParams propensity_score_linear
#' @return returns a function that accepts \code{train_test_list}
#' and this returns a vector of predictions for the test data
#' @author Colman Humphrey
#'
#' @export
propensity_score_linear <- function(use_linear_lm = FALSE) {
    function(train_test_list) {
        train_frame <- as.data.frame(train_test_list[["x_train"]])
        train_frame[["y"]] <- train_test_list[["y_train"]]

        test_frame <- as.data.frame(train_test_list[["x_test"]])

        if (use_linear_lm) {
            train_res <- lm(y ~ ., data = train_frame)
            lin_pred <- predict(train_res, newdata = test_frame,
                                type = "response")
            return(pmax(pmin(lin_pred, 1), 0))
        }

        train_res <- glm(y ~ ., data = train_frame, family = "binomial")
        predict(train_res, newdata = test_frame, type = "response")
    }
}
