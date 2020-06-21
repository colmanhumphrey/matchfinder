#' Computes a simple mean difference in an outcome vector
#' between treatment and control in a paired match
#'
#' Computes the average difference between the treated units
#' and the control units for a match, given as a match list.
#' Optionally confirms that all treated units are indeed treated.
#' @param match_list Typical \code{match_list} object from
#'   \code{{non}bipartite_matches}.
#' @param y_vector The outcome vector.
#' @param treat_vec Default NULL, provide if you want it checked.
#' @return Returns a single number, the mean difference.
#' @author Colman Humphrey
#'
#' @export
match_estimate <- function(match_list,
                           y_vector,
                           treat_vec = NULL) {
    if (!is.null(treat_vec)) {
        stopifnot(length(treat_vec) == length(y_vector))
        stopifnot(all(treat_vec[match_list[["treat_index"]]] == 1L))
        stopifnot(all(treat_vec[match_list[["control_index"]]] == 0L))
    }

    mean(y_vector[match_list[["treat_index"]]] -
        y_vector[match_list[["control_index"]]])
}


#' Takes in training/test data and a prediction
#' function to use, generates brier score
#'
#' @param train_test_list output from \code{predict_prepare}
#' @param match_predict_function function to predict treated unit
#' @param avg logical, default TRUE: should we average (or sum) the briers?
#' @return result from \code{calc_brier}, length-one double
#' @author Colman Humphrey
#'
#' @export
brier_score <- function(train_test_list,
                        match_predict_function = match_predict_xgb,
                        avg = TRUE) {
    calc_brier(match_predict_function(train_test_list),
        train_test_list[["y_test"]],
        avg = avg
    )
}


#' For a match, calculates brier score on a test split of the data
#'
#' @inheritParams permutation_brier
#' @param train_fraction split of data to use for training
#' @return brier score
#' @author Colman Humphrey
#'
#' @export
brier_score_split <- function(x_mat,
                              match_list,
                              design = "cross_all",
                              train_fraction = 0.7,
                              match_predict_function = match_predict_xgb) {
    brier_score(predict_prepare(
        x_mat,
        generate_train_test_split(match_list, train_fraction),
        design
    ),
    match_predict_function,
    avg = TRUE
    )
}


#' For a match, calculates brier score using cross validation
#'
#' @inheritParams permutation_brier
#' @param num_folds how many CV folds to use
#' @return brier score (averaged over all units)
#' @author Colman Humphrey
#'
#' @export
brier_score_cv <- function(x_mat,
                           match_list,
                           design = "cross_all",
                           num_folds = 5,
                           match_predict_function = match_predict_xgb) {
    k_fold_lists <- generate_k_fold_index(
        match_list,
        num_folds
    )

    pred_list <- lapply(k_fold_lists, function(k_fold) {
        train_test_list <- predict_prepare(x_mat,
            k_fold,
            design = design
        )
        c(
            brier_score(train_test_list,
                match_predict_function,
                avg = FALSE
            ),
            length(train_test_list[["y_test"]])
        )
    })

    sums <- Reduce("+", pred_list)
    sums[1] / sums[2]
}


#' For a given match, computes the brier score distribution
#' if the pairing were truly random
#'
#' @param x_mat typical input matrix
#' @param match_list match result
#' @param design see \code{predict_prepare}
#' @param use_cv logical, default TRUE: use CV to get briers? Else
#'   split.
#' @param num_permutations how many permutations to do
#' @param match_predict_function function to predict treated units
#' @param num_folds if using CV, how many folds?
#' @param train_fraction if using split, fraction to train?
#' @return vector of brier scores for random pairings
#' @author Colman Humphrey
#'
#' @export
permutation_brier <- function(x_mat,
                              match_list,
                              design = "cross_all",
                              use_cv = TRUE,
                              num_permutations = 100L,
                              match_predict_function = match_predict_xgb,
                              num_folds = 5,
                              train_fraction = 0.7) {
    if (use_cv) {
        if (!missing(train_fraction)) {
            stop(
                "only set `train_fraction` if not using cross-validation ",
                "(set `use_cv = FALSE`)"
            )
        }

        brier_function <- (function(x_mat,
                                    design,
                                    num_folds,
                                    match_predict_function) {
            function(match_list) {
                brier_score_cv(
                    x_mat,
                    match_list,
                    design,
                    num_folds,
                    match_predict_function
                )
            }
        })(x_mat, design, num_folds, match_predict_function)
    } else {
        if (!missing(num_folds)) {
            stop(
                "only set `num_folds` if using cross-validation ",
                "(set `use_cv = TRUE`)"
            )
        }

        brier_function <- (function(x_mat,
                                    design,
                                    num_folds,
                                    match_predict_function) {
            function(match_list) {
                brier_score_split(
                    x_mat,
                    match_list,
                    design,
                    train_fraction,
                    match_predict_function
                )
            }
        })(x_mat, design, train_fraction, match_predict_function)
    }

    ## ------------------------------------

    full_index <- 1L:length(match_list[["treat_index"]])

    unlist(lapply(1L:num_permutations, function(i) {
        swaps <- full_index %in% sample(full_index,
            size = floor(length(full_index) / 2)
        )

        swap_match <- swap_pairs(
            match_list,
            swaps
        )

        brier_function(swap_match)
    }))
}
