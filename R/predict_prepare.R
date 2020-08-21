#' Prepares design matrices and other inputs needed for predictive methods
#'
#' It is not immediately clear how to take a match (set of matched pairs)
#' and "predict" it - essentially you likely want some prediction
#' function \eqn{f(x, y) -> {0, 1}} (or \eqn{P(1)} etc) where \eqn{x}
#' and \eqn{y} are rows, one of which is treated and the other control.
#' This function should also have a property resembling:
#' \deqn{f(x, y) = 1 - f(y, x)} since there's an inherent symmetry
#' in the pairs: if you knew which of the two came first, there'd
#' be nothing to do. In fact, the definition of \eqn{f} would likely be
#' to predict if (WLOG) \eqn{x} is the treated unit or not.
#'
#' This function gives various potential design matrix shapes
#' and corresponding output vectors that all can be used to solve the
#' above prediction problem. In some scenarios (e.g. if we were content
#' to use linear models), differences alone might be appropriate. For full
#' generality, we may wish to use both full vectors.
#' @param x_mat original data, i.e. row for each unit
#' @param index_list list of indices to use in the matching:
#'   \describe{
#'     \item{\code{treat_train}}{Index of treated units for the training matrix}
#'     \item{\code{control_train}}{Index of control units for the training matrix}
#'     \item{\code{treat_test}}{Index of treated units for the test matrix}
#'     \item{\code{control_test}}{Index of control units for the test matrix}
#'   }
#' @param design The design matrix form do you want to use.
#'   We'll assume that \code{x_mat} is \eqn{n \times k}.
#'   Let \eqn{X_c} be the rows of \code{x_mat} corresponding
#'   to the control units, and \eqn{X_t} the treated units that
#'   pair with \eqn{X_c} (i.e. so the \eqn{i^{\text{th}}} matched
#'   pair will correspond with the \eqn{i^{\text{th}}} row of each),
#'   and thus each matrix will have dimension \eqn{m \times k},
#'   where \eqn{2m \leq n}.
#'   \describe{
#'     \item{\code{"cross_all"}}{
#'   This will form the \eqn{2m \times 2k} matrix, by having each matched pair
#'   give two rows in the design matrix: each of the two
#'   rows will contain the full treated row of the given pair along with the full
#'   control row, but in one row the pair will be (treated, control), and in the other
#'   (control, treated). For these two rows, the output vector will contain a 1 and a 0
#'   respectively. We construct the following inputs:
#'   \deqn{X =
#'     \begin{bmatrix}
#'     X_t& X_c \\
#'     X_c& X_t \\
#'     \end{bmatrix}, Y = [1 1 \ldots 1 0 0 \ldots 0]'
#'   }
#'   For the test data, we'll only (WLOG more or less) give the (treat, control)
#'   version, for consistency in evaluation.
#'   }
#'      \item{\code{"cross_random"}}{
#'   This is similar to the above, except instead of forming all rows in both
#'   ways, we choose (for each pair) between (treated, control) and
#'   (control, treated) randomly. So
#'   \deqn{Y_i =
#'   \begin{cases}
#'   1 & \text{with probability} 0.5 \\
#'   0 & \text{with probability} 0.5
#'   \end{cases}}
#'   And (\eqn{X_i} here meaning the \eqn{i^{\text{th}}} row of the design matrix):
#'   \deqn{X_i =
#'   \begin{cases}
#'   ((X_t)_i, (X_c)_i) & Y_i = 1 \\
#'   ((X_c)_i, (X_t)_i) & Y_i = 0
#'   \end{cases}}
#'   Giving an \eqn{m \times 2k} design matrix.
#'   }
#'     \item{\code{"differences_random"}}{For
#'   each pair, our design row will be
#'   the row in \code{x_mat} corresponding to
#'   one unit minus the row corresponding to the
#'   other. Similar to the above, we'll choose
#'   the outcome randomly and decide which unit will
#'   be subtracted from the other based on that, and we'll end up with
#'   an \eqn{m \times k} matrix.
#'   That is:
#'   \deqn{Y_i =
#'   \begin{cases}
#'   1 & \text{with probability} 0.5 \\
#'   0 & \text{with probability} 0.5
#'   \end{cases}}
#'   And:
#'   \deqn{X_i =
#'   \begin{cases}
#'   (X_t)_i - (X_c)_i & Y_i = 1 \\
#'   (X_c)_i - (X_t)_i & Y_i = 0
#'   \end{cases}}
#'   Note that if your \code{x_mat} contains an intercept, it'll
#'   now be a column of zeroes.
#'   }
#'   \item{\code{"differences_plain"}}{
#'   In some cases, you may wish to have the differences and add
#'   customization. So here we essentially do the same as above
#'   but with \eqn{Y_i = 1} for all \eqn{i}, therefore
#'   all differences are returned as \eqn{(X_t)_i - (X_c)_i}
#'   (and will also be \eqn{m \times k}).
#'   }
#'   }
#' @return list:
#'   \describe{
#'     \item{\code{x_train}}{Design matrix for training}
#'     \item{\code{x_test}}{Design matrix for testing}
#'     \item{\code{y_train}}{Training outcome vector}
#'     \item{\code{y_test}}{Test outcome vector}
#' }
#' @author Colman Humphrey
#'
#' @export
predict_prepare <- function(x_mat,
                            index_list,
                            design = c(
                                "cross_all",
                                "cross_random",
                                "differences_random",
                                "differences_plain"
                            )) {
    design <- match.arg(design)

    stopifnot(length(index_list[["treat_train"]]) ==
        length(index_list[["control_train"]]) &&
        length(index_list[["treat_test"]]) ==
            length(index_list[["control_test"]]))

    ## ------------------------------------

    treat_train_ind <- index_list[["treat_train"]]
    control_train_ind <- index_list[["control_train"]]

    treat_test_ind <- index_list[["treat_test"]]
    control_test_ind <- index_list[["control_test"]]

    ## ------------------------------------

    x_treat_train <- x_mat[treat_train_ind, , drop = FALSE]
    x_control_train <- x_mat[control_train_ind, , drop = FALSE]

    x_treat_test <- x_mat[treat_test_ind, , drop = FALSE]
    x_control_test <- x_mat[control_test_ind, , drop = FALSE]

    train_index <- seq_len(length(treat_train_ind))
    test_index <- seq_len(length(treat_test_ind))

    ## y_train_full <-

    ## ------------------------------------

    if (design %in% c("cross_random", "differences_random")) {
        ## basically exactly half in each

        treat_left_train <- train_index %in%
            sample(length(treat_train_ind), size = length(treat_train_ind) / 2L)

        treat_left_test <- test_index %in%
            sample(length(treat_test_ind), size = length(treat_test_ind) / 2L)
    } else {
        ## cross_all has all, differences_random just has left
        treat_left_train <- rep(TRUE, times = length(treat_train_ind))
        treat_left_test <- rep(TRUE, times = length(treat_test_ind))
    }

    if (design == "cross_all") {
        treat_right_train <- rep(TRUE, times = length(treat_train_ind))
        treat_right_test <- rep(FALSE, times = length(treat_test_ind))

        train_order <- seq_len(length(treat_train_ind) * 2L)
        test_order <- seq_len(length(treat_test_ind))

        y_train <- rep(c(1L, 0L),
            each = length(treat_train_ind)
        )
        y_test <- rep(1L, times = length(treat_test_ind))
    } else {
        treat_right_train <- !treat_left_train
        treat_right_test <- !treat_left_test

        train_order <- order(c(
            train_index[treat_left_train],
            train_index[treat_right_train]
        ))
        test_order <- order(c(
            test_index[treat_left_test],
            test_index[treat_right_test]
        ))

        y_train <- treat_left_train * 1L
        y_test <- treat_left_test * 1L
    }

    if (design %in% c("cross_all", "cross_random")) {
        ## these will "cross" the rows, and use all
        gen_design_mat <- function(x, y) {
            cbind(x, y)
        }
    } else {
        ## differences regime
        gen_design_mat <- function(x, y) {
            x - y
        }
    }

    ## ------------------------------------

    x_tc_train <- gen_design_mat(
        x_treat_train[treat_left_train, , drop = FALSE],
        x_control_train[treat_left_train, , drop = FALSE]
    )
    x_ct_train <- gen_design_mat(
        x_control_train[treat_right_train, , drop = FALSE],
        x_treat_train[treat_right_train, , drop = FALSE]
    )

    x_tc_test <- gen_design_mat(
        x_treat_test[treat_left_test, , drop = FALSE],
        x_control_test[treat_left_test, , drop = FALSE]
    )
    x_ct_test <- gen_design_mat(
        x_control_test[treat_right_test, , drop = FALSE],
        x_treat_test[treat_right_test, , drop = FALSE]
    )

    ## combining to form full input
    list(
        x_train = rbind(x_tc_train, x_ct_train)[train_order, , drop = FALSE],
        x_test = rbind(x_tc_test, x_ct_test)[test_order, , drop = FALSE],
        y_train = y_train,
        y_test = y_test
    )
}


#' Constructs an \code{index_list} from a \code{match_list} and a logical
#' index for training data
#'
#' @param match_list see \code{bipartite_matches} etc
#' @param train_index logical index, same length as
#'   \code{match_list[["treat_index"]]} (and
#'   \code{match_list[["control_index"]]})
#' @return returns an \code{index_list} object, see
#'   \code{predict_prepare}
#' @author Colman Humphrey
#'
#' @export
index_list_from_match <- function(match_list,
                                  train_index) {
    list(
        treat_train = match_list[["treat_index"]][train_index],
        control_train = match_list[["control_index"]][train_index],
        treat_test = match_list[["treat_index"]][!train_index],
        control_test = match_list[["control_index"]][!train_index]
    )
}


#' Creates an \code{index_list} from a \code{match_list}, splitting
#' according to \code{train_fraction}
#'
#' @param match_list typical \code{match_list} entry
#' @param train_fraction fraction (between 0 and 1) to
#'   use for training data (and the rest for test)
#' @author Colman Humphrey
#'
#' @export
generate_train_test_split <- function(match_list,
                                      train_fraction = 0.7) {
    stopifnot(train_fraction >= 0 && train_fraction <= 1)

    length_index <- length(match_list[["treat_index"]])
    all_index <- 1L:length_index

    train_index <- all_index %in%
        sample(all_index, size = floor(train_fraction * length_index))
    index_list_from_match(
        match_list,
        train_index
    )
}


#' Constructs a k-fold list of \code{index_list} objects for a given
#' \code{match_list}
#'
#' @param match_list typical \code{match_list} entry
#' @param num_folds how many folds you want, default 5
#' @author Colman Humphrey
#'
#' @export
generate_k_fold_index <- function(match_list,
                                  num_folds = 5L) {
    length_index <- length(match_list[["treat_index"]])
    all_index <- 1L:length_index

    fold_res <- fold_indexing(length_index, num_folds)

    lapply(fold_res, function(inds) {
        index_list_from_match(
            match_list,
            !(all_index %in% inds)
        )
    })
}
