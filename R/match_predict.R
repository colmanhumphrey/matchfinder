#' Function factory to predict treatment / control pairs using
#' xgboost
#'
#' The returned function takes in training and test data
#' (output from \code{predict_prepare}), trains an
#' xgboost model on the training, predicts on the test, and
#' returns the test vector
#' @param nrounds training rounds for the xgb algorithm
#' @param nthread threads to use for fitting, default 1...
#' @param params list of params to pass to xgboost,
#'   most likely something like \code{eta} and \code{max.depth}
#' @return returns a function that takes in a \code{train_test_list}
#'   from \code{predict_prepare}; this function returns a
#'   vector of predictions for the test data
#' @author Colman Humphrey
#'
#' @export
match_predict_xgb <- function(nrounds = 50,
                              nthread = 1,
                              params = list(eta = 0.1, max.depth = 4),
                              ...) {
    if (!requireNamespace("xgboost")) {
        stop("xgboost not installed", call. = FALSE)
    }

    function(train_test_list) {
        xgb_train <- xgboost::xgb.DMatrix(
                                  train_test_list[["x_train"]],
                                  label = train_test_list[["y_train"]]
                              )

        xgb_test <- xgboost::xgb.DMatrix(
                                 train_test_list[["x_test"]],
                                 label = train_test_list[["y_test"]]
                             )

        train_res <- xgboost::xgb.train(
                                  params = params,
                                  data = xgb_train,
                                  nrounds = nrounds,
                                  nthread = nthread,
                                  objective = "binary:logistic",
                                  ...
                              )

        predict(train_res, newdata = xgb_test)
    }
}


#' Function factory to predict treatment / control pairs using
#' \code{glm} (binomial) or \code{lm}
#'
#' The returned function takes in training and test data
#' (output from \code{predict_prepare}), trains a
#' "linear" model on the training, predicts on the test, and
#' returns the test vector. NOTE: the glm model will
#' fail if the data is too "tricky", so be sure to check.
#' Also the pure \code{lm} model will be fast - so if either
#' of these cases is what you need, set \code{use_linear_lm}
#' to TRUE. In this case we'll bound the probabilities
#' returned at 0 and 1.
#' @param use_linear_lm logical, default FALSE; if you want
#'   to use regular \code{lm} instead of \code{glm}.
#' @return returns a function that takes in a \code{train_test_list}
#'   from \code{predict_prepare}; this function returns a
#'   vector of predictions for the test data
#' @author Colman Humphrey
#'
#' @export
match_predict_linear <- function(use_linear_lm = FALSE) {
    function(train_test_list) {
        train_frame <- as.data.frame(train_test_list[["x_train"]])
        train_frame[["y"]] <- train_test_list[["y_train"]]

        test_frame <- as.data.frame(train_test_list[["x_test"]])

        if (use_linear_lm) {
            ## we know the intercept should be 0.5 by symmetry
            ## but typically we don't suffer at all from not setting it
            train_res <- lm(y ~ ., data = train_frame)
            lin_pred <- predict(train_res, newdata = test_frame,
                                type = "response")
            return(pmax(pmin(lin_pred, 1), 0))
        }

        ## might fail if using cross methods
        ## because of the symmetry, we know the intercept
        ## is supposed to be zero so we can help the
        ## model out here (and it can indeed help!)
        train_res <- glm(y ~ . + 0, data = train_frame, family = "binomial")
        predict(train_res, newdata = test_frame, type = "response")
    }
}
