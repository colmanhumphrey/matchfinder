#' Function to predict treatment / control pairs using
#' xgboost
#'
#' This function takes in training and test data, trains an
#' xgboost model on the training, predicts on the test, and
#' returns the test vector
#' @param train_test_list output from \code{predict_prepare}
#' @param nrounds training rounds for the xgb algorithm
#' @param nthread threads to use for fitting, default 1...
#' @param params list of params to pass to xgboost,
#'   most likely something like \code{eta} and \code{max.depth}
#' @return returns a vector of predictions for the test data
#' @author Colman Humphrey
#'
#' @export
match_predict_xgb <- function(train_test_list,
                              nrounds = 50,
                              nthread = 1,
                              params = list(eta = 0.1, max.depth = 4),
                              ...) {
    if (!requireNamespace("xgboost")) {
        stop("xgboost not installed", call. = FALSE)
    }

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


#' Function to predict treatment / control pairs using
#' \code{glm} (binomial) or \code{lm}
#'
#' This function takes in training and test data, trains a
#' "linear" model on the training, predicts on the test, and
#' returns the test vector. NOTE: the glm model will
#' fail if the data is too "tricky", so be sure to check.
#' Also the pure \code{lm} model will be fast
#' @param train_test_list output from \code{predict_prepare}
#' @param use_lm_approx logical, default FALSE; if
#' @return returns a vector of predictions for the test data
#' @author Colman Humphrey
match_predict_linear <- function(train_test_list,
                                 use_lm_approx = FALSE) {
    train_frame <- as.data.frame(train_test_list[["x_train"]])
    train_frame[["y"]] <- train_test_list[["y_train"]]

    test_frame <- as.data.frame(train_test_list[["x_test"]])

    if (use_lm_approx) {
        train_res <- lm(y ~ ., data = train_frame)
        lin_pred <- predict(train_res, newdata = test_frame, type = "response")
        return(pmax(pmin(lin_pred, 1), 0))
    }

    ## might fail if using cross methods
    train_res <- glm(y ~ . + 0, data = train_frame, family = "binomial")
    predict(train_res, newdata = test_frame, type = "response")
}
