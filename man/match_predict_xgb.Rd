% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_predict.R
\name{match_predict_xgb}
\alias{match_predict_xgb}
\title{Function factory to predict treatment / control pairs using
xgboost}
\usage{
match_predict_xgb(
  nrounds = 50,
  nthread = 1,
  params = list(eta = 0.1, max.depth = 4),
  ...
)
}
\arguments{
\item{nrounds}{training rounds for the xgb algorithm}

\item{nthread}{threads to use for fitting, default 1...}

\item{params}{list of params to pass to xgboost,
most likely something like \code{eta} and \code{max.depth}}
}
\value{
returns a function that takes in a \code{train_test_list}
from \code{predict_prepare}; this function returns a
vector of predictions for the test data
}
\description{
The returned function takes in training and test data
(output from \code{predict_prepare}), trains an
xgboost model on the training, predicts on the test, and
returns the test vector
}
\author{
Colman Humphrey
}
