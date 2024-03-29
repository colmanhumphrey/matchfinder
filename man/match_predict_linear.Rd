% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_predict.R
\name{match_predict_linear}
\alias{match_predict_linear}
\title{Function factory to predict treatment / control pairs using
\code{glm} (binomial) or \code{lm}}
\usage{
match_predict_linear(use_linear_lm = FALSE)
}
\arguments{
\item{use_linear_lm}{logical, default FALSE; if you want
to use regular \code{lm} instead of \code{glm}.}
}
\value{
returns a function that takes in a \code{train_test_list}
from \code{predict_prepare}; this function returns a
vector of predictions for the test data
}
\description{
The returned function takes in training and test data
(output from \code{predict_prepare}), trains a
"linear" model on the training, predicts on the test, and
returns the test vector. NOTE: the glm model will
fail if the data is too "tricky", so be sure to check.
Also the pure \code{lm} model will be fast - so if either
of these cases is what you need, set \code{use_linear_lm}
to TRUE. In this case we'll bound the probabilities
returned at 0 and 1.
}
\author{
Colman Humphrey
}
