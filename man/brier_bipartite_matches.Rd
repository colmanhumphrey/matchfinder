% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_weighted_matches.R
\name{brier_bipartite_matches}
\alias{brier_bipartite_matches}
\title{Computes all matches, then gets the brier scores for each. Reorder by
number of sinks.}
\usage{
brier_bipartite_matches(
  x_mat,
  cov_x,
  weight_list,
  treat_vec,
  match_method = c("with_replacement", "optimal", "greedy"),
  n_sinks = 0L,
  caliper_list = gen_caliper_list(),
  propensity_list = match_propensity_list(NULL),
  sqrt_mahal = TRUE,
  tol_val = NULL,
  design = "cross_all",
  num_folds = 5,
  match_predict_function = match_predict_xgb(),
  silent = !interactive()
)
}
\arguments{
\item{x_mat}{input/design matrix (already rank-adjusted etc)}

\item{cov_x}{the (potentially rank-adjusted) covariance of \code{x_mat}.
This means it's possible that \code{cov(x_mat)} is not equal to
\code{cov_x}; see \code{covariance_with_ranks} for more details.}

\item{weight_list}{list of weight vectors. See \code{generate_random_weights} to
automatically generate a reasonable set of vectors.}

\item{treat_vec}{Logical (or 1/0) vector, indicating treatment (or control).}

\item{n_sinks}{how many potential matches to not bother with
NOTE: you can do this as a vector, but not for optimal matching.}

\item{caliper_list}{Optional, see \code{gen_caliper_list}. Provide
this to force matches that are close on some metric.}

\item{tol_val}{For optimal matches, you can set a tolerance
to be within optimality of, which can be zero for perfect optimality.
Default 1e-4 is reasonable in many cases.}

\item{design}{see \code{predict_prepare}}

\item{num_folds}{how many CV folds to use}

\item{match_predict_function}{function to predict treated units}
}
\value{
List of matches within sink values,
and brier scores for each.
}
\description{
Computes all matches, then gets the brier scores for each. Reorder by
number of sinks.
}
\author{
Colman Humphrey
}
