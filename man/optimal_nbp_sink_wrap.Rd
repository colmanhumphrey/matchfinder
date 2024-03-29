% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrap_matches.R
\name{optimal_nbp_sink_wrap}
\alias{optimal_nbp_sink_wrap}
\title{Given a vector of sink values, produces optimal NBP matches for each}
\usage{
optimal_nbp_sink_wrap(
  dist_mat,
  tolerance_vec,
  n_sinks = NULL,
  weight_vec = NULL
)
}
\arguments{
\item{dist_mat}{Matrix of pairwise distances.}

\item{n_sinks}{Vector of sinks per match. Note that this is NOT the same
for each \code{match_method}: for both \code{"with_replacement"}
and \code{"greedy"}, it subtracts one pair. But for optimal matching,
it removes one full unit. \code{"greedy"} is the weird one here, but it
wouldn't truly make sense given the code to try and replicate the optimal
implementation. Be safe and don't use greedy...}

\item{weight_vec}{Default \code{NULL}: optionally supply the weight vector
used to generate \code{dist_mat} and it'll be returned in the
\code{match_list} generated from this function}
}
\description{
Given a vector of sink values, produces optimal NBP matches for each
}
\author{
Colman Humphrey
}
\keyword{internal}
