% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_eval.R
\name{match_estimate}
\alias{match_estimate}
\title{Computes a simple mean difference in an outcome vector
between treatment and control in a paired match}
\usage{
match_estimate(match_list, y_vector, treat_vec = NULL)
}
\arguments{
\item{match_list}{Typical \code{match_list} object from
\code{{non}bipartite_matches}.}

\item{y_vector}{The outcome vector.}

\item{treat_vec}{Default NULL, provide if you want it checked.}
}
\value{
Returns a single number, the mean difference.
}
\description{
Computes the average difference between the treated units
and the control units for a match, given as a match list.
Optionally confirms that all treated units are indeed treated.
}
\author{
Colman Humphrey
}