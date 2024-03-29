% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonbipartite_matches.R
\name{add_nbp_sinks}
\alias{add_nbp_sinks}
\title{This function gives a distance matrix that nbpmatching likes}
\usage{
add_nbp_sinks(dist_mat, n_sinks = 0L)
}
\arguments{
\item{dist_mat}{A symmetric distance matrix, e.g. result of
\code{weighted_mahal}}

\item{n_sinks}{How many potential matches to throw away?}
}
\description{
Note: if your \code{dist_mat} is even (e.g. number of rows/cols),
you should supply an even number of sinks, and odd if odd -
else you'll get a ghost value added and a warning. This is true
even for adding zero sinks to an odd.
}
\keyword{internal}
