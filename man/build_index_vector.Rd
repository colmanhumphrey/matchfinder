% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{build_index_vector}
\alias{build_index_vector}
\title{Takes a vector, returns a named list of repeated elements
along with indices.}
\usage{
build_index_vector(vec)
}
\arguments{
\item{vec}{Vector of potentially repeating elements that we
want the index locations of}
}
\value{
List with named elements, names from the input \code{vec},
each element is a list of indices
}
\description{
We have a vector of potentially repeating elements.
We want all elements with repeats, along with
the index of all locations of that value.
So e.g. if \code{vec = c("a", "d", "f", "d", "f")},
we'll get a list with element \code{"d"}
being a vector with values \code{c(2, 4)}
and \code{"f"} having \code{c(3, 5)}.
}
\author{
Colman Humphrey
}
