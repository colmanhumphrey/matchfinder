% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonbipartite_matches.R
\name{gen_tolerance_list}
\alias{gen_tolerance_list}
\title{Wrapper to unify tolerance input, along with minor checks}
\usage{
gen_tolerance_list(
  tolerance_vec = NULL,
  tolerance_min = NULL,
  tolerance_max = NULL
)
}
\arguments{
\item{tolerance_vec}{Default NULL; numeric "continuous treatment"
vector that we use to form
non-bipartite matches: units i and j can be matched if
\eqn{\mid \texttt{tolerance_vec}[i] - \texttt{tolerance_vec}[j]\mid >
  \texttt{tolerance_min}}.}

\item{tolerance_min}{See above for what this does - blocks
matches that are too close on \code{tolerance_vec}. Something like minimum
relevant difference to be a "treatment" effect. Default \code{NULL} gives
zero, i.e. only blocks equality of \code{tolerance_vec}.}

\item{tolerance_max}{Optionally we may want to also ensure
our "treatment" values aren't too far apart. E.g. we may think our
assumptions are reasonable valid for small differences
in the tolerance vector, but not for large. Or another way:
we're asking for say marginal effects: how bad is one extra beer a
day?}
}
\value{
Either \code{NULL}, or a list with the same names as the input,
with validated values.
}
\description{
Note that we will always block equal equality on the tolerance
vec, hence why \code{tolerance_min} being \code{NULL} is equivalent
to zero. For more sophisticated behaviour, you may want to
directly control the caliper list or weighted distance matrix
in the functions that generate matches.
}
\author{
Colman Humphrey
}
