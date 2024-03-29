% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{calc_brier}
\alias{calc_brier}
\title{Simple brier score calc}
\usage{
calc_brier(predict, outcome, avg = TRUE)
}
\arguments{
\item{predict}{\eqn{P(Y = 1)} for each value}

\item{outcome}{result (in \eqn{\{0, 1\}})}

\item{avg}{logica, default TRUE, do you want the mean?}
}
\value{
length one double: the total brier sum, or the avg (default)
}
\description{
Simple brier score calc
}
\author{
Colman Humphrey
}
\keyword{internal}
