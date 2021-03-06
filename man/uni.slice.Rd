% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uni.slice.R
\name{uni.slice}
\alias{uni.slice}
\title{Univariate slice sampling (Radford M. Neal (2008))}
\usage{
uni.slice(
  x0,
  g,
  w = 1,
  m = Inf,
  lower = -Inf,
  upper = +Inf,
  gx0 = NULL,
  uni.slice.calls = 0,
  uni.slice.evals = 0
)
}
\arguments{
\item{x0}{Initial point}

\item{g}{Function returning the log of the probability density (plus constant)}

\item{w}{Size of the steps for creating interval, defaults to 1}

\item{m}{Limit on steps, defaults to \eqn{\infty}}

\item{lower}{Lower bound on support of the distribution, defaults to \eqn{-\infty}}

\item{upper}{Upper bound on support of the distribution, defaults to \eqn{+\infty}}

\item{gx0}{Value of \eqn{g(x0)}, if known, defaults to \code{not known}}
}
\description{
Performs a slice sampling update from an initial point to a new point that leaves invariant the distribution with the specified log density function.
}
\examples{
uni.slice()
}
\keyword{uni.slice()}
