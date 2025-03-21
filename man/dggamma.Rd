% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggamma.R
\name{dggamma}
\alias{dggamma}
\title{Generalized Gamma Density Function}
\usage{
dggamma(
  x,
  mu,
  sigma,
  nu,
  return_log = FALSE,
  negate_loglik = FALSE,
  control = expansionControl()
)
}
\arguments{
\item{x}{Numeric vector of values where the density is evaluated.}

\item{mu}{Positive numeric, location parameter.}

\item{sigma}{Positive numeric, scale parameter.}

\item{nu}{Numeric, shape parameter.}

\item{return_log}{Logical, whether to return the log-density. Default is FALSE.}

\item{negate_loglik}{Logical, whether to negate the log-likelihood. Default is FALSE.}

\item{control}{A list of control parameters, typically created using \code{expansionControl()}.}
}
\value{
A numeric vector of density values.
}
\description{
Computes the density of the generalized gamma distribution using a series expansion.
}
\details{
The control parameters are \code{epsilon} and \code{kappa}. \code{epsilon} is a length two numeric vector
indicating the intervals \code{c(-epsilon[2],epsilon[1])} and \code{c(epsilon[1],epsilon[2])} where the
transition is made from the standard generalized gamma log-lik to its approximation. \code{kappa}
determines the maximal power of \code{nu} to consider in the Taylor and Stirling expansions formulas.
}
\examples{
dggamma(x = c(1, 2, 3), mu = 2, sigma = 1, nu = 0.5)
dggamma(x = 1:5, mu = 3, sigma = 2, nu = -0.2, return_log = TRUE)
}
