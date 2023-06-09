% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{truncnorm_lg}
\alias{truncnorm_lg}
\title{Sample from a truncated normal distribution. Samples are drawn
componentwise, so each component of the vector is allowed its own
mean, standard deviation, and upper and lower limits. The components
are assumed to be independent.}
\usage{
truncnorm_lg(y_lower, y_upper, mu, sigma, u_rand)
}
\arguments{
\item{y_lower}{\code{n x p} matrix of lower endpoints}

\item{y_upper}{\code{n x p} matrix of upper endpoints}

\item{mu}{\code{n x p} matrix of conditional expectations}

\item{sigma}{\code{p x 1} vector of conditional standard deviations}

\item{u_rand}{\code{n x p} matrix of uniform random variables}
}
\value{
z_star \code{n x p} draw from the truncated normal distribution
}
\description{
Sample from a truncated normal distribution. Samples are drawn
componentwise, so each component of the vector is allowed its own
mean, standard deviation, and upper and lower limits. The components
are assumed to be independent.
}
\note{
This function uses \code{Rcpp} for computational efficiency.
}
