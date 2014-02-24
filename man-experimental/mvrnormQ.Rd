\name{mvrnormQ}
\alias{mvrnormQ}
\title{samples from a multivariate normal given a mean and precision}
\usage{
  mvrnormQ(n, mu, Q, tol = 1e-06, return.density = FALSE)
}
\arguments{
  \item{n}{the number of samples to take}

  \item{mu}{the mean vector}

  \item{Q}{the precision matrix}

  \item{tol}{the tolerance for what is positive definate}

  \item{return.density}{whether to return the log density
  of the simulation.  Only makes sence when n is 1.}
}
\value{
  a pointer to the environment in which summaries of the
  data reside
}
\description{
  samples from a multivariate normal given a mean and
  precision
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

