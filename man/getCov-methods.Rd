\name{getCov}
\alias{getCov}
\title{Get covariance matrix}
\usage{
getCov(n, model, tau)
}
\arguments{
  \item{n}{integer giving the size of the random feild}

  \item{model}{chatacter giving the name of the GMRF}

  \item{tau}{numeric giving the multiplier of the structure
  matrix for the model}
}
\value{
a covariance matrix
}
\description{
Returns the covariance matrix of the specified Gaussian
markov random field model
}

