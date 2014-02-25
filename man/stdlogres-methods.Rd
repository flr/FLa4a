\docType{methods}
\name{stdlogres}
\alias{stdlogres}
\alias{stdlogres-methods}
\alias{stdlogres,FLQuant,FLQuant-method}
\title{Standardized log residuals}
\arguments{
  \item{obs}{a \code{FLQuant} object with the observations}

  \item{fit}{a \code{FLQuant} object with the fitted value}
}
\value{
  a \code{FLQuant} with stardardized log residuals
}
\description{
  Method o compute the standardized residuals in the log
  scale for the a4a stock assessment framework. Meaning
  indices and catch-at-age residuals.
}
\examples{
library(FLa4a)
data(ple4)
data(ple4.index)
fit. <- a4a(stock=ple4, indices=FLIndices(ple4.index))
stdlogres(catch.n(ple4), catch.n(fit.))
}

