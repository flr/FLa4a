\docType{methods}
\name{stdlogres}
\alias{get}
\alias{K}
\alias{Method}
\alias{stdlogres}
\alias{stdlogres,stdlogres-method}
\alias{to}
\alias{values}
\title{Standardized log residuals}
\arguments{
  \item{obs}{a \code{FLQuant} object with the observations}

  \item{fit}{a \code{FLQuant} object with the fitted value}
}
\value{
  a \code{FLQuant} with stardardized log residuals
}
\description{
  Standardized log residuals
}
\examples{
library(FLa4a)
data(ple4)
data(ple4.index)
fmodel <- ~factor(age) + factor(year)
qmodel <- list(~factor(age))
fit. <- a4a(stock=ple4, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(ple4.index), fit ="assessment", wkdir="test")
stdlogres(catch.n(ple4), catch.n(fit.))
}
\author{
  EJ \email{ernesto.jardim@jrc.ec.europa.eu}
}

