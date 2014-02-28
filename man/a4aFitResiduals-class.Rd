\docType{class}
\name{a4aFitResiduals-class}
\alias{a4aFitResiduals}
\alias{a4aFitResiduals-class}
\alias{a4aFitResiduals-methods}
\alias{residuals}
\alias{residuals,a4aFit-method}
\title{S4 class \code{a4aFitResiduals}}
\description{
  The \code{a4aFitResiduals} class extends \code{FLQuants}
  to store residuals of the a4a stock assessment fit. By
  default these should be log residuals of catches and
  indices.
}
\examples{
library(FLa4a)
data(ple4)
data(ple4.index)
fit. <- sca(stock=ple4, indices=FLIndices(ple4.index))
flqs <- residuals(fit., ple4, FLIndices(idx=ple4.index))
}

