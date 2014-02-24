\docType{class}
\name{a4aFitResiduals-class}
\alias{a4aFitResiduals-class}
\alias{residuals}
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
fmodel <- ~factor(age) + factor(year)
qmodel <- list(~factor(age))
fit. <- a4a(stock=ple4, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(ple4.index), fit ="assessment", wkdir="test")
flqs <- residuals(fit., ple4, FLIndices(idx=ple4.index))
}

