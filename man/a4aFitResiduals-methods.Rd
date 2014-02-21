\docType{methods}
\name{a4aFit}
\alias{a4aFit}
\alias{a4aFit,a4aFit-method}
\title{Computes log residuals of catches and indices}
\arguments{
  \item{a4aFit}{object}

  \item{...}{Additional argument list that might not ever
  be used.}
}
\value{
  an FLQuant
}
\description{
  Computes log residuals of catches and indices
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
\seealso{
  \code{\link{print}} and \code{\link{cat}}
}

