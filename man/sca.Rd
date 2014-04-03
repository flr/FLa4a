\docType{methods}
\name{sca}
\alias{sca}
\alias{sca,FLStock,FLIndex-method}
\alias{sca,FLStock,FLIndices-method}
\alias{sca-methods}
\title{statistical catch-at-age method}
\arguments{
  \item{stock}{a \code{FLStock} object}

  \item{indices}{a \code{FLIndices} object}

  \item{fmodel}{a formula object depicting the model for
  log fishing mortality at age}

  \item{qmodel}{a list of formula objects depicting the
  models for log survey catchability at age}

  \item{srmodel}{a formula object depicting the model for
  log recruitment}

  \item{fit}{Character with type of fit, 'MP' or
  'assessment', the first doesn't require the hessian to be
  computed, while the former does.}
}
\value{
  a \code{a4aFit} or \code{a4aFitSA} object with the fit
  results.
}
\description{
  User interface to the statistical catch-at-age method of
  the a4a stock assessment framework.
}
\examples{
data(ple4)
data(ple4.index)
fit1 <- sca(ple4, FLIndices(ple4.index))
plot(ple4 + fit1)
fit2 <- sca(ple4, FLIndices(ple4.index), fmodel=~factor(age) + factor(year))
plot(ple4 + fit2)
fit3 <- sca(ple4, FLIndices(ple4.index), fmodel=~factor(age) + factor(year), qmodel=list(~s(age, k=4)))
plot(ple4 + fit3)
fit4 <- sca(ple4, FLIndices(ple4.index), fmodel=~factor(age) + factor(year), qmodel=list(~s(age, k=4)), srmodel=~s(year, k=45))
plot(ple4 + fit4)
AIC(fit1, fit2, fit3, fit4)
BIC(fit1, fit2, fit3, fit4)
}

