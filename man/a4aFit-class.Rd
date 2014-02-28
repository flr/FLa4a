\docType{class}
\name{a4aFit-class}
\alias{a4aFit}
\alias{a4aFit-class}
\alias{a4aFit-methods}
\alias{a4aFit,missing-method}
\alias{catch.n}
\alias{catch.n,a4aFit-method}
\alias{clock}
\alias{clock,a4aFit-method}
\alias{fitSumm}
\alias{harvest}
\alias{harvest,a4aFit-method}
\alias{index}
\alias{index,a4aFit-method}
\alias{logLik}
\alias{logLik,a4aFit-method}
\alias{show}
\alias{show,a4aFit-method}
\alias{stock.n}
\alias{stock.n,a4aFit-method}
\title{S4 class \code{a4aFit}}
\usage{
  va4aFit(object)
}
\description{
  The \code{a4aFit} class was built to store the a4a stock
  assessment fits.
}
\section{Slots}{
  \describe{ \item{call}{The function call}

  \item{clock}{Information on call duration}

  \item{fitSumm}{Fit summary}

  \item{stock.n}{Estimates of stock numbers-at-age}

  \item{harvest}{Estimates of fishing mortality at age}

  \item{catch.n}{Estimates of catch numbers-at-age}

  \item{index}{Estimates of indices at age}

  }
}

\section{Accesors}{
  All slots in the class have accessor and replacement
  methods defined that allow retrieving and substituting
  individual slots.

  The values passed for replacement need to be of the class
  of that slot. A numeric vector can also be used when
  replacing FLQuant slots, and the vector will be used to
  substitute the values in the slot, but not its other
  attributes.
}

\section{Constructor}{
  A construction method exists for this class that can take
  named arguments for any of its slots. All slots are then
  created to match the requirements of the class validity.
  If an unnamed \code{FLQuant} object is provided, this is
  used for sizing but not stored in any slot.
}
\examples{
data(ple4)
data(ple4.index)
obj <- sca(stock=ple4, indices=FLIndices(ple4.index))
obj
clock(obj)
fitSumm(obj)
flq <- stock.n(obj)
is(flq)
flq <- catch.n(obj)
is(flq)
flq <- harvest(obj)
is(flq)
flq <- index(obj)
is(flq)
logLik(obj)
}

