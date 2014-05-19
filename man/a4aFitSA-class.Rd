\docType{class}
\name{a4aFitSA-class}
\alias{a4aFit,a4aFitSA-method}
\alias{a4aFitSA}
\alias{a4aFitSA,a4aFit-method}
\alias{a4aFitSA,missing-method}
\alias{a4aFitSA-class}
\alias{a4aFitSA-methods}
\alias{a4aFitSAs}
\alias{a4aFitSAs,a4aFitSA-method}
\alias{a4aFitSAs,list-method}
\alias{a4aFitSAs,missing-method}
\alias{a4aFitSAs-class}
\alias{m,a4aFitSA-method}
\alias{pars,a4aFitSA-method}
\alias{show,a4aFitSA-method}
\alias{submodels,a4aFitSA-method}
\title{S4 class \code{a4aFitSA}}
\description{
  The \code{a4aFitSA} class extends \code{a4aFit} to store
  information about the parameters of the model.
}
\section{Slots}{
  \describe{ \item{SCAPars}{An object of class
  \code{SCAPars} with information about model parameters} }
}

\section{Accessors}{
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
obj <- sca(stock=ple4, indices=FLIndices(ple4.index), fit="assessment")
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
is(pars(obj))
}

