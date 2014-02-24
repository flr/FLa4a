\docType{class}
\name{a4aFitMCMC-class}
\alias{a4aFit}
\alias{a4aFitMCMC}
\alias{a4aFitMCMC-class}
\title{S4 class \code{a4aFitMCMC}}
\description{
  The \code{a4aFitMCMC} class extends \code{a4aFit} and was
  built to store the MCMC simulations of an a4a stock
  assessment fit.
}
\section{Slots}{
  \describe{ \item{pars}{Matrix with model parameters}

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

