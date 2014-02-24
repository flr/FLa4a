\docType{class}
\name{a4aFitExt-class}
\alias{a4aFitExt}
\alias{a4aFitExt-class}
\alias{baseLvlPars}
\alias{designMatrix}
\alias{L}
\alias{Sigma}
\title{S4 class \code{a4aFitExt}}
\description{
  The \code{a4aFitExt} calss extends \code{a4aFitSA} and
  was built to store extra slots like the hessian, which
  maybe needed for advanced users.
}
\section{Slots}{
  \describe{ \item{Sigma}{The variance-covariance matrix.}

  \item{L}{Likelihood}

  \item{designMatrix}{The design matrix of the model}

  \item{baseLvlPars}{...}

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

