\docType{class}
\name{submodel}
\alias{covar}
\alias{covar,submodel-method}
\alias{model}
\alias{model,submodel-method}
\alias{pars}
\alias{pars-methods}
\alias{pars,submodel-method}
\alias{submodel}
\alias{submodel-class}
\alias{submodel-methods}
\alias{submodel,missing-method}
\title{Submodel class}
\description{
  Class definition (slots), constructors, accessors,
  replacement (when relevant) and common methods.
}
\section{Slot}{
  \describe{

  \item{\code{Mod}}{\code{formula} describing the model}

  \item{\code{params}}{\code{FLPar} with model parameters}

  \item{\code{vcov}}{\code{array} with variance covariance
  Paramaters related with the variance model}

  \item{\code{centering}}{\code{numeric} value used for
  centering the data}

  \item{\code{distr}}{a character with the parameters
  statistical distribution, it must match a known
  distribution for R, \emph{e.g.} "norm" for gaussian, so
  that \code{rnorm} can be called} }
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

