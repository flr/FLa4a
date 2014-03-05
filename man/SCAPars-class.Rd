\docType{class}
\name{SCAPars}
\alias{fCovar}
\alias{fCovar-methods}
\alias{fCovar,SCAPars-method}
\alias{fFrml}
\alias{fFrml-methods}
\alias{fFrml,SCAPars-method}
\alias{fmodel}
\alias{fmodel-methods}
\alias{fmodel,SCAPars-method}
\alias{fPars}
\alias{fPars-methods}
\alias{fPars,SCAPars-method}
\alias{n1model}
\alias{n1model-methods}
\alias{n1model,SCAPars-method}
\alias{qCovar}
\alias{qCovar-methods}
\alias{qCovar,SCAPars-method}
\alias{qFrml}
\alias{qFrml-methods}
\alias{qFrml,SCAPars-method}
\alias{qmodel}
\alias{qmodel-methods}
\alias{qmodel,SCAPars-method}
\alias{qPars}
\alias{qPars-methods}
\alias{qPars,SCAPars-method}
\alias{SCAPars}
\alias{SCAPars-class}
\alias{SCAPars-methods}
\alias{SCAPars,missing-method}
\alias{srCovar}
\alias{srCovar-methods}
\alias{srCovar,SCAPars-method}
\alias{srFrml}
\alias{srFrml-methods}
\alias{srFrml,SCAPars-method}
\alias{srmodel}
\alias{srmodel-methods}
\alias{srmodel,SCAPars-method}
\alias{srPars}
\alias{srPars-methods}
\alias{srPars,SCAPars-method}
\alias{stkmodel}
\alias{stkmodel-methods}
\alias{stkmodel,SCAPars-method}
\alias{vCovar}
\alias{vCovar-methods}
\alias{vCovar,SCAPars-method}
\alias{vFrml}
\alias{vFrml-methods}
\alias{vFrml,SCAPars-method}
\alias{vmodel}
\alias{vmodel-methods}
\alias{vmodel,SCAPars-method}
\alias{vPars}
\alias{vPars-methods}
\alias{vPars,SCAPars-method}
\title{Model parameters class}
\description{
  Class definition (slots), constructors, accessors,
  replacement (when relevant) and common methods.
}
\section{Slot}{
  \describe{

  \item{\code{stkmodel}}{Parameters related with the stock
  dynamics}

  \item{\code{qmodel}}{Paramaters related with catchability
  of the tunning fleets}

  \item{\code{vmodel}}{Paramaters related with the variance
  model}

  }
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

