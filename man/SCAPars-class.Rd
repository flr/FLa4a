\docType{class}
\name{SCAPars}
\alias{SCAPars}
\alias{SCAPars,missing-method}
\alias{SCAPars-class}
\alias{SCAPars-methods}
\alias{fCovar}
\alias{fCovar,SCAPars-method}
\alias{fCovar-methods}
\alias{fFrml}
\alias{fFrml,SCAPars-method}
\alias{fFrml-methods}
\alias{fPars}
\alias{fPars,SCAPars-method}
\alias{fPars-methods}
\alias{fmodel}
\alias{fmodel,SCAPars-method}
\alias{fmodel-methods}
\alias{m}
\alias{m,SCAPars-method}
\alias{n1model}
\alias{n1model,SCAPars-method}
\alias{n1model-methods}
\alias{qCovar}
\alias{qCovar,SCAPars-method}
\alias{qCovar-methods}
\alias{qFrml}
\alias{qFrml,SCAPars-method}
\alias{qFrml-methods}
\alias{qPars}
\alias{qPars,SCAPars-method}
\alias{qPars-methods}
\alias{qmodel}
\alias{qmodel,SCAPars-method}
\alias{qmodel-methods}
\alias{srCovar}
\alias{srCovar,SCAPars-method}
\alias{srCovar-methods}
\alias{srFrml}
\alias{srFrml,SCAPars-method}
\alias{srFrml-methods}
\alias{srPars}
\alias{srPars,SCAPars-method}
\alias{srPars-methods}
\alias{srmodel}
\alias{srmodel,SCAPars-method}
\alias{srmodel-methods}
\alias{stkmodel}
\alias{stkmodel,SCAPars-method}
\alias{stkmodel-methods}
\alias{vCovar}
\alias{vCovar,SCAPars-method}
\alias{vCovar-methods}
\alias{vFrml}
\alias{vFrml,SCAPars-method}
\alias{vFrml-methods}
\alias{vPars}
\alias{vPars,SCAPars-method}
\alias{vPars-methods}
\alias{vmodel}
\alias{vmodel,SCAPars-method}
\alias{vmodel-methods}
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

