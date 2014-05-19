\docType{class}
\name{a4aStkParams}
\alias{a4aStkParams}
\alias{a4aStkParams,missing-method}
\alias{a4aStkParams-class}
\alias{a4aStkParams-methods}
\alias{distr,a4aStkParams-method}
\alias{distr-methods}
\alias{distr<-,a4aStkParams,character-method}
\alias{distr<--methods}
\alias{fMod}
\alias{fMod,a4aStkParams-method}
\alias{fMod-methods}
\alias{fMod<-}
\alias{fMod<-,a4aStkParams,formula-method}
\alias{fMod<--methods}
\alias{n1Mod}
\alias{n1Mod,a4aStkParams-method}
\alias{n1Mod-methods}
\alias{n1Mod<-}
\alias{n1Mod<-,a4aStkParams,formula-method}
\alias{n1Mod<--methods}
\alias{params,a4aStkParams-method}
\alias{params<-,a4aStkParams,FLPar-method}
\alias{rMod}
\alias{rMod,a4aStkParams-method}
\alias{rMod-methods}
\alias{rMod<-}
\alias{rMod<-,a4aStkParams,formula-method}
\alias{rMod<--methods}
\alias{vcov,a4aStkParams-method}
\alias{vcov<-,a4aStkParams-method}
\title{Stock parameters class}
\description{
  Class definition (slots), constructors, accessors,
  replacement (when relevant) and common methods.
}
\section{Slot}{
  \describe{ \item{\code{fMod}}{F submodel \code{formula}}
  \item{\code{n1Mod}}{First year N \code{formula}}
  \item{\code{srMod}}{Stock-recruitment submodel
  \code{formula}} \item{\code{params}}{\code{FLPar} with
  parameters} \item{\code{vcov}}{\code{array} with
  variance-covariance} \item{\code{centering}}{Centering
  values \code{numeric}} \item{\code{distr}}{Statistical
  distribution \code{character}} \item{\code{m}}{Natural
  mortality \code{FLQuant}} \item{\code{units}}{Data units
  \code{character}} }
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

