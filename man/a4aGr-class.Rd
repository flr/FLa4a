\docType{class}
\name{m}
\alias{a4aGr}
\alias{a4aGr-class}
\alias{a4aGr-methods}
\alias{a4aGr,missing-method}
\alias{distr}
\alias{distr-methods}
\alias{distr,a4aGr-method}
\alias{distr<-}
\alias{distr<--methods}
\alias{distr<-,a4aGr,character-method}
\alias{grInvMod}
\alias{grInvMod-methods}
\alias{grInvMod,a4aGr-method}
\alias{grInvMod<-}
\alias{grInvMod<--methods}
\alias{grInvMod<-,a4aGr,formula-method}
\alias{grMod}
\alias{grMod-methods}
\alias{grMod,a4aGr-method}
\alias{grMod<-}
\alias{grMod<--methods}
\alias{grMod<-,a4aGr,formula-method}
\alias{m}
\alias{m,a4aStkParams-method}
\alias{params}
\alias{params,a4aGr-method}
\alias{params<-}
\alias{params<--methods}
\alias{params<-,a4aGr,FLPar-method}
\alias{vcov}
\alias{vcov,a4aGr-method}
\alias{vcov<-}
\alias{vcov<-,a4aGr,numeric-method}
\title{Individual growth class}
\description{
  Class definition (slots), constructors, accessors,
  replacement (when relevant) and common methods.
}
\section{Slot}{
  \describe{

  \item{\code{grMod}}{the formula of the growth model,
  \emph{e.g.} von bertallanffy}

  \item{\code{grInvMod}}{the formula of the inverse of the
  growth model having length as the independent variable}

  \item{\code{params}}{a FLPar object with the parameters
  of the model. Must match the equations in the models}

  \item{\code{vcov}}{an array with the variance covariance
  matrix of the parameters}

  \item{\code{distr}}{a character with the parameters
  statistical distribution, it must match a known
  distribution for R, \emph{e.g.} "norm" for gaussian, so
  that \code{rnorm} can be called}

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
mm <- matrix(NA, ncol=3, nrow=3)
diag(mm) <- c(50, 0.001,0.001)
mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","ano-1","ano")), vcov=mm, distr="norm")
}

