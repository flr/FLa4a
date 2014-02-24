\name{a4aGr}
\alias{a4aGr}
\alias{grInvMod}
\alias{grInvMod<-}
\alias{grMod}
\alias{grMod<-}
\title{Individual growth class}
\description{
  Class definition (slots), constructors, accessors,
  replacement (when relevant) and common methods.
}
\section{Slot}{
  \describe{ \item{\code{grMod}}{the formula of the growth
  model, \emph{e.g.} von bertallanffy}
  \item{\code{grInvMod}}{the formula of the inverse of the
  growth model having length as the independent variable}
  \item{\code{params}}{a FLPar object with the parameters
  of the model. Must match the equations in the models}
  \item{\code{vcov}}{an array with the variance covariance
  matrix of the parameters} \item{\code{distr}}{a character
  with the parameters statistical distribution, it must
  match a known distribution for R, \emph{e.g.} "norm" for
  gaussian, so that \code{rnorm} can be called} }
}

\section{Constructor}{
  \describe{ \item{grMod}{a \code{formula} with the growth
  model (length~age)} \item{grInvMod}{a \code{formula} with
  the inverse growth model (age~length)} \item{params}{a
  \code{FLPar} object with the parameters of the models}
  \item{vcov}{a \code{array} with the variance covariance
  matrix of the parameters} \item{distr}{a \code{character}
  with the distribution of the parameters} }
}
\examples{
mm <- matrix(NA, ncol=3, nrow=3)
diag(mm) <- c(50, 0.001,0.001)
mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","ano-1","ano")), vcov=mm, distr="norm")
}

