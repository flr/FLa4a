\name{a4aM}
\alias{a4aM}
\alias{a4aM-class}
\alias{a4aM-methods}
\alias{a4aM,missing-method}
\alias{level}
\alias{level-methods}
\alias{level,a4aM-method}
\alias{level<-}
\alias{level<--methods}
\alias{level<-,a4aM-method}
\alias{shape}
\alias{shape-methods}
\alias{shape,a4aM-method}
\alias{shape<-}
\alias{shape<--methods}
\alias{shape<-,a4aM-method}
\alias{show}
\alias{show,a4aM-method}
\alias{trend}
\alias{trend-methods}
\alias{trend,a4aM-method}
\alias{trend<-}
\alias{trend<--methods}
\alias{trend<-,a4aM-method}
\title{Natural mortality class}
\description{
  Class definition (slots), constructors, accessors,
  replacement (when relevant) and common methods.
}
\section{Slot}{
  \describe{ \item{\code{shape}}{the shape of M by age}
  \item{\code{level}}{the mean level of M in a range of
  ages, which will be used to scale the \code{shape}}
  \item{\code{trend}}{the yearly trend in M} }
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
mod1 <- FLModelSim(model=~exp(-age-0.5))
mod2 <- FLModelSim(model=~1.5*k, params=FLPar(k=0.4))
m1 <- a4aM(shape=mod1, level=mod2)
}

