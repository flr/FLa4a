\name{a4aM}
\alias{a4aM}
\title{Natural mortality class}
\arguments{
  \item{shape}{a \code{FLModelSim} object}

  \item{level}{a \code{FLModelSim} object}

  \item{trend}{a \code{FLModelSim} object}
}
\value{
  a \code{a4aM} object
}
\description{
  Natural mortality class

  Constructor for the natural mortality class using
  \code{FLModelSim} objects.
}
\section{Slot}{
  \itemize{ \item \code{shape} the shape of M by age.
  \item \code{level} the mean level of M in a range of
  ages, which will be used to scale the \code{shape}.
  \item \code{trend} the yearly trend in M. '}
}
\examples{
mod1 <- FLModelSim(model=~exp(-age-0.5))
mod2 <- FLModelSim(model=~1.5*k, params=FLPar(k=0.4))
m1 <- a4aM(shape=mod1, level=mod2)
}
\author{
  EJ \email{ernesto.jardim@jrc.ec.europa.eu}

  EJ \email{ernesto.jardim@jrc.ec.europa.eu}
}

