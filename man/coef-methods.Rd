\docType{methods}
\name{coef}
\alias{coef}
\alias{coef,FLa4aFit-method}
\alias{coef,FLa4aFitSA-method}
\alias{coef<-}
\alias{predict}
\alias{simulate}
\alias{vcov}
\alias{vcov<-}
\title{Calculate the median accross iterations}
\arguments{
  \item{object}{an FLQuant with iters}

  \item{...}{Additional argument list that might not ever
  be used.}

  \item{object}{an FLQuant with iters}

  \item{...}{Additional argument list that might not ever
  be used.}
}
\value{
  an FLQuant

  an FLQuant
}
\description{
  Calculate the median accross iterations

  Calculate the median accross iterations
}
\examples{
data(ple4)
genFLQuant(harvest(ple4), method = "ac")
data(ple4)
genFLQuant(harvest(ple4), method = "ac")
}
\seealso{
  \code{\link{print}} and \code{\link{cat}}

  \code{\link{print}} and \code{\link{cat}}
}

