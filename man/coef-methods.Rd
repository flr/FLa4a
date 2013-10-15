\docType{methods}
\name{coef}
\alias{coef}
\alias{coef,FLa4aFit-method}
\title{Calculate the median accross iterations}
\arguments{
  \item{object}{an FLQuant with iters}

  \item{...}{Additional argument list that might not ever
  be used.}
}
\value{
  an FLQuant
}
\description{
  Calculate the median accross iterations
}
\examples{
data(ple4)
genFLQuant(harvest(ple4), method = "ac")
}
\seealso{
  \code{\link{print}} and \code{\link{cat}}
}

