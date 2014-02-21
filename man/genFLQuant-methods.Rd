\docType{methods}
\name{genFLQuant}
\alias{genFLQuant}
\alias{genFLQuant,FLQuant-method}
\title{Methods to genetate FLStock objects}
\arguments{
  \item{object}{an FLQuant}

  \item{rec}{an FLQuant}

  \item{catch.n}{an FLQuant}

  \item{harvest}{an FLQuant}

  \item{...}{Additional argument list that might not ever
  be used.}
}
\value{
  an FLQuant
}
\description{
  Some additional details about this S4 generic and its
  methods. The extra blank line between this section and
  the title is critical for roxygen2 to differentiate the
  title from the description section.
}
\examples{
data(ple4)
sim.F <- genFLQuant(harvest(ple4), method = "ac")
}
\seealso{
  \code{\link{print}} and \code{\link{cat}}
}

