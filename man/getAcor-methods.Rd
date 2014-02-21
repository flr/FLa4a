\docType{methods}
\name{getAcor}
\alias{getAcor}
\alias{getAcor,FLQuant-method}
\title{Methods to add uncertainty to FLQuant objects}
\arguments{
  \item{object}{an FLQuant}

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
genFLQuant(harvest(ple4), method = "ac")
}
\seealso{
  \code{\link{print}} and \code{\link{cat}}
}

