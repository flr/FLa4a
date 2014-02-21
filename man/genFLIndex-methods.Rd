\docType{methods}
\name{genFLIndex}
\alias{genFLIndex}
\alias{genFLIndex,FLQuant-method}
\title{Methods to create index FLQuants from stock.n}
\arguments{
  \item{object}{an FLQuant containing stock numbers at age}

  \item{...}{Additional argument list that might not ever
  be used.}
}
\value{
  an FLQuant?
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
  \code{\link{stock.n}} and \code{\link{FLQuant}}
}

