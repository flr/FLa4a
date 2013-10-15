\name{getYidx}
\alias{getYidx}
\title{Gets the FLQuant's numeric id for a vector of "years".
For internal use, not very interesting for users}
\arguments{
  \item{object}{a \code{FLQuant} object}

  \item{year}{a \code{vector} of years}
}
\value{
  \code{numeric vector} that can be used to subset the
  \code{FLQuant}
}
\description{
  Gets the FLQuant's numeric id for a vector of "years".
  For internal use, not very interesting for users
}
\examples{
data(ple4)
flq <- catch(ple4)
getYidx(flq, 2000:2004)
flq[, getYidx(flq, 2000:2004)]
}
\author{
  EJ \email{ernesto.jardim@jrc.ec.europa.eu}
}

