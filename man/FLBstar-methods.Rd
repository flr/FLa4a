\docType{methods}
\name{FLblockboot}
\alias{FLblockboot}
\alias{FLBstar}
\alias{FLBstar-methods}
\alias{FLBstar,FLQuant-method}
\title{Block bootstrap method}
\arguments{
  \item{object}{an FLQuant}
}
\value{
  an FLQuant containing the optimal block lengths for each
  year
}
\description{
  Calculate the optimal block length for stationary block
  bootstrap
}
\examples{
data(ple4)
FLBstar(catch.n(ple4))
}

