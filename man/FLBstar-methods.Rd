\docType{methods}
\name{FLBstar}
\alias{FLBstar}
\alias{FLBstar,FLQuant-method}
\title{Calculate the optimal block length for stationary block bootstrap}
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

