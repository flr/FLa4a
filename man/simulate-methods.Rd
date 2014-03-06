\name{simulate}
\alias{simulate}
\alias{simulate-methods}
\alias{simulate,a4aFitSA-method}
\alias{simulate,a4aStkParams-method}
\alias{simulate,SCAPars-method}
\alias{simulate,submodel-method}
\alias{simulate,submodels-method}
\title{Simulation methods for SCA}
\description{
  Simulation methods for a4a stock assessment fits.
}
\examples{
data(ple4)
data(ple4.index)
fmodel <- ~factor(age) + factor(year)
qmodel <- list(~factor(age))
fit1 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
}

