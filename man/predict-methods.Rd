\docType{methods}
\name{predict}
\alias{predict}
\alias{predict,SCAPars-method}
\alias{predict,a4aFitSA-method}
\alias{predict,a4aStkParams-method}
\alias{predict,submodel-method}
\alias{predict,submodels-method}
\title{predict methods for SCA}
\usage{
\S4method{predict}{a4aFitSA}(object)

\S4method{predict}{SCAPars}(object)

\S4method{predict}{a4aStkParams}(object)

\S4method{predict}{submodels}(object, ...)

\S4method{predict}{submodel}(object, ...)
}
\arguments{
  \item{object}{a model object for which prediction is
  desired.}

  \item{...}{additional arguments affecting the predictions
  produced.}
}
\description{
Predict methods for a4a stock assessment fits.
}
\examples{
data(ple4)
data(ple4.index)
fmodel <- ~factor(age) + factor(year)
qmodel <- list(~factor(age))
fit1 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
}

