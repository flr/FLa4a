\docType{methods}
\name{simulate}
\alias{simulate}
\alias{simulate,SCAPars-method}
\alias{simulate,a4aFitSA-method}
\alias{simulate,a4aStkParams-method}
\alias{simulate,submodel-method}
\alias{simulate,submodels-method}
\alias{simulate-methods}
\title{Simulation methods for SCA}
\usage{
simulate(object, nsim = 1, seed = NULL, ...)

\S4method{simulate}{a4aFitSA}(object, nsim = 1, iter = NULL)

\S4method{simulate}{SCAPars}(object, nsim = 1, iter = NULL)

\S4method{simulate}{a4aStkParams}(object, nsim = 1, iter = NULL)

\S4method{simulate}{submodels}(object, nsim = 1, iter = NULL)

\S4method{simulate}{submodel}(object, nsim = 1, iter = NULL)
}
\arguments{
  \item{object}{an object representing a fitted model.}

  \item{nsim}{number of response vectors to simulate.
  Defaults to \code{1}.}

  \item{seed}{an object specifying if and how the random
  number generator should be initialized
  (\sQuote{seeded}).\cr For the "lm" method, either
  \code{NULL} or an integer that will be used in a call to
  \code{set.seed} before simulating the response vectors.
  If set, the value is saved as the \code{"seed"} attribute
  of the returned value.  The default, \code{NULL} will not
  change the random generator state, and return
  \code{\link{.Random.seed}} as the \code{"seed"}
  attribute, see \sQuote{Value}.  }

  \item{...}{additional optional arguments.}
}
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

