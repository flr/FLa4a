\docType{methods}
\name{vcov}
\alias{vcov}
\alias{vcov,SCAPars-method}
\alias{vcov,a4aFitSA-method}
\alias{vcov,submodel-method}
\alias{vcov,submodels-method}
\alias{vcov<-,SCAPars,numeric-method}
\alias{vcov<-,a4aFitSA,numeric-method}
\alias{vcov<-,a4aStkParams,numeric-method}
\alias{vcov<-,submodel,numeric-method}
\alias{vcov<-,submodels,numeric-method}
\title{variance-covariance matrix}
\usage{
\S4method{vcov}{a4aFitSA}(object)

\S4method{vcov}{SCAPars}(object)

\S4method{vcov}{submodels}(object)

\S4method{vcov}{submodel}(object)

\S4method{vcov}{a4aFitSA,numeric}(object, ...) <- value

\S4method{vcov}{SCAPars,numeric}(object, ...) <- value

\S4method{vcov}{a4aStkParams,numeric}(object, ...) <- value

\S4method{vcov}{submodels,numeric}(object, ...) <- value

\S4method{vcov}{submodel,numeric}(object, ...) <- value
}
\arguments{
  \item{object}{ a fitted model object, typically.
  Sometimes also a \code{\link{summary}()} object of such a
  fitted model.  }

  \item{...}{ additional arguments for method functions.
  For the \code{\link{glm}} method this can be used to pass
  a \code{dispersion} parameter.  }
}
\description{
Methods to extract and replace the variance-covariance
matrix
}

