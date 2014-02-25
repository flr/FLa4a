\docType{methods}
\name{genFLStock}
\alias{genFLStock}
\alias{genFLStock-methods}
\alias{genFLStock,FLStock,FLQuant,FLQuant,missing-method}
\alias{genFLStock,FLStock,FLQuant,missing,FLQuant-method}
\alias{genFLStock,FLStock,missing,FLQuant,FLQuant-method}
\title{Methods to generate FLStock objects}
\arguments{
  \item{object}{an FLStock}

  \item{R}{an FLQuant with iterations or missing}

  \item{C}{an FLQuant with iterations or missing}

  \item{F}{an FLQuant with iterations or missing}

  \item{...}{Additional argument list that might not ever
  be used.}
}
\value{
  an FLStock
}
\description{
  This method computes the \code{FLStock} slots
  consistently with the information provided by the
  \code{FLQuant}. It requires two of the triplet R/C/F to
  compute the third consistent with Baranov and survival's
  equations.
}

