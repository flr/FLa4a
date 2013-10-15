\name{a4a}
\alias{a4a}
\title{The user interface to the a4a fitting routine.}
\usage{
  a4a(fmodel = ~s(age, k = 3) + factor(year),
    qmodel = lapply(seq(length(indices)), function(i) ~1),
    srmodel = ~factor(year), stock, indices, covar = NULL,
    wkdir = NULL, verbose = FALSE, fit = "MP")
}
\arguments{
  \item{fmodel}{a formula object depicting the model for
  log fishing mortality at age}

  \item{qmodel}{a list of formula objects depicting the
  models for log survey catchability at age}

  \item{rmodel}{a formula object depicting the model for
  log recruitment}

  \item{stock}{an FLStock object containing catch and stock
  information}

  \item{indices}{an FLIndices object containing survey
  indices}

  \item{fmodel.extra}{a named list of expressions used to
  add covariates to the model for log fishing mortality at
  age}

  \item{qmodel.extra}{a named list of expressions used to
  add covariates the models for log survey catchability at
  age}

  \item{wkdir}{used to set a working directory for the admb
  optimiser.  If wkdir is set all admb files are saved to
  this folder otherwise they are deleted.}

  \item{verbose}{if true admb fitting information is
  printed to the screen}
}
\value{
  an FLa4aFit object
}
\description{
  The user interface to the a4a fitting routine.
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}

  Ernesto Jardim \email{ernesto.jardim@jrc.ec.europa.eu}
}

