\docType{methods}
\name{a4aInternal}
\alias{a4aInternal}
\title{stock assessment model advanced method}
\usage{
a4aInternal(stock, indices, fmodel = ~s(age, k = 3) + factor(year),
  qmodel = lapply(seq(length(indices)), function(i) ~1),
  srmodel = ~factor(year), n1model = ~factor(age),
  vmodel = lapply(seq(length(indices) + 1), function(i) ~1),
  covar = missing, wkdir = missing, verbose = FALSE, fit = "assessment",
  niters = 1000, center = TRUE)
}
\arguments{
  \item{fmodel}{a formula object depicting the model for
  log fishing mortality at age}

  \item{qmodel}{a list of formula objects depicting the
  models for log survey catchability at age}

  \item{srmodel}{a formula object depicting the model for
  log recruitment}

  \item{n1model}{a formula object depicting the model for
  the first year of catch data}

  \item{vmodel}{a list of formula objects depicting the
  models for log survey and log fishing mortality variance}

  \item{stock}{an FLStock object containing catch and stock
  information}

  \item{indices}{an FLIndices object containing survey
  indices}

  \item{covar}{a list with covariates}

  \item{wkdir}{used to set a working directory for the admb
  optimiser.  If wkdir is set all admb files are saved to
  this folder otherwise they are deleted.}

  \item{verbose}{if true admb fitting information is
  printed to the screen}

  \item{fit}{Character with type of fit, 'MP' or
  'assessment', the first doesn't require the hessian to be
  computed, while the former does.}

  \item{niters}{number of iterations to be simulated}

  \item{center}{\code{logical} should data be centered
  before estimating or not}
}
\value{
an \code{a4aFit} object if fit is "MP" or an
\code{a4aFitSA} if fit is "assessment"
}
\description{
The advanced user interface to the a4a fitting routine.
}

