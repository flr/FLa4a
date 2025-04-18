% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fittingFunctions.R
\docType{methods}
\name{scas}
\alias{scas}
\title{Run several stock assessments in a single run}
\usage{
scas(
  stocks,
  indicess,
  fmodel = missing,
  qmodel = missing,
  srmodel = missing,
  n1model = missing,
  vmodel = missing,
  combination.all = FALSE,
  ...
)
}
\arguments{
\item{stocks}{an \code{FLStocks} object containing catch and stock information}

\item{indicess}{a list of \code{FLIndices} objects containing survey indices}

\item{fmodel}{a list of \code{fmodel} objects, each with a formula object depicting the model for log fishing mortality at age}

\item{qmodel}{a list of \code{qmodel} objects, each with a list of formula objects depicting the models for log survey catchability at age}

\item{srmodel}{a list of \code{srmodel} objects, each with a formula object depicting the model for log recruitment}

\item{n1model}{a list of \code{n1model} objects, each with a formula object depicting the model for the first year of catch data}

\item{vmodel}{a list of \code{vmodel} objects, each with a list of formula objects depicting the models for log survey and log fishing mortality variance}

\item{combination.all}{bolean parameter (default is FALSE) to define if a full factorial across all stocks, indices, and submodel is run or just a sequence of runs.}

\item{...}{all other arguments to be passed to \code{sca}}

\item{stock}{an \code{FLStocks} object, each component with a \code{FLStock} object containing catch and stock information}
}
\value{
an \code{a4aFits} or \code{a4aFitSAs} or \code{a4aFitMCMCs} depending on the argument \code{fit}
}
\description{
Internal method to run several stock assessment fits with different stocks, indices and submodels
}
