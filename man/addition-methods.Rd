% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addition-methods.R
\name{addition}
\alias{addition}
\alias{+,FLStock,a4aFit-method}
\alias{+,FLIndices,a4aFit-method}
\alias{+,FLStocks,a4aFits-method}
\title{+ methods}
\usage{
\S4method{+}{FLStock,a4aFit}(e1, e2)

\S4method{+}{FLStocks,a4aFits}(e1, e2)
}
\arguments{
\item{e1}{the original \code{FLStocks} object}

\item{e2}{a \code{a4aFits} object from where the new \code{FLStock} slots will be extracted.}
}
\description{
Update \code{FLStock} and \code{FLIndex} objects with stock assessment results.

Update \code{FLStocks} objects with multiple stock assessment results in a \code{a4aFits}.
}
\details{
If both objects have the same number of iterations, the \code{FLStock} slots will be replaced by the \code{a4aFit} slots, in the case of 1 iter, or \code{a4aFitSA} slots, in the case of n iters. If one of the objects has 1 iter and the other n, the method will simulate using the fit results from the \code{a4aFitSA} object to update the slots of the \code{FLStock} object.

If both objects have the same number of iterations, the \code{FLStocks} slots will be replaced by the \code{a4aFits} slots, in the case of 1 iter, or \code{a4aFitSA} slots, in the case of n iters. If one of the objects has 1 iter and the other n, the method will simulate using the fit results from the \code{a4aFitSA} object to update the slots of the \code{FLStock} object.
}
