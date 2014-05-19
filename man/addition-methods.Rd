\name{addition}
\alias{+,FLIndices,a4aFit-method}
\alias{+,FLStock,a4aFit-method}
\alias{+,FLStock,a4aFitSA-method}
\title{Addition methods}
\description{
  Update \code{FLStock} and \code{FLIndex} objects with
  stock assessment results
}
\details{
  If both objects have the same number of iterations, the
  \code{FLStock} slots will be replaced by the
  \code{a4aFit} slots, in case of 1 iter, or
  \code{a4aFitSA} slots, in case of n iters. If one of the
  objects has 1 iter and the other n, the method will
  simulate using the fit results from the \code{a4aFitSA}
  object, to update the slots of the \code{FLStock} object.
}

