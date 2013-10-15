\name{plotIters}
\alias{plotIters}
\title{Adds the location of the a4a executable to the search path.  This function is called when FLa4a is attached}
\usage{
  plotIters(object, nsamples = 5, ylab = "", by = "quant",
    col = 1, ...)
}
\arguments{
  \item{libname}{required}

  \item{pkgname}{required}
}
\value{
  a pointer to the environment in which summaries of the
  data reside
}
\description{
  Adds the location of the a4a executable to the search
  path.  This function is called when FLa4a is attached
}
\note{
  \code{.onAttach} is intended to be used internally
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

