\name{geta4aLatticeOptions}
\alias{geta4aLatticeOptions}
\title{Adds the location of the a4a executable to the search path.  This function is called when FLa4a is attached}
\usage{
  geta4aLatticeOptions(...)
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

