\name{a4a.dir}
\alias{a4a.dir}
\title{returns the location on the file system of the ADMB executable}
\usage{
  a4a.dir()
}
\arguments{
  \item{stock}{an FLStock object containing catch and stock
  information}

  \item{index}{an FLIndex object containing survey indices}
}
\value{
  a pointer to the environment in which summaries of the
  data reside
}
\description{
  returns the location on the file system of the ADMB
  executable
}
\note{
  \code{a4a.dir} is intended to be used internally and
  based on the same function in INLA
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

