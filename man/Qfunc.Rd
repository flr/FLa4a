\name{Qfunc}
\alias{Qfunc}
\title{Makes a precision matrix for a specific GMRF model}
\usage{
  Qfunc(n, type = "rw1", weights = 1)
}
\arguments{
  \item{n}{the size of the GMRF}

  \item{type}{a character giving the name of a GMRF (only
  rw1 and rw2 are currently available)}

  \item{weights}{numeric to supply weights if GMRF is to
  have a different variance through the feild}
}
\value{
  the covariance matrix of a GMRF
}
\description{
  Makes a precision matrix for a specific GMRF model
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

