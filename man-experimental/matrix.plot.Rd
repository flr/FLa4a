\name{matrix.plot}
\alias{matrix.plot}
\title{Extracts the catch and survey data, makes useful summaries
and places them in an environment}
\usage{
  matrix.plot(x, xlim = 0.5 + c(0, di[2]),
    ylim = 0.5 + c(di[1], 0), xlab = "", ylab = "",
    cols = c("red", "blue"), text = FALSE, xmin = 1,
    ymin = 1, ...)
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
  Extracts the catch and survey data, makes useful
  summaries and places them in an environment
}
\note{
  \code{extractData} is intended to be used internally
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

