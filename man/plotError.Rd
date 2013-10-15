\name{plotError}
\alias{plotError}
\title{Extracts the catch and survey data, makes useful summaries
and places them in an environment}
\usage{
  plotError(x, y, y.se, FUN = function(x) x, ylab = "",
    xlab = "", add = FALSE, cols)
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

