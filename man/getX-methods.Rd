\name{getX}
\alias{getX}
\alias{getX,formula-method}
\alias{getX-methods}
\title{Get model matrix}
\usage{
getX(object, ...)
}
\arguments{
  \item{formula}{a formula object}

  \item{df}{the data.frame to build the model matrix
  against}
}
\value{
a matrix.
}
\description{
uses the user specified formula to build a model matrix
}
\note{
\code{getX} is intended to be used internally
}

