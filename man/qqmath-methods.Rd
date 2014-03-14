\docType{methods}
\name{qqmath}
\alias{qqmath}
\alias{qqmath,a4aFitResiduals,missing-method}
\title{qqplot of standardized log residuals}
\usage{
qqmath(x, data, ...)
}
\arguments{
  \item{x}{a \code{a4aFitResiduals} object with the
  standardized residuals}

  \item{...}{Additional argument list that might not ever
  be used.}
}
\value{
a \code{qqplot} with stardardized log residuals
}
\description{
Method to plot qqplots of standardized residuals
}
\examples{
data(ple4)
data(ple4.index)
obj <- sca(ple4, FLIndices(ple4.index))
}

