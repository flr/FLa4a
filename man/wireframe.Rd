% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\docType{methods}
\name{wireframe plot for FLQuant}
\alias{wireframe plot for FLQuant}
\alias{wireframe,FLQuant-method}
\alias{wireframe,FLQuant,missing-method}
\title{wireframe plot for FLQuant}
\usage{
\S4method{wireframe}{FLQuant,missing}(x, y, screen = list(x = -90, y = -45), ...)
}
\arguments{
\item{x}{a \code{FLQuant}}

\item{y}{missing}

\item{screen}{list with numeric components 'x','y' and 'z' to change the 3D perspective}

\item{...}{additional argument list for the lattice engine}
}
\value{
a 3D surface plot
}
\description{
Method to 3D plot \code{FLQuant} objects.
}
\examples{
data(ple4)
wireframe(harvest(ple4))
}
