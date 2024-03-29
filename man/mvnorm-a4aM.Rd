% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/a4aM-methods.R
\name{mvnorm}
\alias{mvnorm}
\alias{mvrnorm,numeric,a4aM,missing,missing,missing,missing-method}
\title{natural mortality}
\usage{
\S4method{mvrnorm}{numeric,a4aM,missing,missing,missing,missing}(n = 1, mu)
}
\arguments{
\item{n}{the number of iterations to be generated}

\item{mu}{an \code{a4aM} object}
}
\value{
an \code{a4aM} object with n iterations
}
\description{
Method to simulate multivariate normal parameters for an \code{a4aM} object.
}
\examples{
mod1 <- FLModelSim(model=~exp(-age-0.5))
mod2 <- FLModelSim(model=~k^0.66*t^0.57, params=FLPar(matrix(c(0.4,10,0.5,11),
 ncol=2, dimnames=list(params=c("k","t"), iter=1:2))),
 vcov=array(c(0.004,0.,0.,0.001,0.006,0.,0.,0.003), dim=c(2,2,2)))
mod3 <- FLModelSim(model=~1+b*v, params=FLPar(b=0.05))
mObj <- a4aM(shape=mod1, level=mod2, trend=mod3,
 range=c(min=0,max=15,minyear=2000,maxyear=2003,minmbar=0,maxmbar=0))
mObj <- mvrnorm(100, mObj)
# Generate 100 iterations with no trend over time
  m(mObj, v=c(1,1,1,1))
# Generate replicates based on iteration-specific multivariate distributions
# (as defined by params() and vcov())
  params(mod2)
  vcov(mod2)
  m1<-mvrnorm(mod2)
  c(params(m1))
# Generate replicates based on a single multivariate distribution (here the
# median of params() and vcov() is used)
  mvrnorm(2,mod2)
  m2<-mvrnorm(2,mod2)
  c(params(m2))
}
