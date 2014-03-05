\name{mvnorm}
\alias{mvnorm}
\alias{mvrnorm,numeric,a4aM,missing,missing,missing,missing-method}
\title{natural mortality}
\arguments{
  \item{n}{the number of simulations to be generated}

  \item{mu}{a \code{a4aM} object}
}
\value{
  a \code{a4aM} object with n iterations
}
\description{
  Method to simulate multivariate normal parameters for a
  \code{a4aM} object.
}
\examples{
mod1 <- FLModelSim(model=~exp(-age-0.5))
mod2 <- FLModelSim(model=~k^0.66*t^0.57, params=FLPar(matrix(c(0.4,10,0.5,11), ncol=2, dimnames=list(params=c("k","t"), iter=1:2))), vcov=array(c(0.004, 0.00,0.00, 0.001), dim=c(2,2,2)))
mod3 <- FLModelSim(model=~1+b*v, params=FLPar(b=0.05))
mObj <- a4aM(shape=mod1, level=mod2, trend=mod3, range=c(min=0,max=15,minyear=2000,maxyear=2003,minmbar=0,maxmbar=0))
mObj <- mvrnorm(100, mObj)
m(mObj, v=c(1,1,1,1))
}

