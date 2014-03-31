\name{getK}
\alias{getK}
\alias{getK,a4aGr-method}
\alias{getK-methods}
\title{get K}
\usage{
getK(object, ...)
}
\arguments{
  \item{object}{a \code{a4aGr} object}
}
\value{
a \code{vector} with K values
}
\description{
Method to get K values
}
\examples{
mm <- matrix(NA, ncol=3, nrow=3)
diag(mm) <- c(50, 0.001,0.001)
mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","ano-1","ano")), vcov=mm, distr="norm")
vbObj <- mvrnorm(100,vbObj)
getK(vbObj)
}

