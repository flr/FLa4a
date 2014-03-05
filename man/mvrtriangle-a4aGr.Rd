\name{mvrtriangle}
\alias{mvrtriangle}
\alias{mvrtriangle,numeric,a4aGr-method}
\title{mvrtriangle}
\arguments{
  \item{n}{the number of simulations to be generated}

  \item{mu}{a \code{a4aGr} object}
}
\value{
  a \code{a4aGr} object with n iterations
}
\description{
  Method to simulate multivariate parameters with triangle
  marginals and elliptic copulas for \code{a4aGr} objects.
}
\examples{
mm <- matrix(NA, ncol=3, nrow=3)
diag(mm) <- c(50, 0.001,0.001)
mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","ano-1","ano")), vcov=mm, distr="norm")
pars <- list(list(a=90, b=125, c=120), list(a=0.2, b=0.4), list(a=0, b=0.4, c=0.1))
vbObj <- mvrtriangle(10000, vbObj, paramMargins=pars)
splom(data.frame(t(params(vbObj)@.Data)), pch=".")
}

