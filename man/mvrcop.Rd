\name{mvrcop}
\alias{mvrcop}
\title{Simulation using copula models}
\arguments{
  \item{n}{the number of iterations}

  \item{object}{the \code{FLModelSim} object}

  \item{...}{arguments to be passed to the rMvdc and copula
  methods}

  \item{n}{the number of iterations}

  \item{object}{the \code{FLModelSim} object}

  \item{...}{arguments to be passed to the rMvdc and copula
  methods}
}
\value{
  a \code{FLModelSim} object with n groups of parameters

  a \code{FLModelSim} object with n groups of parameters
}
\description{
  Simulates the model parameters using a self defined
  copula model and margins.

  Simulates the model parameters using self defined copulas
  and margins
}
\examples{
mm <- matrix(NA, ncol=3, nrow=3)
diag(mm) <- c(100, 0.001,0.001)
mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.1,0.0003)
vb <- FLModelSim(model=~linf*(1-exp(-k*(t-t0))), params=FLPar(linf=120, k=0.3, t0=0.1, units=c("cm","ano-1","ano")), vcov=mm, distr="norm")
pars <- list(list(a=90, b=125, c=120), list(a=0.2, b=0.4), list(a=0, b=0.4, c=0.1))
vbSim <- mvrcop(10000, vb, copula="archmCopula", family="clayton", param=2, margins="triangle", paramMargins=pars)
boxplot(t(predict(vbSim, t=0:20+0.5)))
splom(data.frame(t(params(vbSim)@.Data)), pch=".")
mm <- matrix(NA, ncol=3, nrow=3)
diag(mm) <- c(50, 0.001,0.001)
mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","ano-1","ano")), vcov=mm, distr="norm")
pars <- list(list(a=90, b=125, c=120), list(a=0.2, b=0.4), list(a=0, b=0.4, c=0.1))
vbObj <- mvrcop(10000, vbObj, copula="archmCopula", family="clayton", param=2, margins="triangle", paramMargins=pars)
splom(data.frame(t(params(vbObj)@.Data)), pch=".")
}
\author{
  EJ \email{ernesto.jardim@jrc.ec.europa.eu}

  EJ \email{ernesto.jardim@jrc.ec.europa.eu}
}
\seealso{
  \code{\link{rMvdc}}, \code{\link{copula}}
}

