\name{mvrtriangle}
\alias{mvrtriangle}
\title{Simulation with a copula model and triangle distributions}
\arguments{
  \item{n}{the number of iterations}

  \item{object}{the \code{FLModelSim} object}

  \item{...}{arguments to be passed to the rMvdc and copula
  methods}

  \item{n}{the number of simulations to be generated}

  \item{mu}{a \code{a4aGr} object}
}
\value{
  a \code{FLModelSim} object with n groups of parameters

  a \code{a4aGr} object with n iterations
}
\description{
  Simulates the model parameters using elliptic copulas and
  triangle margins.

  Method to simulate multivariate parameters with triangle
  marginals and elliptic copulas
}
\examples{
mm <- matrix(NA, ncol=3, nrow=3)
diag(mm) <- c(100, 0.001,0.001)
mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.1,0.0003)
vb <- FLModelSim(model=~linf*(1-exp(-k*(t-t0))), params=FLPar(linf=120, k=0.3, t0=0.1, units=c("cm","ano-1","ano")), vcov=mm, distr="norm")
par(mfrow=c(3,1))
set.seed(1)
vbSim <- mvrnorm(10000, vb)
mm <- predict(vbSim, t=0:20+0.5)
boxplot(t(mm), main="normal")
set.seed(1)
vbSim1 <- mvrtriangle(10000, vb)
mm1 <- predict(vbSim1, t=0:20+0.5)
boxplot(t(mm1), main="triangle")
set.seed(1)
pars <- list(list(a=90, b=125, c=120), list(a=0.2, b=0.4), list(a=0, b=0.4, c=0.1))
vbSim2 <- mvrtriangle(10000, vb, paramMargins=pars)
mm2 <- predict(vbSim2, t=0:20+0.5)
boxplot(t(mm2), main="triangle2")
splom(data.frame(t(params(vbSim)@.Data)), pch=".")
splom(data.frame(t(params(vbSim1)@.Data)), pch=".")
splom(data.frame(t(params(vbSim2)@.Data)), pch=".")
mm <- matrix(NA, ncol=3, nrow=3)
diag(mm) <- c(50, 0.001,0.001)
mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","ano-1","ano")), vcov=mm, distr="norm")
pars <- list(list(a=90, b=125, c=120), list(a=0.2, b=0.4), list(a=0, b=0.4, c=0.1))
vbObj <- mvrtriangle(10000, vbObj, paramMargins=pars)
splom(data.frame(t(params(vbObj)@.Data)), pch=".")
}
\author{
  EJ \email{ernesto.jardim@jrc.ec.europa.eu}

  EJ \email{ernesto.jardim@jrc.ec.europa.eu}
}

