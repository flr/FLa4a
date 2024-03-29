% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/a4aGr-methods.R
\name{mvrtriangle for a4aGr}
\alias{mvrtriangle for a4aGr}
\alias{mvrtriangle,numeric,a4aGr-method}
\title{mvrtriangle}
\usage{
\S4method{mvrtriangle}{numeric,a4aGr}(n = 1, object, ...)
}
\arguments{
\item{n}{the number of iterations}

\item{object}{object of relevant class (see signature of method)}

\item{...}{additional argument list that might never be used}
}
\value{
an \code{a4aGr} object with n iterations
}
\description{
Method to generate multivariate parameters with elliptical copulas and triangular marginals for \code{a4aGr} objects.
}
\details{
The method is essentially a special case of \code{mvrcop}, where the copula is of type "ellipCopula" and family "t", and where the marginals are triangular.
}
\examples{
# Set up the a4aGr object and parameters for the marginals
mm <- matrix(NA, ncol=3, nrow=3)
diag(mm) <- c(50, 0.001,0.001)
mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
md <- ~linf*(1-exp(-k*(t-t0)))
imd <- ~t0-1/k*log(1-len/linf)
prs <- FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","yr^-1","yr"))
vbObj <- a4aGr(grMod=md, grInvMod=imd, params=prs, vcov=mm, distr="norm")
pars <- list(list(a=50, b=100, c=58.5), list(a=0.06, b=0.2, c=0.086), list(a=0, b=0.005, c=0.001))

# Note that mvrtriangle is a special case of mvrcop
set.seed(1)
vbObj1 <- mvrtriangle(10000, vbObj, paramMargins=pars, dispstr="ex", param=0)
set.seed(1)
vbObj2 <- mvrcop(10000, vbObj, copula="ellipCopula", family="t", 
	param=0, margins="triangle", paramMargins=pars)
all.equal(vbObj2, vbObj1)
}
