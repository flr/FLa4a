\docType{methods}
\name{l2a}
\alias{l2a}
\alias{l2a,FLIndex,a4aGr-method}
\alias{l2a,FLQuant,a4aGr-method}
\alias{l2a,FLStockLen,a4aGr-method}
\alias{l2a-methods}
\title{Method to convert age to length}
\usage{
l2a(object, model, ...)

\S4method{l2a}{FLStockLen,a4aGr}(object, model, plusgroup = "missing", ...)

\S4method{l2a}{FLIndex,a4aGr}(object, model, ...)
}
\arguments{
  \item{object}{a \code{FLQuant} object}

  \item{model}{a \code{a4aGr} object}

  \item{stat}{the aggregation statistic, must be \"mean\"
  or \"sum\"}
}
\value{
a \code{FLQuant} object
}
\description{
Method to convert age to length
}
\examples{
# red fish
# M=0.05; Linf=58.5, k=0.086
mm <- matrix(NA, ncol=3, nrow=3)
diag(mm) <- c(50, 0.001,0.001)
mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","ano-1","ano")), vcov=mm, distr="norm")
data(rfLen)
cth <- catch.n(rfLen.stk)
# both with iter=1
cthA1 <- l2a(cth, vbObj)
# both with iter=n
cthA2 <- l2a(propagate(cth,10), mvrnorm(10, vbObj))
# mod: iter=n, data: iter=1
cthA3 <- l2a(cth, mvrnorm(10, vbObj))
# mod: iter=1, data: iter=n
cthA4 <- l2a(propagate(cth,10), vbObj)
# converting a stock object
rfAge.stk <- l2a(rfLen.stk, vbObj)
rfAge.stk <- l2a(rfLen.stk, mvrnorm(10, vbObj))
rfAge.stk <- l2a(propagate(rfLen.stk, 10), vbObj)
# converting a index object
rfAge.idx <- l2a(rfTrawl.idx, vbObj)
rfAge.idx <- l2a(rfTrawl.idx, mvrnorm(10, vbObj))
rfAge.idx <- l2a(propagate(rfTrawl.idx, 10), vbObj)
}

