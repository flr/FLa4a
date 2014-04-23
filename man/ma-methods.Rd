\name{ma}
\alias{ma}
\alias{ma,a4aFitSAs-method}
\alias{ma-methods}
\title{Model averaging (experimental)}
\arguments{
  \item{object}{a \code{a4aFits} object with the fits to be
  averaged across}

  \item{stock}{a \code{stock} object with the original data
  used for fitting}

  \item{fun}{a \code{function} that will be used to extract
  the values for weighting, for now it must be "AIC", "BIC"
  or "LogLik"}

  \item{nsim}{a \code{numeric} with the number os
  simulations to be drawn}
}
\value{
  a \code{FLStock} object with iterations defined by
  \code{nsim}
}
\description{
  Method to average across a set of models. This is still
  experimental. Use with care.
}
\examples{
data(ple4)
data(ple4.indices)
f1 <- sca(ple4, ple4.indices, fmodel=~ factor(age) + s(year, k=20), qmodel=list(~ s(age, k = 4), ~ s(age, k = 4), ~ s(age, k = 3)), fit = "assessment")
f2 <- sca(ple4, ple4.indices, fmodel=~ factor(age) + s(year, k=20), qmodel=list(~ s(age, k = 4)+year, ~ s(age, k = 4), ~ s(age, k = 3)), fit = "assessment")
stock.sim <- ma(a4aFits(list(f1, f2)), ple4, AIC, nsim = 1000)
stks <- FLStocks(f1=ple4+f1, f2=ple4+f2, ma=stock.sim)
flqs <- lapply(stks, ssb)
flqs <- lapply(flqs, iterMedians)
xyplot(data~year, groups=qname, data=flqs, type="l")
plot(stks)
}

