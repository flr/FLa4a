\name{a4aInternal}
\alias{a4aInternal}
\title{The user interface to the a4a fitting routine.}
\usage{
  a4aInternal(fmodel = ~s(age, k = 3) + factor(year),
    qmodel = lapply(seq(length(indices)), function(i) ~1),
    srmodel = ~factor(year), n1model = ~factor(age),
    vmodel = lapply(seq(length(indices) + 1), function(i) ~1),
    stock, indices, covar = NULL, wkdir = NULL,
    verbose = FALSE, fit = "assessment", niters = 1000)
}
\arguments{
  \item{fmodel}{a formula object depicting the model for
  log fishing mortality at age}

  \item{qmodel}{a list of formula objects depicting the
  models for log survey catchability at age}

  \item{rmodel}{a formula object depicting the model for
  log recruitment}

  \item{stock}{an FLStock object containing catch and stock
  information}

  \item{indices}{an FLIndices object containing survey
  indices}

  \item{fmodel.extra}{a named list of expressions used to
  add covariates to the model for log fishing mortality at
  age}

  \item{qmodel.extra}{a named list of expressions used to
  add covariates the models for log survey catchability at
  age}

  \item{wkdir}{used to set a working directory for the admb
  optimiser.  If wkdir is set all admb files are saved to
  this folder otherwise they are deleted.}

  \item{verbose}{if true admb fitting information is
  printed to the screen}
}
\value{
  an FLa4aFit object
}
\description{
  The user interface to the a4a fitting routine.
}
\examples{
data(ple4)
data(ple4.indices)

# define sub models
fmodel <- ~ s(age, k=4) + factor(year)
# three surveys - try 'length(ple4.indices)' - so three formulas
qmodel <- list(~ s(age, k=4), ~ s(age, k=4), ~ age)
fit <- a4a(fmodel = fmodel, qmodel = qmodel, stock = ple4, indices = ple4.indices, fit = "assessment")

# add the fit to the stock, but with fitting error
fitstk <- propagate(ple4, 1000) + fit

doOne <- function(p, object) cbind(as.data.frame(apply(object, 1:5, quantile, p)), p = p[1])
dat <- do.call(rbind, lapply(c(0.025, 0.5, 0.975), doOne, object = ssb(fitstk)))
xyplot(I(data/1000) ~ year, group = p, data = dat,
       type = c("l","g"), lwd = c(1,2,1), lty = c(2,1,2), col = 1,
       ylab = "SSB (thousand tonnes)", xlab = "Year")
# end of example
}
\author{
  Colin Millar \email{colin.millar@jrc.ec.europa.eu}

  Ernesto Jardim \email{ernesto.jardim@jrc.ec.europa.eu}
}

