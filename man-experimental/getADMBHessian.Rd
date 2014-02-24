\name{getADMBHessian}
\alias{getADMBHessian}
\title{Reads the hessian file from any ADMB fit.  Used here with the a4a model.}
\usage{
  getADMBHessian(wkdir)
}
\arguments{
  \item{wkdir}{the location of the admb output}
}
\value{
  a list with the following elements
}
\description{
  Reads the hessian file from any ADMB fit.  Used here with
  the a4a model.
}
\note{
  \code{getADMBHessian} is intended to be used internally
}
\examples{
# load some data
data(ple4)
data(ple4.indices)
# choose a working directory
wkdir <- tempfile()
# do an 'assessment' fit wth default settings (not recomended!) and keep results in wkdir
fit <- a4a(stock = ple4, indices = ple4.indices, wkdir = wkdir, fit = "assessment")

hessInfo <- getADMBHessian(wkdir)
str(hessInfo)

# calculate covariance matrix
Sigma <- solve(hessInfo $ hes)

# plot correlation matrix of parameters
Cor <- cov2cor(Sigma)
diag(Cor) <- 0
colors <- colorRampPalette(c("blue","grey99","red"))(100)
image(Matrix(Cor), main = "correlations between parameter estimates", lwd = 0, col.regions = colors)
# end of example
}
\author{
  not Colin Millar \email{colin.millar@jrc.ec.europa.eu}
}

