#' @title Assorted methods needed by FLa4a
#' @docType methods
#' @name methods
#' @rdname assorted-methods
#' @aliases is.empty
#' @section is.empty:
#' Method \code{is.empty} checks if object is empty. Takes any object returns a logical, \code{TRUE} if object is of length 0.
#' @examples
#' is.empty(list())
#' is.empty(list(a=2))

is.empty <- function(object) {
	length(object) == 0
}

#' @rdname assorted-methods
#' @aliases pars2dim,FLPar-method
#' @section pars2dim:
#' Checks that the name of the second dimension in params is "iter". For internal use, not very interesting for users. It takes a \code{FLPar} object and returns a \code{logical}
#' @examples
#' pars2dim(FLPar())
#' pars2dim(FLPar(array(dim=c(1,1,1))))
setMethod("pars2dim", "FLPar", function(object) {

	dnm <- dimnames(object)
	names(dnm)[2]=="iter"

})

#' @rdname assorted-methods
#' @aliases getYidx getYidx-methods getYidx,FLQuant-method
#' @section getYidx:
#' Gets the FLQuant's numeric id for a vector of "years". For internal use, not very interesting for users. It takes a \code{FLQuant} and a \code{vector} of years and returns a \code{numeric vector} that can be used to subset the \code{FLQuant}.
#' @examples
#' data(ple4)
#' flq <- catch(ple4)
#' getYidx(flq, 2000:2004)
#' flq[, getYidx(flq, 2000:2004)]

setGeneric("getYidx", function(object, ...) standardGeneric("getYidx"))
setMethod("getYidx", "FLQuant", function(object, year) {
	yrs <- dimnames(object)[[2]]
	if(sum(year>0)>0){
		idx <- match(as.character(year), yrs, nomatch=0)
		if(sum(idx)>0){
			idx
		} else {
			year
		}
	} else {
		length(yrs)+year+1
	} 

})

#' @rdname assorted-methods
#' @aliases niters niters-methods niters,FLModelSim-method
#' @section niters:
#' Compute number of iterations. Takes an object of any \code{FLR} class and returns a \code{numeric}.
setGeneric("niters", function(object, ...) standardGeneric("niters"))
setMethod("niters", "FLModelSim", function(object){
	dim(params(object))[2]
})

#' @rdname assorted-methods
#' @aliases geta4aLatticeOptions
#' @section geta4aLatticeOptions:
#' Set lattice options for a4a plots.
#' @examples
#' # set up grapical options
#'  lattice.options(default.theme = geta4aLatticeOptions())
geta4aLatticeOptions <- function(...) 
{
  opts <- standard.theme(color = FALSE)
  opts <- modifyList(opts, list(strip.background = list(col = grey(0.95))))
  opts <- modifyList(opts, list(plot.symbol = list(pch = 16, col = "black", cex = 0.7)))
  opts <- modifyList(opts, list(axis.line = list(col = grey(0.8))))
  opts <- modifyList(opts, list(strip.border = list(col = grey(0.8))))

  opts
}

#' @rdname assorted-methods
#' @aliases plotIters
#' @section plotIters:
#' Plot iterations.
plotIters <- function(object, nsamples = 5, ylab = "", by = "quant", col = 1, ...) {

  by <- match.arg(by, c("year","quant"))
  quant(object) <- "quant"

  x <- setdiff(c("year","quant"), by)

  doOne <- function(p, object) cbind(as.data.frame(apply(object, 1:5, quantile, p)), p = p)

  dat <- do.call(rbind, lapply(c(0.025, 0.5, 0.975), doOne, object = object))

  dat $ by <- dat[[by]]
  dat $ x <- dat[[x]]


  p1 <- xyplot(data ~ x | factor(by), group = p, data = dat,
               type = c("l","g"), scales = list(y = list(relation = "free")), 
               lwd = c(1,2,1), lty = c(2,1,2), col = col,
                ylab = ylab, xlab = x, as.table = TRUE, ...)
  if (nsamples > 0) {
  
    nsamples <- min(nsamples, min(9, dim(object)[6]))
    dat2 <- do.call(rbind, lapply(sample(dims(object) $ iter, nsamples), function(i) cbind(as.data.frame(object[,,,,,i]), p = i)))
    dat2 $ by <- dat2[[by]]
    dat2 $ x <- dat2[[x]]          

    p2 <- xyplot(data ~ x | factor(by), group = p, data = dat2,
                 type = "l", col =  brewer.pal(max(3, nsamples), "Set1"), lty = 1)
               
    p1 + p2
  } else {
    p1
  }
}

#' @title Get ADMB Hessian
#' @name getADMBHessian
#' @rdname getADMBHessian
#' @description Reads the hessian file from any ADMB fit.  Used here with the a4a model.
#' @param wkdir the location of the admb output
#' @return a list with the following elements
#' @note \code{getADMBHessian} is intended to be used internally
#' @aliases getADMBHessian
#' @examples
#' # load some data
#' data(ple4)
#' data(ple4.indices)
#' # choose a working directory
#' wkdir <- tempfile()
#' # do an 'assessment' fit wth default settings (not recomended!) and keep results in wkdir
#' fit <- a4aSCA(stock = ple4, indices = ple4.indices, wkdir = wkdir, fit = "assessment")
#' hessInfo <- getADMBHessian(wkdir)
#' str(hessInfo)
#' # calculate covariance matrix
#' Sigma <- solve(hessInfo $ hes)
#' # plot correlation matrix of parameters
#' Cor <- cov2cor(Sigma)
#' diag(Cor) <- 0
#' colors <- colorRampPalette(c("blue","grey99","red"))(100)
#' image(Matrix(Cor), main = "correlations between parameter estimates", lwd = 0, col.regions = colors)

getADMBHessian <- function(wkdir) {
## This function reads in all of the information contained in the
## admodel.hes file. Some of this is needed for relaxing the covariance
## matrix, and others just need to be recorded and rewritten to file so ADMB
## "sees" what it’s expecting.
  filename <- file(paste0(wkdir,"/admodel.hes"), "rb")
  on.exit(close(filename))
  num.pars <- readBin(filename, "integer", 1)
  hes.vec <- readBin(filename, "numeric", num.pars^2)
  hes <- matrix(hes.vec, ncol=num.pars, nrow=num.pars)
  hybrid_bounded_flag <- readBin(filename, "integer", 1)
  scale <- readBin(filename, "numeric", num.pars)

  list(num.pars = num.pars, hes = hes, hybrid_bounded_flag = hybrid_bounded_flag, scale = scale)
}


