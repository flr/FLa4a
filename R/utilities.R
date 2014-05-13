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

#' @rdname assorted-methods
#' @aliases dims,a4aStkParams-method
#' @section dims:
#' Extracts the dims of the parameters.
setMethod("dims", "a4aStkParams", function(obj) {
  dim(obj@params)
})
