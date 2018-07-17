globalVariables("qname")
#' @title Assorted methods needed by FLa4a
#' @docType methods
#' @name assorted methods
#' @description Assorted methods needed by FLa4a
#' @rdname assorted-methods
#' @aliases getYidx getYidx-methods
#' @template bothargs
#' @section getYidx:
#' Gets an FLQuant's numeric id for a vector of "years". For internal use and not very interesting for users. It takes an \code{FLQuant} object and \code{vector} of years and returns a \code{numeric vector} that can be used to subset the \code{FLQuant}.
#' @examples
#' #Example use of getYidx:
#' data(ple4)
#' flq <- catch(ple4)
#' getYidx(flq, 2000:2004)
#' flq[, getYidx(flq, 2000:2004)]

setGeneric("getYidx", function(object, ...) standardGeneric("getYidx"))
#' @rdname assorted-methods
#' @param year \code{numeric} with year to be extracted
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

#' @rdname pars2dim-methods
#' @section pars2dim:
#' Checks that the name of the second dimension in params is "iter". For internal use and not very interesting for users. It takes an \code{FLPar} object and returns a \code{logical}.
#' @examples
#' #Example use of pars2dim:
#' pars2dim(FLPar())
#' pars2dim(FLPar(array(dim=c(1,1,1))))
setMethod("pars2dim", "FLPar", function(object) {

	dnm <- dimnames(object)
	names(dnm)[2]=="iter"

})

#' @rdname assorted-methods
#' @aliases is.empty
#' @section is.empty:
#' Method \code{is.empty} checks if an object is empty. It takes any object and returns a logical, \code{TRUE}, if the object is of length 0.
#' @examples
#' #Example use of is.empty:
#' is.empty(list())
#' is.empty(list(a=2))
is.empty <- function(object) {
	length(object) == 0
}

#' @rdname assorted-methods
#' @aliases niters niters-methods
#' @section niters:
#' Compute number of iterations. Takes an object of any \code{FLR} class and returns a \code{numeric}.
#' @examples
#' #Example use of niters:
#' mm <- matrix(NA, ncol=3, nrow=3)
#' diag(mm) <- c(50, 0.001,0.001)
#' mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
#' md <- ~linf*(1-exp(-k*(t-t0)))
#' imd <- ~t0-1/k*log(1-len/linf)
#' prs <- FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","yr^-1","yr"))
#' vbObj <- a4aGr(grMod=md, grInvMod=imd, params=prs, vcov=mm, distr="norm")
#' # Generate 100 sample sets
#' vbObj <- mvrnorm(100,vbObj)
#' niters(vbObj)

setGeneric("niters", function(object, ...) standardGeneric("niters"))
#' @rdname assorted-methods
setMethod("niters", "FLModelSim", function(object){
	dim(params(object))[2]
})
#' @rdname assorted-methods
setMethod("niters", "a4aGr", function(object){
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
#' # do an 'assessment' fit with default settings (not recomended!) and keep results in wkdir
#' fit <- sca(stock=ple4,indices=ple4.indices,wkdir=wkdir,fit="assessment", useADMB=TRUE)
#' hessInfo <- getADMBHessian(wkdir)
#' str(hessInfo)
#' # calculate covariance matrix
#' Sigma <- solve(hessInfo$hes)

getADMBHessian <- function(wkdir) {
## This function reads in all of the information contained in the
## admodel.hes file. Some of this is needed for relaxing the covariance
## matrix, and others just need to be recorded and rewritten to file so ADMB
## 'sees' what it's expecting.
  filename <- file(paste0(wkdir,"/admodel.hes"), "rb")
  on.exit(close(filename))
  num.pars <- readBin(filename, "integer", 1)
  hes.vec <- readBin(filename, "numeric", num.pars^2)
  hes <- matrix(hes.vec, ncol=num.pars, nrow=num.pars)
  hybrid_bounded_flag <- readBin(filename, "integer", 1)
  scale <- readBin(filename, "numeric", num.pars)

  list(num.pars = num.pars, hes = hes, hybrid_bounded_flag = hybrid_bounded_flag, scale = scale)
}


#' @rdname getADMBHessian
#' @aliases getADMBCovariance
getADMBCovariance <- function(wkdir) {
## This function reads in all of the information contained in the
## admodel.hes file. Some of this is needed for relaxing the covariance
## matrix, and others just need to be recorded and rewritten to file so ADMB
## 'sees' what it's expecting.
  filename <- file(paste0(wkdir,"/admodel.cov"), "rb")
  on.exit(close(filename))
  num.pars <- readBin(filename, "integer", 1)
  cov.vec <- readBin(filename, "numeric", num.pars^2)
  cov <- matrix(cov.vec, ncol=num.pars, nrow=num.pars)
  hybrid_bounded_flag <- readBin(filename, "integer", 1)
  scale <- readBin(filename, "numeric", num.pars)
  list(num.pars = num.pars, cov = cov, hybrid_bounded_flag = hybrid_bounded_flag, scale = scale)
}

#' @rdname assorted-methods
#' @param obj an object
#' @section dims:
#' Extracts the dims of the parameters.
#' @examples
#' #Example use of dims:
#' dims(FLPar())

setMethod("dims", "a4aStkParams", function(obj) {
  dim(obj@coefficients)
})

#' @title plot for fitted catch-at-age
#' @name plot for fitted catch-at-age
#' @docType methods
#' @rdname plotc
#' @aliases plot,a4aFit,FLStock-method
#' @description Method to plot fitted versus observed catch numbers-at-age. Note the yaxis doesn't has a scale. The visual is about the difference between the two lines, not about the value of each line, which in any case would be very difficult to assess visually.
#' @param x an \code{a4aFit} object with the fitted values
#' @param y an \code{FLStock} object with the observed values
#' @param ... additional argument list that might never be used
#' @return a \code{plot} with fitted and observed catch numbers-at-age
#' @examples
#' data(ple4)
#' data(ple4.index)
#' obj <- sca(ple4, FLIndices(ple4.index))
#' plot(obj, ple4)

setMethod("plot", c("a4aFit", "FLStock"), function(x, y, ...){
	args <- list()
	args$data <- as.data.frame(FLQuants(fit=catch.n(x), obs=as(catch.n(y), "FLQuant")))
	args$data$qname <- factor(args$data$qname, levels=c("obs", "fit"))
	args$x <- data~age|factor(year)
	args$type=c("l")
	args$groups <- quote(qname)
	args$ylab="numbers"
	args$xlab="age"
	args$scales=list(y=list(relation="free", draw=FALSE))
	args$auto.key <- list(points=FALSE, lines=TRUE, columns=2)
	args$par.settings=list(
		superpose.line=list(col=c("gray70", "black"), lty=1, lwd=c(2,1)), 
		strip.background=list(col="gray90"), 
		strip.border=list(col="black"))
	args$main="fitted and observed catch-at-age"
	do.call("xyplot", args)
})

##' @title testing
#' @name plot for fitted indices-at-age
#' @docType methods
#' @rdname ploti
#' @aliases plot,a4aFit,FLIndices-method
#' @description Method to plot fitted versus observed indices-at-age. Note the yaxis doesn't has a scale. The visual is about the difference between the two lines, not about the value of each line, which in any case would be very difficult to assess visually.
#' @param x an \code{a4aFit} object with the fitted values
#' @param y an \code{FLIndices} object with the observed values
#' @param ... additional argument list that might never be used
#' @return a \code{plot} with fitted and observed indices-at-age
#' @examples
#' data(ple4)
#' data(ple4.index)
#' obj <- sca(ple4, FLIndices(ple4.index))
#' plot(obj, FLIndices(ple4.index))

setMethod("plot", c("a4aFit", "FLIndices"), function(x, y, ...){
	v <- unlist(lapply(y, is, "FLIndexBiomass"))
	if(sum(v)>0){
		warning("Biomass indices will be removed, can't plot age based estimates.")
		y <- y[!v]
	}	
	args <- list()
	dfx <- as.data.frame(index(x))
	dfy <- as.data.frame(lapply(y, index))
	dfx$src="fit"
	dfy$src="obs"
	df0 <- rbind(dfx, dfy)
	df0$src <- factor(df0$src, levels=c("obs", "fit"))
	args$x <- data~age|factor(year)
	args$type=c("l")
	args$groups <- quote(src)
	args$ylab="numbers"
	args$xlab=""
	args$scales=list(y=list(relation="free", draw=FALSE))
	args$auto.key=list(lines=TRUE, points=FALSE, columns=2)
	args$par.settings=list(
		superpose.line=list(col=c("gray70", "black"), lty=1, lwd=c(2,1)), 
		strip.background=list(col="gray90"), 
		strip.border=list(col="black"))
	args$main="fitted and observed index-at-age"
	if(length(index(x))>1){
		for(i in names(y)){
			args$data <- subset(df0, qname==i)
			args$layout <- c(0,length(unique(args$data$year)))
			args$main <- paste(i, " fitted and observed index-at-age", sep=":")
			print(do.call("xyplot", args))
		}
	} else {
		args$data <- df0 
		args$layout <- c(0,length(unique(args$data$year)))
		do.call("xyplot", args)
	}
})


##' @title wireframe plot for FLQuant
#' @name wireframe plot for FLQuant
#' @docType methods
#' @rdname wireframe
#' @aliases wireframe,FLQuant-method wireframe,FLQuant,missing-method
#' @description Method to 3D plot \code{FLQuant} objects.
#' @param x a \code{FLQuant} 
#' @param y missing
#' @param screen list with numeric components 'x','y' and 'z' to change the 3D perspective 
#' @param ... additional argument list for the lattice engine 
#' @return a 3D surface plot
#' @examples
#' data(ple4)
#' wireframe(harvest(ple4))

setMethod("wireframe", c("FLQuant", "missing"), function(x, y, screen=list(x = -90, y=-45), ...){
	args <- list()
	args$data <- as.data.frame(x)
	args$x <- data~age*year
	args$screen = screen
	args$par.settings = list(axis.line = list(col = 'transparent'), box.3d=list(col="transparent")) 
	args$scales=list(col=1)
	do.call("wireframe", args)
})

#====================================================================
# get tpl
#====================================================================

#' @title Get TPL with ADMB code 
#' @name getTPL
#' @docType methods
#' @rdname getTPL
##' @aliases getTPL
#' @description Function to get the a4a TPL file with ADMB code and copy into a specific folder.
#' @param dir folder where the a4a.tpl file will be copied to.
#' @return file a4a.tpl
#' @examples
#' getTPL("myfolder")

getTPL <- function(dir){
	to <- paste0(dir,"/a4a.tpl", sep="")
	from <- system.file("admb/a4a.tpl", package = "FLa4a")
	file.copy(from, to)
}

