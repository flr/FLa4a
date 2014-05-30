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

#' @title Retrospective analysis
#' @name ra
#' @rdname ra
#' @description Runs a retrospective analysis of the fit.
#' @param stock A \code{FLStock} object.
#' @param indices A \code{FLIndices} object with the abundance and biomass indices.
#' @param n A \code{numeric} with the number of years to go back in the retrospectiva analysis.
#' @param ... To pass the model specifications to fit the model.
#' @return A \code{FLStocks} object with the \code{n} stocks resulting from removing one year and refitting the model.
#' @note The default method will use the default submodels of \code{sca}.
#' @aliases ra ra-methods ra,a4aFit-method
#' @examples
#' # load some data
#' data(ple4)
#' data(ple4.indices)
#' retro <- ra(ple4, ple4.indices, 5)
#' plot(retro)

setGeneric("ra", function(stock, indices, ...) standardGeneric("ra"))
setMethod("ra", c("FLStock","FLIndices"), function(stock, indices, n, ...){

	args <- list(...)
	stkdms <- dims(stock)
	
	if(is.null(args$qmodel)){
	  inddms <- lapply(indices, dims)
	  ka <- sapply(inddms, "[[", "age")
	  ka <- pmin(pmax(2, floor(.7 * ka)), 7)
	  args$qmodel <- lapply(ka, function(i) if (i == 2) ~ age else formula(paste("~ s(age, k = ", ka, ")")))
	}

	if(is.null(args$n1model)){
		if (stkdms$age > 10) {
			args$n1model <- ~ s(age, k = 10) 
	  	} else {
			args$n1model <- ~ factor(age)
	  	}  
	}
	
	if(is.null(args$vmodel)){
		args$vmodel  <- lapply(seq(length(indices) + 1), function(i) ~ 1)
		args$vmodel[[1]] <- ~ s(age, k = 3)
	}

	mxyr <- range(stock)['maxyear']
	
	lst <- list()
	length(lst) <- n

	for(i in 1:n){
		# the fmodel must be dealt inside the loop because the
		# number of years change and so does ky
		if(is.null(args$fmodel)){
			# this must be in agreement with sca default
			ka <- stkdms$age
			ka <- if (ka < 3) ka else min(max(3, floor(.5 * ka)), 6)
			ky <- floor(.5 * (stkdms$year - i + 1))
			if (ka >= 3) {
				args$fmodel <- formula(paste("~ te(age, year, k = c(", ka,",", ky,"), bs = 'tp')"))
			} else {
				args$fmodel <- formula(paste("~ age + s(year, k = ", ky,")"))
			}
			args$stock <- window(stock, end=mxyr-i+1)
			args$indices <- window(indices, end=mxyr-i+1)
			lst[[i]] <- args$stock + do.call("a4aSCA", args)
			args$fmodel <- NULL
		} else {
			args$stock <- window(stock, end=mxyr-i+1)
			args$indices <- window(indices, end=mxyr-i+1)
			lst[[i]] <- args$stock + do.call("a4aSCA", args)
		}
	}
	FLStocks(lst)
})


#' @title plot of fitted catch numbers-at-age
#' @name plot
#' @docType methods
#' @rdname plotc
#' @aliases plot,a4aFit,FLStock-method
#' @description Method to plot fitted versus observed catch numbers-at-age
#' @param x a \code{a4aFit} object with the fitted values
#' @param y a \code{FLStock} object with the observed values
#' @param ... Additional argument list that might not ever be used.
#' @return a \code{plot} with fitted and observed catch number-at-age
#' @template runsca
#' @examples
#' plot(fit, ple4)

setMethod("plot", c("a4aFit", "FLStock"), function(x, y, ...){
	args <- list()
	args$data <- as.data.frame(FLQuants(fit=catch.n(x), obs=catch.n(y)))
	args$x <- data~age|factor(year)
	args$type=c("l")
	args$groups <- quote(qname)
	args$ylab="numbers"
	args$xlab=""
	args$scales=list(y="free")
#	args$panel=function(x,y,...){
#		panel.abline(h=0, col.line="gray80")
#		panel.xyplot(x,y,...)
#		}
#	args$par.settings=list(
#		superpose.symbol=list(col="gray50", pch=19, cex=0.2), 
#		superpose.line=list(col=1, lty=1, lwd=2), 
#		strip.background=list(col="gray90"), 
#		strip.border=list(col="gray90"), 
#		box.rectangle=list(col="gray90"))
#	args$main="log residuals of catch and abundance indices"
	do.call("xyplot", args)
})

#' @title plot of fitted indices-at-age
#' @name plot
#' @docType methods
#' @rdname ploti
#' @aliases plot,a4aFit,FLIndices-method
#' @description Method to plot fitted versus observed indices-at-age
#'
#' @param x a \code{a4aFit} object with the fitted values
#' @param y a \code{FLIndices} object with the observed values
#' @param ... Additional argument list that might not ever be used.
#' @return a \code{plot} with fitted and observed indices-at-age
#' @template runsca
#' @examples
#' plot(fit, ple4.indices)

setMethod("plot", c("a4aFit", "FLIndices"), function(x, y, ...){
	args <- list()
	dfx <- as.data.frame(index(x))
	dfy <- as.data.frame(lapply(y, index))
	dfx$src="fit"
	dfy$src="obs"
	df0 <- rbind(dfx, dfy)
	args$x <- data~age|factor(year)*qname
	args$type=c("l")
	args$groups <- quote(src)
	args$ylab="numbers"
	args$xlab=""
	args$scales=list(y="free")
#	args$panel=function(x,y,...){
#		panel.abline(h=0, col.line="gray80")
#		panel.xyplot(x,y,...)
#		}
#	args$par.settings=list(
#		superpose.symbol=list(col="gray50", pch=19, cex=0.2), 
#		superpose.line=list(col=1, lty=1, lwd=2), 
#		strip.background=list(col="gray90"), 
#		strip.border=list(col="gray90"), 
#		box.rectangle=list(col="gray90"))
#	args$main="log residuals of catch and abundance indices"

	if(length(index(x))>1){
		for(i in names(y)){
			x11()
			args$data <- subset(df0, qname==i)
			args$layout <- c(0,length(unique(args$data$year)))
			print(do.call("xyplot", args))
		}
	} else {
		args$data <- df0 
		do.call("xyplot", args)
	}
})





