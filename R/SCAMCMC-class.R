#' @title MCMC settings class
#' @docType class
#' @name SCAMCMC
#' @rdname SCAMCMC-class
#' @template ClassDescription
#' @section Slot: 
#' \describe{
#'
#'	\item{\code{mcmc N}}{Run N MCMC iterations}
#'	\item{\code{mcsave N}}{Save every N th MCMC iteration}
#'	\item{\code{mcscale N}}{Rescale step size for first N iterations}
#'	\item{\code{mcmult N}}{Rescale the covariance matrix}
#'	\item{\code{mcrb N}}{Reduce high parameter correlations}
#'	\item{\code{mcprobe X}}{Use a fat-tailed proposal distribution}
#'	\item{\code{mcdiag}}{Use a diagonal covariance matrix}
#'	\item{\code{mcnoscale}}{Do not scale the algorithm during}
#'	\item{\code{mcu}}{Use a uniform distribution as proposal distribution}
#'	\item{\code{hybrid}}{Use the hybrid method}
#'	\item{\code{hynstep N}}{Mean number of steps for the leapfrog method}
#'	\item{\code{hyeps X}}{The stepsize for the leapfrog method [X numeric and > 0]}
#'
#' }
#' @aliases SCAMCMC-class
setClass("SCAMCMC",
        representation(
			mcmc	= "numeric",
			mcsave	= "numeric",
			mcscale	= "numeric",
			mcmult	= "numeric",
			mcrb	= "numeric",
			mcprobe	= "numeric",
			mcseed	= "numeric",
			mcdiag	= "logical",
			mcnoscale	= "logical",
			mcu	= "logical",
			hybrid	= "logical",
			hynstep = "numeric",
			hyeps	= "numeric"
		),
        prototype = prototype(
			mcmc	= 10000,
			mcsave	= 100,
			mcscale	= NaN,
			mcmult	= NaN,
			mcrb	= NaN,
			mcprobe	= NaN,
			mcseed	= NaN,
			mcdiag	= FALSE,
			mcnoscale	= FALSE,
			mcu	= FALSE,
			hybrid	= FALSE,
			hynstep = NaN,
			hyeps	= NaN
        ),
        validity = function(object) {
			# if hybrid mcsave must be 1
			if(object@hybrid & object@mcsave!=1)
				return("to use the hybrid method mcsave must be 1")
			# Everything is fine
			return(TRUE)}
)

#' @rdname SCAMCMC-class
#' @template Accessors
#' @template Constructors
#' @aliases SCAMCMC SCAMCMC-methods SCAMCMC,missing-method
setGeneric("SCAMCMC", function(object, ...) standardGeneric("SCAMCMC"))
setMethod("SCAMCMC", signature(object="missing"),
  function(...) {
    # empty
  	if(missing(...)){
	  	new("SCAMCMC")
    # or not
  	} else {
      args <- list(...)
	  args$Class <- 'SCAMCMC'
      do.call("new", args)
	  }
  }
)

# accessors

##' @rdname SCAMCMC-class
##' @aliases mcmc mcmc-methods mcmc,SCAMCMC-method
#setGeneric("mcmc", function(object, ...) standardGeneric("mcmc"))
#setMethod("mcmc", "SCAMCMC", function(object) object @ mcmc)

#' @rdname SCAMCMC-class
#' @aliases SCAMCMC SCAMCMC-methods SCAMCMC,missing-method

setGeneric("getADMBCallArgs", function(object, ...) standardGeneric("getADMBCallArgs"))
setMethod("getADMBCallArgs", signature(object="SCAMCMC"),
  function(object, ...) {
	slts <- getSlots("SCAMCMC")
	lslts <- slts[slts=="logical"]
	nslts <- slts[slts!="logical"]
	callargs <- c("")
	for(i in names(lslts)){
		if(isTRUE(slot(object, i))) callargs <- paste(callargs, " -", i, sep="")
	}

	for(i in names(nslts)){
		if(!is.na(slot(object, i))) callargs <- paste(callargs, " -", i, " ", slot(object, i), sep="")
	}
	callargs
  }
)

#' @rdname SCAMCMC-class
#' @aliases getADMBCallArgs getADMBCallArgs-methods SCAMCMC,missing-method

setGeneric("getN", function(object, ...) standardGeneric("getN"))
setMethod("getN", signature(object="SCAMCMC"),
  function(object, ...) {
	floor(object@mcmc/object@mcsave)
  }
)





















