
# validity: right side of the inverse model must have a "len" variable
valida4aGr <- function(object){

	pars <- object@params	
	v <- c(pars2dim(pars), grep("len", deparse(object@grInvMod)))

	if(sum(!v)>0) return("Object is not valid. Check that params have 2 dimensions with the second being named \"iter\" and that the inverse growth model has a term called \"len\".")
	return(TRUE)
}

#' @title Individual growth class
#' @docType class
#' @name a4aGr
#' @rdname a4aGr-class
#' @template ClassDescription
#' @section Slot: 
#' \describe{
#'
#'	\item{\code{grMod}}{the formula of the growth model, \emph{e.g.} von bertallanffy}
#'
#'	\item{\code{grInvMod}}{the formula of the inverse of the growth model having length as the independent variable}
#'
#'	\item{\code{params}}{a FLPar object with the parameters of the model. Must match the equations in the models}
#'
#'	\item{\code{vcov}}{an array with the variance covariance matrix of the parameters}
#'
#'	\item{\code{distr}}{a character with the parameters statistical distribution, it must match a known distribution for R, \emph{e.g.} "norm" for gaussian, so that \code{rnorm} can be called}
#'
#' }
#' @aliases a4aGr-class

setClass("a4aGr",
        representation(
				"FLComp",
                grMod = "formula",
                grInvMod = "formula",
				params="FLPar",
				vcov="array",
				distr="character"),
        prototype = prototype(
				name	= character(0),
				desc	= character(0),
				range	= c(min=0, max=0, plusgroup=0, minyear=0, maxyear=0),
                grMod = ~1,
                grInvMod = ~1,
				params=FLPar(),
				vcov=array(),
				distr="norm"),
				validity=valida4aGr
)

#' @rdname a4aGr-class
#' @aliases a4aGr a4aGr-methods a4aGr,missing-method
#' @template Accessors
#' @template Constructors
#' @examples
#' mm <- matrix(NA, ncol=3, nrow=3)
#' diag(mm) <- c(50, 0.001,0.001)
#' mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
#' vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","ano-1","ano")), vcov=mm, distr="norm")

setGeneric("a4aGr", function(object, ...) standardGeneric("a4aGr"))
setMethod("a4aGr", signature(object="missing"),
  function(...) {
    # empty
  	if(missing(...)){
	  	new("a4aGr")
    # or not
  	} else {
      args <- list(...)
	  args$Class <- 'a4aGr'
      do.call("new", args)
	  }
  }
)

#' @rdname a4aGr-class
#' @aliases grMod grMod-methods grMod,a4aGr-method
setGeneric("grMod", function(object, ...) standardGeneric("grMod"))
setMethod("grMod", "a4aGr", function(object) object@grMod)

#' @rdname a4aGr-class
#' @aliases grMod<- grMod<--methods grMod<-,a4aGr,formula-method
setGeneric("grMod<-", function(object,value) standardGeneric("grMod<-"))
setReplaceMethod("grMod", signature("a4aGr","formula"), function(object, value){
	object@grMod <- value
	object
})

#' @rdname a4aGr-class
#' @aliases grInvMod grInvMod-methods grInvMod,a4aGr-method
setGeneric("grInvMod", function(object, ...) standardGeneric("grInvMod"))
setMethod("grInvMod", "a4aGr", function(object) object@grInvMod)

#' @rdname a4aGr-class
#' @aliases grInvMod<- grInvMod<--methods grInvMod<-,a4aGr,formula-method
setGeneric("grInvMod<-", function(object,value) standardGeneric("grInvMod<-"))
setReplaceMethod("grInvMod", signature("a4aGr","formula"), function(object, value){
	object@grInvMod <- value
	object
})

#' @rdname a4aGr-class
#' @aliases params,a4aGr-method
setMethod("params", "a4aGr", function(object) object@params)

#' @name a4aGr
#' @rdname a4aGr-class
#' @aliases params<-,a4aGr,FLPar-method
setReplaceMethod("params", signature("a4aGr","FLPar"), function(object, value){
	object@params <- value
	object
})

#' @rdname a4aGr-class
#' @aliases distr,a4aGr-method
setMethod("distr", "a4aGr", function(object) object@distr)

#' @name distr<- for a4aGr
#' @rdname a4aGr-class
#' @aliases distr<-,a4aGr,character-method
setReplaceMethod("distr", signature("a4aGr","character"), function(object, value){
	object@distr <- value
	object
})

#' @rdname a4aGr-class
#' @aliases vcov,a4aGr-method
setMethod("vcov", "a4aGr", function(object) object@vcov)

#' @name vcov<- for a4aGr
#' @rdname a4aGr-class
#' @aliases vcov<-,a4aGr,numeric-method
setReplaceMethod("vcov", signature(object = "a4aGr", value = "numeric"), function(object, value){
	object@vcov <- value
	object
})

#' @name niters for a4aGr
#' @rdname niters
#' @title number of iterations
#' @description method to extract from \code{a4aGr} objects the number of iterations.
#' @param object a \code{a4aGr} object
#' @return a \code{numeric} object
#' @aliases niters,a4aGr-method
setMethod("niters", "a4aGr", function(object){
	dim(params(object))[2]
})


