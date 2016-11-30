
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
#'	\item{\code{grMod}}{the formula for the growth model, \emph{e.g.} von Bertallanffy}
#'	\item{\code{grInvMod}}{the formula for the inverse of the growth model, having length as the independent variable}
#'	\item{\code{params}}{an FLPar object with the parameters of the model; must match the equations in the models}
#'	\item{\code{vcov}}{an array with the variance covariance matrix of the parameters}
#'	\item{\code{distr}}{a character with the parameters' statistical distribution; it must match a known distribution for R (\emph{e.g.} "norm" for gaussian), so that \code{rnorm} can be called}
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
#' @template bothargs
#' @aliases a4aGr a4aGr-methods
#' @template Accessors
#' @template Constructors
#' @examples
#' mm <- matrix(NA, ncol=3, nrow=3)
#' diag(mm) <- c(50, 0.001,0.001)
#' mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
#' vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","yr^-1","yr")), vcov=mm, distr="norm")
setGeneric("a4aGr", function(object, ...) standardGeneric("a4aGr"))
#' @rdname a4aGr-class
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
#' @aliases grMod grMod-methods
setGeneric("grMod", function(object, ...) standardGeneric("grMod"))
#' @rdname a4aGr-class
setMethod("grMod", "a4aGr", function(object) object@grMod)

#' @rdname a4aGr-class
#' @param value the new object
#' @aliases grMod<- grMod<--methods
setGeneric("grMod<-", function(object,value) standardGeneric("grMod<-"))
#' @rdname a4aGr-class
setReplaceMethod("grMod", signature("a4aGr","formula"), function(object, value){
	object@grMod <- value
	object
})

#' @rdname a4aGr-class
#' @aliases grInvMod grInvMod-methods
setGeneric("grInvMod", function(object, ...) standardGeneric("grInvMod"))
#' @rdname a4aGr-class
setMethod("grInvMod", "a4aGr", function(object) object@grInvMod)

#' @rdname a4aGr-class
#' @aliases grInvMod<- grInvMod<--methods
setGeneric("grInvMod<-", function(object,value) standardGeneric("grInvMod<-"))
#' @rdname a4aGr-class
setReplaceMethod("grInvMod", signature("a4aGr","formula"), function(object, value){
	object@grInvMod <- value
	object
})

#' @rdname a4aGr-class
setMethod("params", "a4aGr", function(object) object@params)

#' @rdname a4aGr-class
setReplaceMethod("params", signature("a4aGr","FLPar"), function(object, value){
	object@params <- value
	object
})

#' @rdname a4aGr-class
setMethod("distr", "a4aGr", function(object) object@distr)

#' @rdname a4aGr-class
setReplaceMethod("distr", signature("a4aGr","character"), function(object, value){
	object@distr <- value
	object
})

#' @rdname a4aGr-class
setMethod("vcov", "a4aGr", function(object) object@vcov)

#' @rdname a4aGr-class
setReplaceMethod("vcov", signature(object = "a4aGr", value = "numeric"), function(object, value){
	object@vcov <- value
	object
})
