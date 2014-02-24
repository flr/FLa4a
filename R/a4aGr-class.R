
# validity: right side of the inverse model must have a "len" variable
valida4aGr <- function(object){

	pars <- object@params	
	v <- c(pars2dim(pars), grep("len", deparse(object@grInvMod)))

	if(sum(!v)>0) return("Object is not valid. Check that params have 2 dimensions with the second being named \"iter\" and that the inverse growth model has a term called \"len\".")
	return(TRUE)
}

#' @title Individual growth class
#' @name a4aGr
#' @rdname a4aGr-class
#' @template ClassDescription
#' @section Slot: 
#' \describe{
#'	\item{\code{grMod}}{the formula of the growth model, \emph{e.g.} von bertallanffy}
#'	\item{\code{grInvMod}}{the formula of the inverse of the growth model having length as the independent variable}
#'	\item{\code{params}}{a FLPar object with the parameters of the model. Must match the equations in the models}
#'	\item{\code{vcov}}{an array with the variance covariance matrix of the parameters}
#'	\item{\code{distr}}{a character with the parameters statistical distribution, it must match a known distribution for R, \emph{e.g.} "norm" for gaussian, so that \code{rnorm} can be called}
#' }
#' @alias a4aGr-class

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
#' @alias a4aGr a4aGr-methods a4aGr,missing-method
#' @section Constructor: 
#' \describe{
#'	\item{grMod}{a \code{formula} with the growth model (length~age)}
#' 	\item{grInvMod}{a \code{formula} with the inverse growth model (age~length)}
#' 	\item{params}{a \code{FLPar} object with the parameters of the models}
#' 	\item{vcov}{a \code{array} with the variance covariance matrix of the parameters}
#' 	\item{distr}{a \code{character} with the distribution of the parameters}
#' }
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
#' @alias grMod grMod-methods grMod,a4aGr-method
setGeneric("grMod", function(object, ...) standardGeneric("grMod"))
setMethod("grMod", "a4aGr", function(object) object@grMod)

#' @rdname a4aGr-class
#' @alias grMod<- grMod<--methods grMod<-,a4aGr-method
setGeneric("grMod<-", function(object,value) standardGeneric("grMod<-"))
setReplaceMethod("grMod", "a4aGr", function(object, value){
	if(all.equal(is(value), is(object@grMod))) object@grMod <- value
	object
})

#' @rdname a4aGr-class
#' @alias grInvMod grInvMod-methods grInvMod,a4aGr-method
setGeneric("grInvMod", function(object, ...) standardGeneric("grInvMod"))
setMethod("grInvMod", "a4aGr", function(object) object@grInvMod)

#' @rdname a4aGr-class
#' @alias grInvMod<- grInvMod<--methods grInvMod<-,a4aGr-method
setGeneric("grInvMod<-", function(object,value) standardGeneric("grInvMod<-"))
setReplaceMethod("grInvMod", "a4aGr", function(object, value){
	if(all.equal(is(value), is(object@grInvMod))) object@grInvMod <- value
	object
})

#' @rdname a4aGr-class
#' @alias params,a4aGr-method
setMethod("params", "a4aGr", function(object) object@params)

#' @rdname a4aGr-class
#' @alias params<- params<-methods params<-,a4aGr-method
setGeneric("params<-", function(object, value) standardGeneric("params<-"))
setReplaceMethod("params", "a4aGr", function(object, value){
	if(all.equal(is(value), is(object@params))) object@params <- value
	object
})

#' @rdname a4aGr-class
#' @alias distr distr-methods distr,a4aGr-method
setGeneric("distr", function(object, ...) standardGeneric("distr"))
setMethod("distr", "a4aGr", function(object) object@distr)

#' @rdname a4aGr-class
#' @alias distr<- distr<--methods distr<-,a4aGr-method
setGeneric("distr<-", function(object, value) standardGeneric("distr<-"))
setReplaceMethod("distr", "a4aGr", function(object, value){
	if(all.equal(is(value), is(object@distr))) object@distr <- value
	object
})

#' @rdname a4aGr-class
#' @alias vcov,a4aGr-method
setMethod("vcov", "a4aGr", function(object) object@vcov)

#' @rdname a4aGr-class
#' @alias vcov<-,a4aGr-method
setReplaceMethod("vcov", "a4aGr", function(object, value){
	if(all.equal(is(value), is(object@vcov))) object@vcov <- value
	object
})

#' @name rnglen
#' @rdname rnglen
#' @title length range
#' @description method to extract from \code{a4aGr} objects the length range.
#' @param object a \code{a4aGr} object
#' @return a \code{vector} object
#' @alias rnglen rnglen-methods rnglen,a4aGr-method
setGeneric("rnglen", function(object, ...) standardGeneric("rnglen"))
setMethod("rnglen", "a4aGr", function(object){
	object@range[c("min","max")]
})

#' @name rnglen
#' @rdname rnglenrplc
#' @title length range replacement
#' @description method to replace \code{a4aGr} object's length range.
#' @param object a \code{a4aGr} object
#' @param value a \code{vector} with max and min age range 
#' @return a \code{a4aGr} object
#' @alias rnglen<- rnglen<--methods rnglen<-,a4aGr-method
setGeneric("rnglen<-", function(object,value) standardGeneric("rnglen<-"))
setReplaceMethod("rnglen", "a4aGr", function(object, value){
	object@range[c("min","max")] <- sort(value)
	object
})

#' @name veclen
#' @rdname veclen
#' @title length vector
#' @description method to extract from \code{a4aGr} objects the vector of lengths.
#' @param object a \code{a4aGr} object
#' @return a \code{vector} object
#' @alias veclen veclen-methods veclen,a4aGr-method
setGeneric("veclen", function(object, ...) standardGeneric("veclen"))
setMethod("veclen", "a4aGr", function(object){
	rng <- object@range[c("min","max")]
	rng[1]:rng[2]
})

#' @name niters
#' @rdname niters
#' @title number of iterations
#' @description method to extract from \code{a4aGr} objects the number of iterations.
#' @param object a \code{a4aGr} object
#' @return a \code{numeric} object
#' @alias niters,a4aGr-method
setMethod("niters", "a4aGr", function(object){
	dim(params(object))[2]
})


