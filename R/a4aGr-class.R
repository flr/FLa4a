
# validity: right side of the inverse model must have a "len" variable
valida4aGr <- function(object){

	pars <- object@params	
	v <- c(pars2dim(pars), grep("len", deparse(object@grInvMod)))

	if(sum(!v)>0) return("Object is not valid. Check that params have 2 dimensions with the second being named \"iter\" and that the inverse growth model has a term called \"len\".")
	return(TRUE)
}

#' Individual growth class
#'
#' @section Slot: 
#' \itemize{
#'	\item \code{grMod} the formula of the growth model, \emph{e.g.} von bertallanffy.
#'	\item \code{grInvMod} the formula of the inverse of the growth model having length as the independent variable.
#'	\item \code{params} a FLPar object with the parameters of the model. Must match the equations in the models.
#'	\item \code{vcov} an array with the variance covariance matrix of the parameters.
#'	\item \code{vcov} a character with the parameters statistical distribution, it must match a known distribution for R, \emph{e.g.} "norm" for gaussian, so that \code{rnorm} can be called.
#' '}
#' 
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export

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

# constructor

#' @title Constructor for the \code{a4aGr} class 
#'
#' @description Constructor for the individual growth model class.
#'
#' @param grMod a \code{formula} with the growth model (length~age)
#' @param grInvMod a \code{formula} with the inverse growth model (age~length)
#' @param params a \code{FLPar} object with the parameters of the models
#' @param vcov a \code{array} with the variance covariance matrix of the parameters
#' @param distr a \code{character} with the distribution of the parameters
#' @return a \code{a4aGr} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
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

# accessors

#' @title growth model accessor
#' @description accessor method for \code{a4aGr} object's slot \code{grMod}.
#' @param object a \code{a4aGr} object
#' @return a \code{formula} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export

setGeneric("grMod", function(object, ...) standardGeneric("grMod"))
setMethod("grMod", "a4aGr", function(object) object@grMod)

#' @title set growth model
#' @description set method for \code{a4aGr} object's slot \code{grMod}.
#' @param object a \code{a4aGr} object
#' @param value a \code{formula} object
#' @return a \code{a4aGr} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("grMod<-", function(object,value) standardGeneric("grMod<-"))
setReplaceMethod("grMod", "a4aGr", function(object, value){
	if(all.equal(is(value), is(object@grMod))) object@grMod <- value
	object
})

#' @title growth inverse model accessor
#' @description accessor method for \code{a4aGr} object's slot \code{grInvMod}.
#' @param object a \code{a4aGr} object
#' @return a \code{formula} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("grInvMod", function(object, ...) standardGeneric("grInvMod"))
setMethod("grInvMod", "a4aGr", function(object) object@grInvMod)

#' @title set growth inverse model
#' @description set method for \code{a4aGr} object's slot \code{grInvMod}.
#' @param object a \code{a4aGr} object
#' @param value a \code{formula} object
#' @return a \code{a4aGr} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("grInvMod<-", function(object,value) standardGeneric("grInvMod<-"))
setReplaceMethod("grInvMod", "a4aGr", function(object, value){
	if(all.equal(is(value), is(object@grInvMod))) object@grInvMod <- value
	object
})

#' @title params accessor
#' @description accessor method for \code{a4aGr} object's slot \code{params}.
#' @param object a \code{a4aGr} object
#' @return a \code{FLPar} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setMethod("params", "a4aGr", function(object) object@params)

#' @title set parameters
#' @description set method for \code{a4aGr} object's slot \code{params}.
#' @param object a \code{a4aGr} object
#' @param value a \code{FLPar} object
#' @return a \code{a4aGr} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("params<-", function(object, value) standardGeneric("params<-"))
setReplaceMethod("params", "a4aGr", function(object, value){
	if(all.equal(is(value), is(object@params))) object@params <- value
	object
})

#' @title distribution accessor
#' @description accessor method for \code{a4aGr} object's slot \code{distr}.
#' @param object a \code{a4aGr} object
#' @return a \code{character} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("distr", function(object, ...) standardGeneric("distr"))
setMethod("distr", "a4aGr", function(object) object@distr)

#' @title set distribution
#' @description set method for \code{a4aGr} object's slot \code{distr}.
#' @param object a \code{a4aGr} object
#' @param value a \code{character} object
#' @return a \code{a4aGr} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("distr<-", function(object, value) standardGeneric("distr<-"))
setReplaceMethod("distr", "a4aGr", function(object, value){
	if(all.equal(is(value), is(object@distr))) object@distr <- value
	object
})


setMethod("vcov", "a4aGr", function(object) object@vcov)

setReplaceMethod("vcov", "a4aGr", function(object, value){
	if(all.equal(is(value), is(object@vcov))) object@vcov <- value
	object
})

#' @title length range
#' @description method to extract from \code{a4aGr} objects the length range.
#' @param object a \code{a4aGr} object
#' @return a \code{vector} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("rnglen", function(object, ...) standardGeneric("rnglen"))
setMethod("rnglen", "a4aGr", function(object){
	object@range[c("min","max")]
})

#' @title length range replacement
#' @description method to replace \code{a4aGr} object's length range.
#' @param object a \code{a4aGr} object
#' @param value a \code{vector} with max and min age range 
#' @return a \code{a4aGr} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("rnglen<-", function(object,value) standardGeneric("rnglen<-"))
setReplaceMethod("rnglen", "a4aGr", function(object, value){
	object@range[c("min","max")] <- sort(value)
	object
})

#' @title length vector
#' @description method to extract from \code{a4aGr} objects the vector of lengths.
#' @param object a \code{a4aGr} object
#' @return a \code{vector} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("veclen", function(object, ...) standardGeneric("veclen"))
setMethod("veclen", "a4aGr", function(object){
	rng <- object@range[c("min","max")]
	rng[1]:rng[2]
})

setMethod("niters", "a4aGr", function(object){
	dim(params(object))[2]
})


