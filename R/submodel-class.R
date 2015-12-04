#' @title Submodel class
#' @docType class
#' @name submodel
#' @rdname submodel-class
#' @template ClassDescription
#' @section Slot: 
#' \describe{
#'
#'	\item{\code{Mod}}{\code{formula} describing the model}
#'
#'	\item{\code{params}}{\code{FLPar} with model parameters}
#'
#'	\item{\code{vcov}}{\code{array} with variance covariance paramaters related to the variance model}
#'
#'	\item{\code{centering}}{\code{numeric} value used for centering the data}
#'
#'	\item{\code{distr}}{a character with the parameters' statistical distribution; it must match a known distribution for R (\emph{e.g.} "norm" for gaussian) so that \code{rnorm} can be called}
#' }
#' @aliases submodel-class
setClass("submodel",
        representation(
				"FLComp",
        Mod       = "formula",
				params    = "FLPar",
				vcov      = "array",
				centering = "numeric",
				distr     = "character"),
        prototype = prototype(
				name	= character(0),
				desc	= character(0),
				range	= c(min=0, max=0, plusgroup=0, minyear=0, maxyear=0),
        Mod = ~1,
				params = FLPar(),
				vcov = array(),
				centering = 0,
				distr = "lnorm",
		validity=function(object) {
			# no unit, area, season fits
			if(length(dim(object@params)) > 2 | length(dim(object@vcov)) > 3){
				return("Params or vcov have unit, area or season. Can't work with that!")
			}
			# Everything is fine
			return(TRUE)})
)


#' @rdname submodel-class
#' @template Accessors
#' @template Constructors
#' @aliases submodel submodel-methods submodel,missing-method
setGeneric("submodel", function(object, ...)
	standardGeneric("submodel"))

setMethod("submodel", signature(object="missing"),
  function(...) {
    # empty
  	if(missing(...)){
	  	new("submodel")
    # or not
  	} else {
      args <- list(...)
	  args$Class <- 'submodel'
      do.call("new", args)
	  }
  }
)

#' @rdname submodel-class
#' @aliases params,submodel-method
setMethod("params", "submodel", function(object) object@params)

#' @rdname submodel-class
#' @aliases model,submodel-method
setMethod("model", "submodel", function(object) object@model)

#' @rdname submodel-class
#' @aliases covar,submodel-method
setMethod("covar", "submodel", function(object) object@covar)




