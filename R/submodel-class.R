# submodel-class - «Short one line description»
# submodel-class

#' a4aFitSA extends \code{"a4aFit"} class.
#'
#' Some details about this class and my plans for it in the body.
#'
#' \describe{
#'    \item{myslot1}{A logical keeping track of something.}
#'
#'    \item{myslot2}{An integer specifying something else.}
#' 
#'    \item{myslot3}{A data.frame holding some data.}
#'  }
#' @name submodel-class
#' @rdname submodel-class
#' @exportClass submodel
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
				distr = "lnorm")
)



# constructor

#' @export
setGeneric("submodel", function(object, ...)
	standardGeneric("submodel"))

#' @export
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


# accessors

#' @export
setGeneric("pars", function(object, ...) standardGeneric("pars"))
#' @export
setMethod("pars", "submodel", function(object) object@pars)

#' @export
setMethod("model", "submodel", function(object) object@model)

#' @export
setMethod("covar", "submodel", function(object) object@covar)


