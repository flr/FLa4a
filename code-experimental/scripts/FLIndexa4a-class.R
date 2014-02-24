# FLIndexa4a-class - «Short one line description»
# FLIndexa4a-class

#' FLIndexa4a extends \code{"FLIndex"} class.
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
#' @name FLIndexa4a-class
#' @rdname FLIndexa4a-class
#' @exportClass FLIndexa4a
setClass("FLIndexa4a",
        representation(
                "FLIndex",
                index.covar = "FLQuant",
        prototype = prototype(
                index.covar = new("FLQuant"))
)

# constructor

setGeneric("FLIndexa4a", function(object, ...) standardGeneric("FLIndexa4a"))

setMethod("FLIndexa4a", signature(object="missing"),
  function(...) {
    # empty
  	if(missing(...)){
	  	new("FLIndexa4a")
    # or not
  	} else {
      args <- list(...)
	  args$Class <- 'FLIndexa4a'
      do.call("new", args)
	  }
  }
)

# accessors

#setMethod("call", "a4aFitSA", function(object) object@call)

setGeneric("index.covar", function(object, ...) standardGeneric("index.covar"))
setMethod("index.covar", "FLIndexa4a", function(object) object@index.covar)

