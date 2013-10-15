# FLStocka4a-class - «Short one line description»
# FLStocka4a-class

#' FLStocka4a extends \code{"FLStock"} class.
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
#' @name FLStocka4a-class
#' @rdname FLStocka4a-class
#' @exportClass FLStocka4a
setClass("FLStocka4a",
        representation(
                "FLStock",
                ssb = "FLQuantExt",
				fbar = "FLQuantExt"),
        prototype = prototype(
                ssb = new("FLQuantExt"),
				fbar = new("FLQuantExt"))
)

# constructor

setGeneric("FLStocka4a", function(object, ...)
	standardGeneric("FLStocka4a"))

setMethod("FLStocka4a", signature(object="missing"),
  function(...) {
    # empty
  	if(missing(...)){
	  	new("FLStocka4a")
    # or not
  	} else {
      args <- list(...)
	  args$Class <- 'FLStocka4a'
      do.call("new", args)
	  }
  }
)

# accessors

#setMethod("call", "a4aFitSA", function(object) object@call)

setMethod("ssb", "FLIndexa4a", function(object) object@ssb)
setMethod("fbar", "FLIndexa4a", function(object) object@fbar)

