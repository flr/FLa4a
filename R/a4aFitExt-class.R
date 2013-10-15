#' a4aFitSA-class - «Short one line description»
#' a4aFitSA-class

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
#' @name a4aFitExt-class
#' @rdname a4aFitExt-class
#' @exportClass a4aFitExt
setClass("a4aFitExt",
        representation(
                "a4aFitSA",
                Sigma = "matrix",
                L = "matrix",
                designMatrix = "list",
                baseLvlPars = "numeric"
                ),
        prototype = prototype(
                Sigma = new("matrix"),
                L = new("matrix"),
                designMatrix = new("list"),
                baseLvlPars = new("numeric"))
)



# constructor

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("a4aFitExt", function(object, ...)	standardGeneric("a4aFitExt"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("a4aFitExt", signature(object="missing"),
  function(...) {
    # empty
  	if(missing(...)) {
	  	new("a4aFitExt")
    # or not
  	} else {
      args <- list(...)
	  args$Class <- 'a4aFitExt'
      do.call("new", args)
	  }
  }
)


setMethod("a4aFitExt", signature(object="a4aFit"),
  function(object, ...) {
    out <- a4aFitExt()
    out @ name    <- object @ name
    out @ desc    <- object @ desc
    out @ range   <- object @ range
    out @ call    <- object @ call
    out @ clock   <- object @ clock
    out @ stock.n <- object @ stock.n
    out @ harvest <- object @ harvest
    out @ catch.n <- object @ catch.n
    out @ index   <- object @ index
    out @ fitSumm   <- object @ fitSumm
    out
  }
)

setMethod("a4aFitExt", signature(object="a4aFitSA"),
  function(object, ...) {
    out <- a4aFitExt(a4aFit(object))
    out @ pars <- object @ pars
    out
  }
)


# accessors

#setMethod("call", "a4aFitSA", function(object) object@call)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("hessian", function(object, ...) standardGeneric("hessian"))
#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("hessian", "a4aFitExt", function(object) object@hessian)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("designMatrix", function(object, ...) standardGeneric("designMatrix"))
#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("designMatrix", "a4aFitExt", function(object) object@designMatrix)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("grad", function(object, ...) standardGeneric("grad"))
#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("grad", "a4aFitExt", function(object) object@grad)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("baseLvlPars", function(object, ...) standardGeneric("baseLvlPars"))
#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("baseLvlPars", "a4aFitExt", function(object) object@baseLvlPars)


