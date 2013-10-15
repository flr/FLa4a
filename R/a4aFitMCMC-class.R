


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
#' @name a4aFitSA-class
#' @rdname a4aFitSA-class
#' @exportClass a4aFitSA
setClass("a4aFitMCMC",
        representation(
                "a4aFit",
                pars    = "matrix"
                ),
        prototype = prototype(
                pars    = new('matrix'))
)

# show method





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
setGeneric("a4aFitMCMC", function(object, ...) standardGeneric("a4aFitMCMC"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("a4aFitMCMC", signature(object="missing"),
  function(...) {
    # empty
  	if(missing(...)){
	  	new("a4aFitMCMC")
    # or not
  	} else {
      args <- list(...)
	  args$Class <- 'a4aFitMCMC'
      do.call("new", args)
	  }
  }
)

setMethod("a4aFit", signature(object="a4aFitMCMC"),
  function(object, ...) {
    out <- a4aFit()
    out @ name    <- object @ name
    out @ desc    <- object @ desc
    out @ range   <- object @ range
    out @ call    <- object @ call
    out @ clock   <- object @ clock
    out @ stock.n <- object @ stock.n
    out @ harvest <- object @ harvest
    out @ catch.n <- object @ catch.n
    out @ index   <- object @ index
    out
  }
)


setMethod("a4aFitMCMC", signature(object="a4aFit"),
  function(object, ...) {
    out <- a4aFitMCMC()
    out @ name    <- object @ name
    out @ desc    <- object @ desc
    out @ range   <- object @ range
    out @ call    <- object @ call
    out @ clock   <- object @ clock
    out @ stock.n <- object @ stock.n
    out @ harvest <- object @ harvest
    out @ catch.n <- object @ catch.n
    out @ index   <- object @ index
    out
  }
)


