#' @title S4 class \code{a4aFitMCMC}
#'
#' @description The \code{a4aFitMCMC} class extends \code{a4aFit} and was built to store the MCMC simulations of an a4a stock assessment fit.
#'
#' @section Slots:
#' \describe{
#'    \item{pars}{Matrix with model parameters}
#'
#'  }
#' @template Accessors
#' @template Constructors
#' @docType class
#' @name a4aFitMCMC-class
#' @rdname a4aFitMCMC-class
#' @alias a4aFitMCMC-class
setClass("a4aFitMCMC",
        representation(
                "a4aFit",
                pars    = "matrix"
                ),
        prototype = prototype(
                pars    = new('matrix'))
)

#' @rdname a4aFitMCMC-class
#' @alias a4aFitMCMC a4aFitMCMC-methods
setGeneric("a4aFitMCMC", function(object, ...) standardGeneric("a4aFitMCMC"))

#' @rdname a4aFitMCMC-class
#' @alias a4aFitMCMC,missing-method
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

#' @rdname a4aFitMCMC-class
#' @alias a4aFit,a4aFitMCMC-method
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

#' @rdname a4aFitMCMC-class
#' @alias a4aFitMCMC,a4aFit-method
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


