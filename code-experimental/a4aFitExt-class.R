#' @title S4 class \code{a4aFitExt}
#'
#' @description The \code{a4aFitExt} calss extends \code{a4aFitSA} and was built to store extra slots like the hessian, which maybe needed for advanced users.
#'
#' @section Slots:
#' \describe{
#'    \item{Sigma}{The variance-covariance matrix.}
#'
#'    \item{L}{Likelihood}
#'
#'    \item{designMatrix}{The design matrix of the model}
#'
#'    \item{baseLvlPars}{...}
#'
#'  }
#' @template Accessors
#' @template Constructors
#' @docType class
#' @name a4aFitExt-class
#' @rdname a4aFitExt-class
#' @alias a4aFitExt-class

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

#' @rdname a4aFitExt-class
#' @alias a4aFitExt a4aFitExt-methods
setGeneric("a4aFitExt", function(object, ...)	standardGeneric("a4aFitExt"))

#' @rdname a4aFitExt-class
#' @alias a4aFitExt,missing-method
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

#' @rdname a4aFitExt-class
#' @alias a4aFitExt,a4aFit-method
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

#' @rdname a4aFitExt-class
#' @alias a4aFitExt,a4aFitSA-method
setMethod("a4aFitExt", signature(object="a4aFitSA"),
  function(object, ...) {
    out <- a4aFitExt(a4aFit(object))
    out @ pars <- object @ pars
    out
  }
)

#' @rdname a4aFitExt-class
#' 
#' @alias designMatrix,a4aFitExt-method
setGeneric("designMatrix", function(object, ...) standardGeneric("designMatrix"))
setMethod("designMatrix", "a4aFitExt", function(object) object@designMatrix)

#' @rdname a4aFitExt-class
#' @alias baseLvlPars,a4aFitExt-method
setGeneric("baseLvlPars", function(object, ...) standardGeneric("baseLvlPars"))
setMethod("baseLvlPars", "a4aFitExt", function(object) object@baseLvlPars)

#' @rdname a4aFitExt-class
#' @alias L,a4aFitExt-method
setGeneric("L", function(object, ...) standardGeneric("L"))
setMethod("L", "a4aFitExt", function(object) object@L)

#' @rdname a4aFitExt-class
#' @alias Sigma,a4aFitExt-method
setGeneric("Sigma", function(object, ...) standardGeneric("Sigma"))
setMethod("Sigma", "a4aFitExt", function(object) object@Sigma)

