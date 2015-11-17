#' @title S4 class \code{a4aFit}
#'
#' @description The \code{a4aFit} class was built to store the a4a stock assessment fits.
#'
#' @section Slots:
#' \describe{
#'    \item{call}{The function call}
#'
#'    \item{clock}{Information on call duration}
#'
#'    \item{fitSumm}{Fit summary}
#'
#'    \item{stock.n}{Estimates of stock numbers-at-age}
#'
#'    \item{harvest}{Estimates of fishing mortality at age}
#'
#'    \item{catch.n}{Estimates of catch numbers-at-age}
#'
#'    \item{index}{Estimates of survey or CPUE indices-at-age}
#'
#'  }
#' @template Accessors
#' @template Constructors
#' @docType class
#' @name a4aFit-class
#' @rdname a4aFit-class
#' @aliases a4aFit-class
#' @examples
#' data(ple4)
#' data(ple4.index)
#'
#' obj <- sca(stock=ple4, indices=FLIndices(ple4.index))
#' obj
#'
#' slotNames(obj)
#' clock(obj)
#' fitSumm(obj)
#'
#' flq <- stock.n(obj)
#' is(flq)
#' flq <- index(obj)
#' is(flq)
#'
#' logLik(obj)
#' AIC(obj)
#' BIC(obj)

setClass("a4aFit",
        representation(
                "FLComp",
                call         = "call",
                clock        = "numeric",
                fitSumm      = "array",
                stock.n      = "FLQuant",
                harvest      = "FLQuant",
                catch.n      = "FLQuant",
                index        = "FLQuants"),
        prototype = prototype(
                call         = new('call'),
                clock        = new('numeric'),
                fitSumm      = new('array'),
                stock.n      = new('FLQuant'),
                harvest      = new('FLQuant'),
                catch.n      = new('FLQuant'),
                index        = new('FLQuants')),
        validity = function(object) {
			# All FLQuant objects must have same dimensions
			Dim <- dim(object@stock.n)
			if ((sum(dim(object@harvest) != Dim) + sum(dim(object@catch.n) != Dim))>=1)
				return("stock.n, catch.n and harvest slots must have same dimensions")
			# Everything is fine
			return(TRUE)}
)

#' @rdname a4aFit-class
#' @aliases a4aFit a4aFit-methods
setGeneric("a4aFit", function(object, ...) standardGeneric("a4aFit"))

#' @rdname a4aFit-class
#' @aliases a4aFit,missing-method
setMethod("a4aFit", signature(object="missing"),
  function(...) {
    # empty
  	if(missing(...)){
	  	new("a4aFit")
    # or not
  	} else {
      args <- list(...)
	  args$Class <- 'a4aFit'
      do.call("new", args)
	  }
  }
)

#' @rdname a4aFit-class
#' @aliases clock,a4aFit-method
setGeneric("clock", function(object, ...) standardGeneric("clock"))
setMethod("clock", "a4aFit", function(object) object@clock)

#' @rdname a4aFit-class
#' @aliases fitSumm fitSumm-methods fitSumm,a4aFit-method
setGeneric("fitSumm", function(object, ...) standardGeneric("fitSumm"))
setMethod("fitSumm", "a4aFit", function(object) object@fitSumm)

#' @rdname a4aFit-class
#' @aliases stock.n,a4aFit-method
setMethod("stock.n", "a4aFit", function(object) object@stock.n)

#' @rdname a4aFit-class
#' @aliases harvest,a4aFit-method
setMethod("harvest", "a4aFit", function(object) object@harvest)

#' @rdname a4aFit-class
#' @aliases catch.n,a4aFit-method
setMethod("catch.n", "a4aFit", function(object) object@catch.n)

#' @rdname a4aFit-class
#' @aliases index,a4aFit-method
setMethod("index", "a4aFit", function(object) object@index)

#' @rdname a4aFit-class
#' @aliases show,a4aFit-method
setMethod("show", signature(object = "a4aFit"),
  function(object) 
  {
    cat("a4a model fit for:", object@name, "\n")
    cat("\nCall:\n", paste(deparse(object @ call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")

    cat("Time used:\n")
    print(object @ clock)
  
 })

#' @rdname a4aFit-class
#' @aliases logLik,a4aFit-method
setMethod("logLik", signature(object = "a4aFit"),
  function(object, ...) 
  {  
    dim2 <- length(dim(object @ fitSumm))
    if (dim2 == 1) {
      val <- -1 * unname(object @ fitSumm["nlogl"])
      attr(val, "nobs") <- unname(object @ fitSumm["nobs"])
      attr(val, "df") <- unname(object@fitSumm["nopar"])
    } else if (dim2 == 2) {
      val <- -1 * unname(object @ fitSumm["nlogl",])
      attr(val, "nobs") <- unname(object @ fitSumm["nobs",])
      attr(val, "df") <- unname(object@fitSumm["nopar",])
    }
    class(val) <- "logLik"
    val
 })


