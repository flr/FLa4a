# a4aFit-class - «Short one line description»
# a4aFit-class

##
va4aFit <- function(object) {

        # All FLQuant objects must have same dimensions
        Dim <- dim(object@stock.n)
        if (dim(object@harvest) != Dim | dim(object@catch.n != Dim))
                return("stock.n, catch.n and harvest slots must have same dimensions")
        # Everything is fine
        return(TRUE)
}

#' a4aFit extends \code{"FLComp"} class.
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
#' @name a4aFit-class
#' @rdname a4aFit-class
#' @exportClass a4aFit
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
        validity = va4aFit
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
#' @rdname a4aFit-methods
#'
#' @examples
#' data(ple4)
setGeneric("a4aFit", function(object, ...) standardGeneric("a4aFit"))

#' Title 
#' @name a4aFit
#' @docType methods
#' @rdname a4aFit-methods
#' @aliases a4aFit,a4aFit-method
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



# show method

#' Hello 
#' @name show
#' @docType methods
#' @rdname show-methods
#' @aliases show,a4aFit-method
setMethod("show", signature(object = "a4aFit"),
  function(object) 
  {
    cat("a4a model fit for:", object@name, "\n")
    cat("\nCall:\n", paste(deparse(object @ call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")

    cat("Time used:\n")
    print(object @ clock)
  
    #frmtQ <- function(x) {
    #  x <- as.matrix(drop(x @ .Data))
    #  if (ncol(x) > 8) {
    #    x[] <- paste(round(x,2))
    #    x <- x[,c(1:4,-2:0 + ncol(x))]
    #    x[,4] <- "..."
    #    colnames(x)[4] <- "..."
    #  }
    #  print(x, quote = FALSE)
    #}
    cat("\nPredicted F-at-age:\n")
      show(harvest(object))
    cat("\nPredicted stock numbers-at-age:\n")
      show(stock.n(object))
 })

#setMethod("show", signature(object = "a4aFit"),
#  function(object) 
#  {

#    cat("Formulas:\n")
#    cat("   log F", format(object@models$fmodel), "\n")
#    sapply(object@models$qmodel, function(m) cat("   log Q", format(m), "\n"))
#    cat("   log R", format(object@models$rmodel), "\n")

#    print(x$formula)
#    cat("Total model degrees of freedom", object@fit.sum["nopar"], "\n")
#    cat("   Residual degrees of freedom", object@fit.sum["nobs"] - object@fit.sum["nopar"], "\n")
#    cat("\nEstimated degrees of freedom:\n")
    #TODO check that nlogl is not 2 x neg loglik
#    cat("\n            AIC: ", AIC(object), sep = "")
#    cat("\n            BIC: ", BIC(object), sep = "")
#    cat("\n log Likelihood: ", c(logLik(object)), "\n", sep = "")
# })


# accessors

#' Title 
#' @name stock.n
#' @docType methods
#' @rdname stock.n-methods
#' @aliases stock.n,a4aFit-method
setMethod("stock.n", "a4aFit", function(object) object@stock.n)

#' Title 
#' @name harvest
#' @docType methods
#' @rdname harvest-methods
#' @aliases harvest,FLa4aFit-method
setMethod("harvest", "a4aFit", function(object) object@harvest)

#' Title 
#' @name catch.n
#' @docType methods
#' @rdname catch.n-methods
#' @aliases catch.n,FLa4aFit-method
setMethod("catch.n", "a4aFit", function(object) object@catch.n)

#' Title 
#' @name catch.n
#' @docType methods
#' @rdname catch.n-methods
#' @aliases catch.n,FLa4aFit-method
setMethod("index", "a4aFit", function(object) object@index)

#' Method to extract the log likelihood 
#' @name logLik
#' @docType methods
#' @rdname logLik-methods
#' @aliases logLik,FLa4aFit-method
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

