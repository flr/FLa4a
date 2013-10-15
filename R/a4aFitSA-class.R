


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
setClass("a4aFitSA",
        representation(
                "a4aFit",
                pars    = "SCAPars"
                ),
        prototype = prototype(
                pars    = new('SCAPars'))
)

# show method

#' Hello 
#' @name show
#' @docType methods
#' @rdname show-methods
#' @aliases show,a4aFitSA-method
setMethod("show", signature(object = "a4aFitSA"),
  function(object) 
  {

    show(a4aFit(object))
    
    cat("\nSubmodels:\n")  
    cat("\t fmodel: "); print( fmodel(pars(object)), showEnv = FALSE)
    cat("\tsrmodel: "); print(srmodel(pars(object)), showEnv = FALSE)
    cat("\tn1model: "); print(n1model(pars(object)), showEnv = FALSE)

    # something to format the qmodel and vmodel    
    printFormList <- function(frmL) {
      if (length(frmL) == 0) return(invisible(NA))
      mods <- lapply(frmL, model)
      mnames <- names(mods)
      maxname <- max(sapply(mnames, nchar))
      for (i in seq(length(mods))) {
        cat("\t   ", mnames[i], ": ", rep(" ", maxname - nchar(mnames[i])), sep = "")
        print(mods[[i]], showEnv = FALSE)
      }
    }

    cat("\t qmodel:\n")
    printFormList(qmodel(pars(object)))
    cat("\t vmodel:\n")
    printFormList(vmodel(pars(object)))
   
 })






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
setGeneric("a4aFitSA", function(object, ...) standardGeneric("a4aFitSA"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("a4aFitSA", signature(object="missing"),
  function(...) {
    # empty
  	if(missing(...)){
	  	new("a4aFitSA")
    # or not
  	} else {
      args <- list(...)
	  args$Class <- 'a4aFitSA'
      do.call("new", args)
	  }
  }
)

setMethod("a4aFit", signature(object="a4aFitSA"),
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
    out @ fitSumm <- object @ fitSumm
    out
  }
)


setMethod("a4aFitSA", signature(object="a4aFit"),
  function(object, ...) {
    out <- a4aFitSA()
    out @ name    <- object @ name
    out @ desc    <- object @ desc
    out @ range   <- object @ range
    out @ call    <- object @ call
    out @ clock   <- object @ clock
    out @ stock.n <- object @ stock.n
    out @ harvest <- object @ harvest
    out @ catch.n <- object @ catch.n
    out @ index   <- object @ index
    out @ fitSumm <- object @ fitSumm
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
setGeneric("fitSumm", function(object, ...) standardGeneric("fitSumm"))

#' Title 
#' @name fitSumm
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("fitSumm", "a4aFitSA", function(object) object@fitSumm)

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
#' @rdname pars-methods
#'
#' @examples
#' data(ple4)
setGeneric("pars", function(object, ...) standardGeneric("pars"))

#' Title 
#' @name pars
#' @docType methods
#' @rdname pars-methods
#' @aliases pars,a4aFitSA-method
setMethod("pars", "a4aFitSA", function(object) object@pars)


