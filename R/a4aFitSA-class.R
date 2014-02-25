#' @title S4 class \code{a4aFitSA}
#'
#' @description The \code{a4aFitSA} class extends \code{a4aFit} to store information about the parameters of the model.
#'
#' @section Slots:
#' \describe{
#'    \item{SCAPars}{An object of class \code{SCAPars} with information about model parameters}
#'
#' @template Accessors
#' @template Constructors
#' @docType class
#' @name a4aFitSA-class
#' @rdname a4aFitSA-class
#' @aliases a4aFitSA-class
#' @template Example-a4aFitSA
setClass("a4aFitSA",
        representation(
                "a4aFit",
                pars    = "SCAPars"
                ),
        prototype = prototype(
                pars    = new('SCAPars'))
)

setMethod("m", signature(object="a4aFitSA"),
  function(object) m(pars(object)))


#' @rdname a4aFitSA-class
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

#' @rdname a4aFitSA-class
#' @aliases a4aFitSA a4aFitSA-methods
setGeneric("a4aFitSA", function(object, ...) standardGeneric("a4aFitSA"))

#' @rdname a4aFitSA-class
#' @aliases a4aFitSA,missing-method
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

#' @rdname a4aFitSA-class
#' @aliases a4aFit,a4aFitSA-method
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

#' @rdname a4aFitSA-class
#' @aliases a4aFitSA,a4aFit-method
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

#' @rdname a4aFitSA-class
#' @aliases pars,a4aFitSA-method
setGeneric("pars", function(object, ...) standardGeneric("pars"))
setMethod("pars", "a4aFitSA", function(object) object@pars)


