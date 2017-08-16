#' @title S4 class \code{a4aFitMCMC}
#'
#' @description The \code{a4aFitMCMC} class extends \code{a4aFitSA} to store information about the MCMC run.
#'
#' @section Slots:
#' \describe{
#'    \item{name}{A character vector for the object name.}
#'    \item{desc}{A textual description of the object contents.}
#'    \item{range}{A named numeric vector with various values of quant and year ranges, plusgroup, fishing mortality ranges, etc.}
#'    \item{call}{The function call}
#'    \item{clock}{Information on call duration}
#'    \item{fitSumm}{Fit summary}
#'    \item{stock.n}{Estimates of stock numbers-at-age}
#'    \item{harvest}{Estimates of fishing mortality at age}
#'    \item{catch.n}{Estimates of catch numbers-at-age}
#'    \item{index}{Estimates of survey or CPUE indices-at-age}
#'    \item{mcmc}{An object of class \code{SCAMCMC} with information about the MCMC run}
#' }
#'
#' @template Accessors
#' @template Constructors
#' @docType class
#' @name a4aFitMCMC-class
#' @rdname a4aFitMCMC-class
#' @aliases a4aFitMCMC-class
#' @template Example-a4aFitSA
a4aFitMCMC <-
  setClass("a4aFitMCMC",
           contains = "a4aFitSA",
           slots = c(mcmc = "SCAMCMC"))

#' @rdname a4aFitMCMC-class
#' @template bothargs
#' @aliases a4aFitMCMC a4aFitMCMC-methods
setGeneric("a4aFitMCMC")

#' @rdname a4aFitMCMC-class
setMethod("a4aFitSA", "a4aFitMCMC",
  function(object, ...) {
    as(object, "a4aFitSA")
  }
)

#' @rdname a4aFitMCMC-class
setMethod("a4aFit", "a4aFitMCMC",
  function(object, ...) {
    as(object, "a4aFit")
  }
)

setMethod("initialize", "a4aFitMCMC",
    function(.Object, ..., mcmc) {
      if (!missing(mcmc)) .Object@mcmc <- mcmc
      .Object <- callNextMethod(.Object, ...)
      .Object
})




#====================================================================
# coerce to coda object
#====================================================================
#' @rdname a4aFitMCMC-class
#' @param x an object to be coerced into mcmc
#' @aliases as.mcmc as.mcmc-methods
setGeneric("as.mcmc", function(x, ...) useAsDefault=coda::as.mcmc)

#' @rdname a4aFitMCMC-class
setMethod("as.mcmc", signature(x="a4aFitMCMC"), function(x, ...) {
		object <- x
		df0 <- t(object@pars@stkmodel@params[drop=T])
		lst <- lapply(object@pars@qmodel, par2mat)
		df1 <- do.call("cbind", lst)		
		lst <- lapply(object@pars@vmodel, par2mat)
		df2 <- do.call("cbind", lst)		
		df0 <- cbind(df0, df1, df2)
		mcmc(df0)
})

#====================================================================
# The methods below should be inherited from a4aFitSA
#====================================================================

##' @rdname a4aFitSA-class
##' @aliases show,a4aFitSA-method
#setMethod("show", signature(object = "a4aFitSA"),
#  function(object) 
#  {

#    show(a4aFit(object))
#    
#    cat("\nSubmodels:\n")  
#    cat("\t fmodel: "); print( fmodel(pars(object)), showEnv = FALSE)
#    cat("\tsrmodel: "); print(srmodel(pars(object)), showEnv = FALSE)
#    cat("\tn1model: "); print(n1model(pars(object)), showEnv = FALSE)

#    # something to format the qmodel and vmodel    
#    printFormList <- function(frmL) {
#      if (length(frmL) == 0) return(invisible(NA))
#      mods <- lapply(frmL, slot, "Mod")
#      mnames <- names(mods)
#      maxname <- max(sapply(mnames, nchar))
#      for (i in seq(length(mods))) {
#        cat("\t   ", mnames[i], ": ", rep(" ", maxname - nchar(mnames[i])), sep = "")
#        print(mods[[i]], showEnv = FALSE)
#      }
#    }

#    cat("\t qmodel:\n")
#    printFormList(qmodel(pars(object)))
#    cat("\t vmodel:\n")
#    printFormList(vmodel(pars(object)))
#   
# })

##' @rdname a4aFitSA-class
##' @aliases submodels,a4aFitSA-method
#setMethod("submodels", signature(object = "a4aFitSA"),
#  function(object) 
#  {
#    cat("\t fmodel: "); print( fmodel(pars(object)), showEnv = FALSE)
#    cat("\tsrmodel: "); print(srmodel(pars(object)), showEnv = FALSE)
#    cat("\tn1model: "); print(n1model(pars(object)), showEnv = FALSE)

#    # something to format the qmodel and vmodel    
#    printFormList <- function(frmL) {
#      if (length(frmL) == 0) return(invisible(NA))
#      mods <- lapply(frmL, slot, "Mod")
#      mnames <- names(mods)
#      maxname <- max(sapply(mnames, nchar))
#      for (i in seq(length(mods))) {
#        cat("\t   ", mnames[i], ": ", rep(" ", maxname - nchar(mnames[i])), sep = "")
#        print(mods[[i]], showEnv = FALSE)
#      }
#    }

#    cat("\t qmodel:\n")
#    printFormList(qmodel(pars(object)))
#    cat("\t vmodel:\n")
#    printFormList(vmodel(pars(object)))
#   
# })
#
##' @rdname a4aFitSA-class
##' @aliases pars pars-method pars,a4aFitSA-method
#setGeneric("pars", function(object) standardGeneric("pars"))
#setMethod("pars", "a4aFitSA", function(object) object@pars)

##' @rdname a4aFitSA-class
##' @aliases m,a4aFitSA-method
#setMethod("m", signature(object="a4aFitSA"), function(object) m(pars(object)))

##' @rdname a4aFitSA-class
##' @aliases wt,a4aFitSA-method
#setMethod("wt", signature(object="a4aFitSA"), function(object) wt(pars(object)))

##====================================================================
## plural class for a4aFitSA (used for model averaging)
##====================================================================

##' @rdname a4aFitSA-class
##' @aliases a4aFitSAs-class

#setClass("a4aFitSAs", 
#	contains="FLComps",
#	validity=function(object){
#		if(!all(unlist(lapply(object, is, 'a4aFitSA'))))
#			return("Components must be a4aFitSA")	
#		return(TRUE)}
#)

##' @rdname a4aFitSA-class
##' @aliases a4aFitSAs a4aFitSAs,list-method
#setGeneric("a4aFitSAs", function(object, ...) standardGeneric("a4aFitSAs"))
#setMethod("a4aFitSAs", signature(object="list"),
#  function(object, ...) {
#    args <- list(...)
#    
#    # names in args, ... 
#    if("names" %in% names(args)) {
#      names <- args[['names']]
#    } else {
#    # ... or in object,
#      if(!is.null(names(object))) {
#        names <- names(object)
#    # ... or in elements, ...
#      } else {
#        names <- unlist(lapply(object, name))
#        # ... or 1:n
#        idx <- names == "NA" | names == ""
#        if(any(idx))
#          names[idx] <- as.character(length(names))[idx]
#      }
#    }

#    # desc & lock
#    args <- c(list(Class="a4aFitSAs", .Data=object, names=names),
#      args[!names(args)%in%'names'])

#    return(
#      do.call('new', args)
#      )

#})

##' @rdname a4aFitSA-class
##' @aliases a4aFitSAs,a4aFitSA-method
#setMethod("a4aFitSAs", signature(object="a4aFitSA"), function(object, ...) {
#    lst <- c(object, list(...))
#    a4aFitSAs(lst)
#})

##' @rdname a4aFitSA-class
##' @aliases a4aFitSAs a4aFitSAs,missing-method
#setMethod("a4aFitSAs", signature(object="missing"),
#  function(...) {
#    # empty
#  	if(missing(...)){
#	  	new("a4aFitSAs")
#    # or not
#  	} else {
#      args <- list(...)
#      object <- args[!names(args)%in%c('names', 'desc', 'lock')]
#      args <- args[!names(args)%in%names(object)]
#      do.call('a4aFitSAs',  c(list(object=object), args))
#	  }
#  }
#)

