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
		df0 <- t(object@pars@stkmodel@coefficients[drop=T])
		lst <- lapply(object@pars@qmodel, par2mat)
		df1 <- do.call("cbind", lst)		
		lst <- lapply(object@pars@vmodel, par2mat)
		df2 <- do.call("cbind", lst)		
		df0 <- cbind(df0, df1, df2)
		mcmc(df0)
})

#====================================================================
# burnin
#====================================================================
#' @rdname a4aFitMCMC-class
#' @param burnin a numeric with the number of iterations to be removed
#' @aliases burnin burnin-methods
setGeneric("burnin", function(object, ...) standardGeneric("burnin"))

#' @rdname a4aFitMCMC-class
setMethod("burnin", signature(object="a4aFitMCMC"), function(object, burnin){
	object@catch.n <- catch.n(object)[,,,,,-c(1:burnin)]
	object@stock.n <- stock.n(object)[,,,,,-c(1:burnin)]
	object@harvest <- harvest(object)[,,,,,-c(1:burnin)]
	object@index <- lapply(index(object), function(x) x[,,,,,-c(1:burnin)])
	object@pars@stkmodel@coefficients <- coefficients(stkmodel(pars(object)))[,-c(1:burnin)]
	object@pars@qmodel@.Data <- lapply(qmodel(pars(object)), function(x){
	    x@coefficients <- x@coefficients[,-c(1:burnin)]
	    x
	}) 
	object@pars@vmodel@.Data <- lapply(vmodel(pars(object)), function(x){
    	x@coefficients <- x@coefficients[,-c(1:burnin)]
    	x
	}) 
	object
})





