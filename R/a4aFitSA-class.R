#' @title S4 class \code{a4aFitSA}
#'
#' @description The \code{a4aFitSA} class extends \code{a4aFit} to store information about the parameters of the model.
#'
#' @section Slots:
#' \describe{
#'    \item{call}{The function call}
#'    \item{clock}{Information on call duration}
#'    \item{fitSumm}{Fit summary}
#'    \item{stock.n}{Estimates of stock numbers-at-age}
#'    \item{harvest}{Estimates of fishing mortality at age}
#'    \item{catch.n}{Estimates of catch numbers-at-age}
#'    \item{index}{Estimates of survey or CPUE indices-at-age}
#'    \item{pars}{an object of class \code{SCAPars} with information about model parameters}
#' }
#'
#' @template Accessors
#' @template Constructors
#' @docType class
#' @name a4aFitSA-class
#' @rdname a4aFitSA-class
#' @aliases a4aFitSA-class
#' @examples
#' data(ple4)
#' data(ple4.index)
#'
#' obj <- sca(stock=ple4, indices=FLIndices(ple4.index), fit="assessment")
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
#'
#' is(pars(obj))

a4aFitSA <-
  setClass("a4aFitSA",
           contains = "a4aFit",
           slots = c(pars = "SCAPars"))

#' @rdname a4aFitSA-class
#' @template bothargs
#' @aliases a4aFitSA a4aFitSA-methods
setGeneric("a4aFitSA")

#' @rdname a4aFitSA-class
setMethod("a4aFit", "a4aFitSA",
  function(object, ...) {
    as(object, "a4aFit")
  }
)

setMethod("initialize", "a4aFitSA",
    function(.Object,
             ...,
             pars = SCAPars()) {
      .Object <- callNextMethod(.Object, ...)
      .Object@pars <- pars
      .Object
})

setValidity("a4aFitSA",
  function(object) {
    # no validation at present
    TRUE
})

#
#  accessor methods
#

#' @rdname a4aFitSA-class
#' @aliases pars pars-methods
setGeneric("pars", function(object) standardGeneric("pars"))
#' @rdname a4aFitSA-class
setMethod("pars", "a4aFitSA", function(object) object@pars)

#' @rdname a4aFitSA-class
setMethod("m", signature(object = "a4aFitSA"), function(object) m(pars(object)))

#' @rdname a4aFitSA-class
setMethod("wt", signature(object = "a4aFitSA"),
  function(object) wt(pars(object)))

#' @rdname a4aFitSA-class
setMethod("qmodel", signature(object = "a4aFitSA"),
  function(object) qmodel(pars(object)))
setMethod("qMod", signature(object = "a4aFitSA"),
  function(object) qmodel(pars(object)))

#' @rdname a4aFitSA-class
setMethod("fmodel", signature(object = "a4aFitSA"),
  function(object) fmodel(pars(object)))
setMethod("fMod", signature(object = "a4aFitSA"),
  function(object) fmodel(pars(object)))

#' @rdname a4aFitSA-class
setMethod("srmodel", signature(object = "a4aFitSA"),
  function(object) srmodel(pars(object)))
setMethod("srMod", signature(object = "a4aFitSA"),
  function(object) srmodel(pars(object)))

#' @rdname a4aFitSA-class
setMethod("n1model", signature(object = "a4aFitSA"),
  function(object) n1model(pars(object)))
setMethod("n1Mod", signature(object = "a4aFitSA"),
  function(object) n1model(pars(object)))

#' @rdname a4aFitSA-class
setMethod("vmodel", signature(object = "a4aFitSA"),
  function(object) vmodel(pars(object)))
setMethod("vMod", signature(object = "a4aFitSA"),
  function(object) vmodel(pars(object)))

#' @rdname a4aFitSA-class
setMethod("stkmodel", signature(object = "a4aFitSA"),
  function(object) stkmodel(pars(object)))



#' @rdname a4aFitSA-class
setMethod("show", "a4aFitSA",
  function(object){
    show(a4aFit(object))

    cat("\nSubmodels:\n")
    show(stkmodel(object))
    cat("\n")
    show(qmodel(object))
    cat("\n")
    show(vmodel(object))
 })


#
#  other methods
#


#' @rdname a4aFitSA-class
#' @param obj the object to be subset
#' @param it iteration to be extracted
setMethod("iter", "a4aFitSA", function(obj, it){
	obj@fitSumm <- obj@fitSumm[,it, drop=FALSE]
	obj@harvest <- iter(obj@harvest, it)
	obj@stock.n <- iter(obj@stock.n, it)
	obj@catch.n <- iter(obj@catch.n, it)
	obj@index <- iter(obj@index, it)
	obj@pars <- iter(obj@pars, 1)
	obj
})
