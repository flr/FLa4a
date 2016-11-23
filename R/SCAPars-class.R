#' @title Model parameters class
#' @docType class
#' @name SCAPars
#' @rdname SCAPars-class
#' @template ClassDescription
#' @section Slot: 
#' \describe{
#'
#'	\item{\code{stkmodel}}{parameters related to stock dynamics}
#'
#'	\item{\code{qmodel}}{paramaters related to catchability of tunning fleets}
#'
#'	\item{\code{vmodel}}{paramaters related to the variance model}
#'
#' }
#' @aliases SCAPars-class
setClass("SCAPars",
        representation(
               stkmodel    = "a4aStkParams",
                 qmodel    = "submodels",
                 vmodel    = "submodels"),
        prototype = prototype(
               stkmodel    = new("a4aStkParams"),
                 qmodel    = new("submodels"),
                 vmodel    = new("submodels"))
)

#' @rdname SCAPars-class
#' @template Accessors
#' @template Constructors
#' @template bothargs
#' @aliases SCAPars SCAPars-methods
setGeneric("SCAPars", function(object, ...) standardGeneric("SCAPars"))
#' @rdname SCAPars-class
setMethod("SCAPars", signature(object="missing"),
  function(...) {
    # empty
  	if(missing(...)){
	  	new("SCAPars")
    # or not
  	} else {
      args <- list(...)
	  args$Class <- 'SCAPars'
      do.call("new", args)
	  }
  }
)

# accessors

#' @rdname SCAPars-class
#' @aliases stkmodel stkmodel-methods
setGeneric("stkmodel", function(object, ...) standardGeneric("stkmodel"))
#' @rdname SCAPars-class
setMethod("stkmodel", "SCAPars", function(object) object @ stkmodel)

#' @rdname SCAPars-class
#' @aliases n1model n1model-methods
setGeneric("n1model", function(object, ...) standardGeneric("n1model"))
#' @rdname SCAPars-class
setMethod("n1model", "SCAPars", function(object) object @ stkmodel @ n1Mod)

#' @rdname SCAPars-class
#' @aliases srmodel srmodel-methods
setGeneric("srmodel", function(object, ...) standardGeneric("srmodel"))
#' @rdname SCAPars-class
setMethod("srmodel", "SCAPars", function(object) object @ stkmodel @ srMod)

#' @rdname SCAPars-class
#' @aliases fmodel fmodel-methods
setGeneric("fmodel", function(object, ...) standardGeneric("fmodel"))
#' @rdname SCAPars-class
setMethod("fmodel", "SCAPars", function(object) object @ stkmodel @ fMod)

#' @rdname SCAPars-class
#' @aliases qmodel qmodel-methods
setGeneric("qmodel", function(object, ...) standardGeneric("qmodel"))
#' @rdname SCAPars-class
setMethod("qmodel", "SCAPars", function(object) object@qmodel)

#' @rdname SCAPars-class
#' @aliases vmodel vmodel-methods
setGeneric("vmodel", function(object, ...) standardGeneric("vmodel"))
#' @rdname SCAPars-class
setMethod("vmodel", "SCAPars", function(object) object@vmodel)

#' @rdname SCAPars-class
#' @aliases srPars srPars-methods
setGeneric("srPars", function(object, ...) standardGeneric("srPars"))
#' @rdname SCAPars-class
setMethod("srPars", "SCAPars", function(object) object@srmodel@pars)

#' @rdname SCAPars-class
#' @aliases srCovar srCovar-methods
setGeneric("srCovar", function(object, ...) standardGeneric("srCovar"))
#' @rdname SCAPars-class
setMethod("srCovar", "SCAPars", function(object) object@srmodel@covar)

#' @rdname SCAPars-class
#' @aliases srFrml srFrml-methods
setGeneric("srFrml", function(object, ...) standardGeneric("srFrml"))
#' @rdname SCAPars-class
setMethod("srFrml", "SCAPars", function(object) object@srmodel@model)

#' @rdname SCAPars-class
#' @aliases fPars fPars-methods
setGeneric("fPars", function(object, ...) standardGeneric("fPars"))
#' @rdname SCAPars-class
setMethod("fPars", "SCAPars", function(object) object@fmodel@pars)

#' @rdname SCAPars-class
#' @aliases fCovar fCovar-methods
setGeneric("fCovar", function(object, ...) standardGeneric("fCovar"))
#' @rdname SCAPars-class
setMethod("fCovar", "SCAPars", function(object) object@fmodel@covar)

#' @rdname SCAPars-class
#' @aliases fFrml fFrml-methods
setGeneric("fFrml", function(object, ...) standardGeneric("fFrml"))
#' @rdname SCAPars-class
setMethod("fFrml", "SCAPars", function(object) object@fmodel@model)

#' @rdname SCAPars-class
#' @aliases qPars qPars-methods
setGeneric("qPars", function(object, ...) standardGeneric("qPars"))
#' @rdname SCAPars-class
setMethod("qPars", "SCAPars", function(object) object@qmodel@pars)

#' @rdname SCAPars-class
#' @aliases qCovar qCovar-methods
setGeneric("qCovar", function(object, ...) standardGeneric("qCovar"))
#' @rdname SCAPars-class
setMethod("qCovar", "SCAPars", function(object) object@qmodel@covar)

#' @rdname SCAPars-class
#' @aliases qFrml qFrml-methods
setGeneric("qFrml", function(object, ...) standardGeneric("qFrml"))
#' @rdname SCAPars-class
setMethod("qFrml", "SCAPars", function(object) object@qmodel@model)

#' @rdname SCAPars-class
#' @aliases vPars vPars-methods
setGeneric("vPars", function(object, ...) standardGeneric("vPars"))
#' @rdname SCAPars-class
setMethod("vPars", "SCAPars", function(object) object@vmodel@pars)

#' @rdname SCAPars-class
#' @aliases vCovar vCovar-methods
setGeneric("vCovar", function(object, ...) standardGeneric("vCovar"))
#' @rdname SCAPars-class
setMethod("vCovar", "SCAPars", function(object) object@vmodel@covar)

#' @rdname SCAPars-class
#' @aliases vFrml vFrml-methods
setGeneric("vFrml", function(object, ...) standardGeneric("vFrml"))
#' @rdname SCAPars-class
setMethod("vFrml", "SCAPars", function(object) object@vmodel@model)

#' @rdname SCAPars-class
setMethod("m", signature(object="SCAPars"), function(object) m(stkmodel(object)))

#' @rdname SCAPars-class
setMethod("wt", signature(object="SCAPars"), function(object) wt(stkmodel(object)))


