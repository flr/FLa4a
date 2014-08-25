#' @title Model parameters class
#' @docType class
#' @name SCAPars
#' @rdname SCAPars-class
#' @template ClassDescription
#' @section Slot: 
#' \describe{
#'
#'	\item{\code{stkmodel}}{Parameters related with the stock dynamics}
#'
#'	\item{\code{qmodel}}{Paramaters related with catchability of the tunning fleets}
#'
#'	\item{\code{vmodel}}{Paramaters related with the variance model}
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
#' @aliases SCAPars SCAPars-methods SCAPars,missing-method
setGeneric("SCAPars", function(object, ...) standardGeneric("SCAPars"))
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
#' @aliases stkmodel stkmodel-methods stkmodel,SCAPars-method
setGeneric("stkmodel", function(object, ...) standardGeneric("stkmodel"))
setMethod("stkmodel", "SCAPars", function(object) object @ stkmodel)

#' @rdname SCAPars-class
#' @aliases n1model n1model-methods n1model,SCAPars-method
setGeneric("n1model", function(object, ...) standardGeneric("n1model"))
setMethod("n1model", "SCAPars", function(object) object @ stkmodel @ n1Mod)

#' @rdname SCAPars-class
#' @aliases srmodel srmodel-methods srmodel,SCAPars-method
setGeneric("srmodel", function(object, ...) standardGeneric("srmodel"))
setMethod("srmodel", "SCAPars", function(object) object @ stkmodel @ srMod)

#' @rdname SCAPars-class
#' @aliases fmodel fmodel-methods fmodel,SCAPars-method
setGeneric("fmodel", function(object, ...) standardGeneric("fmodel"))
setMethod("fmodel", "SCAPars", function(object) object @ stkmodel @ fMod)

#' @rdname SCAPars-class
#' @aliases qmodel qmodel-methods qmodel,SCAPars-method
setGeneric("qmodel", function(object, ...) standardGeneric("qmodel"))
setMethod("qmodel", "SCAPars", function(object) object@qmodel)

#' @rdname SCAPars-class
#' @aliases vmodel vmodel-methods vmodel,SCAPars-method
setGeneric("vmodel", function(object, ...) standardGeneric("vmodel"))
setMethod("vmodel", "SCAPars", function(object) object@vmodel)

#' @rdname SCAPars-class
#' @aliases srPars srPars-methods srPars,SCAPars-method
setGeneric("srPars", function(object, ...) standardGeneric("srPars"))
setMethod("srPars", "SCAPars", function(object) object@srmodel@pars)

#' @rdname SCAPars-class
#' @aliases srCovar srCovar-methods srCovar,SCAPars-method
setGeneric("srCovar", function(object, ...) standardGeneric("srCovar"))
setMethod("srCovar", "SCAPars", function(object) object@srmodel@covar)

#' @rdname SCAPars-class
#' @aliases srFrml srFrml-methods srFrml,SCAPars-method
setGeneric("srFrml", function(object, ...) standardGeneric("srFrml"))
setMethod("srFrml", "SCAPars", function(object) object@srmodel@model)

#' @rdname SCAPars-class
#' @aliases fPars fPars-methods fPars,SCAPars-method
setGeneric("fPars", function(object, ...) standardGeneric("fPars"))
setMethod("fPars", "SCAPars", function(object) object@fmodel@pars)

#' @rdname SCAPars-class
#' @aliases fCovar fCovar-methods fCovar,SCAPars-method
setGeneric("fCovar", function(object, ...) standardGeneric("fCovar"))
setMethod("fCovar", "SCAPars", function(object) object@fmodel@covar)

#' @rdname SCAPars-class
#' @aliases fFrml fFrml-methods fFrml,SCAPars-method
setGeneric("fFrml", function(object, ...) standardGeneric("fFrml"))
setMethod("fFrml", "SCAPars", function(object) object@fmodel@model)

#' @rdname SCAPars-class
#' @aliases qPars qPars-methods qPars,SCAPars-method
setGeneric("qPars", function(object, ...) standardGeneric("qPars"))
setMethod("qPars", "SCAPars", function(object) object@qmodel@pars)

#' @rdname SCAPars-class
#' @aliases qCovar qCovar-methods qCovar,SCAPars-method
setGeneric("qCovar", function(object, ...) standardGeneric("qCovar"))
setMethod("qCovar", "SCAPars", function(object) object@qmodel@covar)

#' @rdname SCAPars-class
#' @aliases qFrml qFrml-methods qFrml,SCAPars-method
setGeneric("qFrml", function(object, ...) standardGeneric("qFrml"))
setMethod("qFrml", "SCAPars", function(object) object@qmodel@model)

#' @rdname SCAPars-class
#' @aliases vPars vPars-methods vPars,SCAPars-method
setGeneric("vPars", function(object, ...) standardGeneric("vPars"))
setMethod("vPars", "SCAPars", function(object) object@vmodel@pars)

#' @rdname SCAPars-class
#' @aliases vCovar vCovar-methods vCovar,SCAPars-method
setGeneric("vCovar", function(object, ...) standardGeneric("vCovar"))
setMethod("vCovar", "SCAPars", function(object) object@vmodel@covar)

#' @rdname SCAPars-class
#' @aliases vFrml vFrml-methods vFrml,SCAPars-method
setGeneric("vFrml", function(object, ...) standardGeneric("vFrml"))
setMethod("vFrml", "SCAPars", function(object) object@vmodel@model)

#' @rdname SCAPars-class
#' @aliases m,SCAPars-method
setMethod("m", signature(object="SCAPars"), function(object) m(stkmodel(object)))

#' @rdname SCAPars-class
#' @aliases wt,SCAPars-method
setMethod("wt", signature(object="SCAPars"), function(object) wt(stkmodel(object)))


