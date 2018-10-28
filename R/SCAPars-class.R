#' @title Model parameters class
#' @docType class
#' @name SCAPars
#' @rdname SCAPars-class
#' @template ClassDescription
#' @section Slot:
#' \describe{
#'  \item{\code{stkmodel}}{parameters related to stock dynamics}
#'  \item{\code{qmodel}}{paramaters related to catchability of tunning fleets}
#'  \item{\code{vmodel}}{paramaters related to the variance model}
#' }
#' @aliases SCAPars-class
setClass("SCAPars",
        slots = c(stkmodel  = "a4aStkParams",
                  qmodel    = "submodels",
                  vmodel    = "submodels")
)

setValidity("SCAPars",
  function(object) {
    # vmodel should be 1 longer than qmodel
    if (length(object@qmodel) + 1 != length(object@vmodel)) {
      "vmodel should be 1 longer than qmodel"
    } else {
      TRUE
    }
})

setMethod("initialize", "SCAPars",
  function(.Object,
           stkmodel = new("a4aStkParams"),
           qmodel = new("submodels", name = "qmodel"),
           vmodel = new("submodels", name = "vmodel")) {
      # initialize FLComp slots
      .Object <- callNextMethod()
      # ensure correct names in q nd v models
      qmodel@name <- "qmodel"
      vmodel@name <- "vmodel"
      # initialize remaining slots
      .Object@stkmodel <- stkmodel
      .Object@qmodel <- qmodel
      .Object@vmodel <- vmodel
      .Object
})




#' @rdname SCAPars-class
#' @template Accessors
#' @template Constructors
#' @template bothargs
#' @aliases SCAPars SCAPars-methods
setGeneric("SCAPars", function(object, ...) standardGeneric("SCAPars"))
#' @rdname SCAPars-class
setMethod("SCAPars", signature(object = "missing"),
  function(...) {
    # empty
    if (missing(...)) {
      new("SCAPars")
    # or not
    } else {
      args <- list(...)
    args$Class <- "SCAPars"
      do.call("new", args)
    }
  }
)

# accessors

#' @rdname SCAPars-class
#' @aliases stkmodel stkmodel-methods
setGeneric("stkmodel", function(object, ...) standardGeneric("stkmodel"))
#' @rdname SCAPars-class
setMethod("stkmodel", "SCAPars", function(object) object@stkmodel)

#' @rdname SCAPars-class
#' @aliases n1model n1model-methods
setGeneric("n1model", function(object, ...) standardGeneric("n1model"))
#' @rdname SCAPars-class
setMethod("n1model", "SCAPars", function(object) n1Mod(stkmodel(object)))
setMethod("n1Mod", "SCAPars", function(object) n1Mod(stkmodel(object)))

#' @rdname SCAPars-class
#' @aliases srmodel srmodel-methods
setGeneric("srmodel", function(object, ...) standardGeneric("srmodel"))
#' @rdname SCAPars-class
setMethod("srmodel", "SCAPars", function(object) srMod(stkmodel(object)))
setMethod("srMod", "SCAPars", function(object) srMod(stkmodel(object)))

#' @rdname SCAPars-class
setMethod("fmodel", "SCAPars", function(object) fmodel(stkmodel(object)))
setMethod("fMod", "SCAPars", function(object) fMod(stkmodel(object)))

#' @rdname SCAPars-class
#' @aliases qmodel qmodel-methods
setGeneric("qmodel", function(object, ...) standardGeneric("qmodel"))
#' @rdname SCAPars-class
setMethod("qmodel", "SCAPars", function(object) object@qmodel)

#' @rdname SCAPars-class
#' @aliases qMod qMod-methods
setGeneric("qMod", function(object, ...) standardGeneric("qMod"))
#' @rdname SCAPars-class
setMethod("qMod", "SCAPars", function(object) qmodel(object))

#' @rdname SCAPars-class
#' @aliases vmodel vmodel-methods
setGeneric("vmodel", function(object, ...) standardGeneric("vmodel"))
#' @rdname SCAPars-class
setMethod("vmodel", "SCAPars", function(object) object@vmodel)

#' @rdname SCAPars-class
#' @aliases vMod vMod-methods
setGeneric("vMod", function(object, ...) standardGeneric("vMod"))
#' @rdname SCAPars-class
setMethod("vMod", "SCAPars", function(object) vmodel(object))


## NOTE:  srPars method not possible with present structure
#' @rdname SCAPars-class
#' @aliases srPars srPars-methods
setGeneric("srPars", function(object, ...) standardGeneric("srPars"))
#' @rdname SCAPars-class
setMethod("srPars", "SCAPars", function(object) object@srmodel@pars)

## NOTE:  srPars method not possible with present structure
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
setMethod("fPars", "SCAPars", function(object) {
  stkmodel_coefficients <- coef(stkmodel(object))
  idx <- grep("fMod:", rownames(stkmodel_coefficients))
  stkmodel_coefficients[idx, ]
})

#' @rdname SCAPars-class
#' @aliases fCovar fCovar-methods
setGeneric("fCovar", function(object, ...) standardGeneric("fCovar"))
#' @rdname SCAPars-class
setMethod("fCovar", "SCAPars", function(object) object@fmodel@covar)

#' @rdname SCAPars-class
#' @aliases fFrml fFrml-methods
setGeneric("fFrml", function(object, ...) standardGeneric("fFrml"))
#' @rdname SCAPars-class
setMethod("fFrml", "SCAPars", function(object) fmodel(object))

#' @rdname SCAPars-class
#' @aliases qPars qPars-methods
setGeneric("qPars", function(object, ...) standardGeneric("qPars"))
#' @rdname SCAPars-class
setMethod("qPars", "SCAPars", function(object) object@qmodel@pars)

#' @rdname SCAPars-class
#' @aliases qCovar qCovar-methods
setGeneric("qCovar", function(object, ...) standardGeneric("qCovar"))
#' @rdname SCAPars-class
setMethod("qCovar", "SCAPars", function(object) vcov(qmodel(object)))

#' @rdname SCAPars-class
#' @aliases qFrml qFrml-methods
setGeneric("qFrml", function(object, ...) standardGeneric("qFrml"))
#' @rdname SCAPars-class
setMethod("qFrml", "SCAPars", function(object) sMod(qmodel(object)))

#' @rdname SCAPars-class
#' @aliases vPars vPars-methods
setGeneric("vPars", function(object, ...) standardGeneric("vPars"))
#' @rdname SCAPars-class
setMethod("vPars", "SCAPars", function(object) object@vmodel@pars)

#' @rdname SCAPars-class
#' @aliases vCovar vCovar-methods
setGeneric("vCovar", function(object, ...) standardGeneric("vCovar"))
#' @rdname SCAPars-class
setMethod("vCovar", "SCAPars", function(object) vcov(vmodel(object)))

#' @rdname SCAPars-class
#' @aliases vFrml vFrml-methods
setGeneric("vFrml", function(object, ...) standardGeneric("vFrml"))
#' @rdname SCAPars-class
setMethod("vFrml", "SCAPars", function(object) object@vmodel@model)

#' @rdname SCAPars-class
setMethod("m", signature(object = "SCAPars"), function(object) m(stkmodel(object)))

#' @rdname SCAPars-class
setMethod("wt", signature(object = "SCAPars"), function(object) wt(stkmodel(object)))





#' @rdname SCAPars-class
#' @param iter the number of iterations to create
#' @param fill.iter should the new iterations be filled with values (TRUE) or NAs (FALSE)
setMethod("propagate",
  signature(object = "SCAPars"),
  function (object, iter, fill.iter = TRUE) {
    # stkmodel
    object@stkmodel <- propagate(object@stkmodel, iter, fill.iter = fill.iter)
    # qmodel
    object@qmodel <- propagate(object@qmodel, iter, fill.iter = fill.iter)
    # vmodel
    object@vmodel <- propagate(object@vmodel, iter, fill.iter = fill.iter)

    object
  }
)


#' @rdname SCAPars-class
#' @param obj the object to be subset
#' @param it iteration to be extracted
setMethod("iter", "SCAPars", function(obj, it){
  obj@stkmodel <- iter(obj@stkmodel, it)
  obj@qmodel <- iter(obj@qmodel, it)
  obj@vmodel <- iter(obj@vmodel, it)
  obj
})
