#' @title Stock parameters class
#' @docType class
#' @name a4aStkParams
#' @rdname a4aStkParams-class
#' @template ClassDescription
#' @section Slot:
#' \describe{
#'  \item{\code{fMod}}{F submodel \code{formula}}
#'  \item{\code{n1Mod}}{first year N \code{formula}}
#'  \item{\code{srMod}}{stock-recruitment submodel \code{formula}}
#'  \item{\code{params}}{\code{FLPar} with parameters}
#'  \item{\code{vcov}}{\code{array} with variance-covariance}
#'  \item{\code{centering}}{centering values \code{numeric}}
#'  \item{\code{distr}}{statistical distribution \code{character}}
#'  \item{\code{m}}{natural mortality \code{FLQuant}}
#'  \item{\code{units}}{data units \code{character}}
#' }
#' @aliases a4aStkParams-class
setClass("a4aStkParams_sim",
  contains = "a4aStkParams",
  slots = c(covariates   = "data.frame")
)

setMethod("initialize", "a4aStkParams_sim",
  function(.Object, 
           ..., 
           covariates = .Object@covariates) 
  {
    .Object <- callNextMethod(.Object, ...)

    # do assignments
    .Object@covariates <- covariates
    
    .Object
  })


#
# Coersion methods
#

# method.skeleton("coerce", "a4aStkParams",  file = stdout())

setMethod("coerce", signature(from = "a4aStkParams_sim", to = "submodels"),
  function (from, to, strict = TRUE) {
    callNextMethod()
  }
)
