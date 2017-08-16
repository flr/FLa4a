#' @title Stock model class
#' @docType class
#' @name a4aStkModel
#' @rdname a4aStkModel-class
#' @template ClassDescription
#' @section Slot:
#' \describe{
#'  \item{\code{fmodel}}{F submodel \code{formula}}
#'  \item{\code{n1model}}{first year N \code{formula}}
#'  \item{\code{srmodel}}{stock-recruitment submodel \code{formula}}
#'  \item{\code{m}}{natural mortality \code{FLQuant}}
#'  \item{\code{wt}}{stock weights \code{FLQuant}}
#'  \item{\code{units}}{data units \code{character}}
#' }
#' @aliases a4aStkModel-class

a4aStkModel <-
  setClass("a4aStkModel",
    contains = c("FLComp", "submodels"),
    slots = c(vcov = "array",
              m = "FLQuant",
              wt = "FLQuant",
              units = "character"))

#' @rdname a4aStkModel-class
#' @aliases a4aStkModel a4aStkModel-methods
#' @template Accessors
#' @template Constructors
setGeneric("a4aStkModel")

setMethod("initialize", "a4aStkModel",
    function(.Object,
             ...,
             vcov,
             m, wt, units) {
      # fill up FLComp slots first
      .Object <- callNextMethod(.Object, ...)
      if (length(.Object) == 0) {
        .Object@.Data <- list(f = submodel(name = "f"),
                              n1 = submodel(name = "n1"),
                              sr = submodel(name = "sr"))
      }
      if (length(.Object) < 3) {
        stop("There must be 3 submodels named f, n1, and sr")
      }
      # if missing set dimensions of of m and wt based on range slot
      if (missing(m) || missing(wt)) {
        flq <- FLQuant(
                  matrix(NA,
                     nrow = .Object@range["max"] - .Object@range["min"] + 1,
                     ncol = .Object@range["maxyear"] - .Object@range["minyear"] + 1),
                     dimnames = list(age = .Object@range["min"]:.Object@range["max"],
                                     year = .Object@range["minyear"]:.Object@range["maxyear"])
                  )
      }
      .Object@m <- if (missing(m)) flq else m
      .Object@wt <- if (missing(wt)) flq else wt
      # throw error if range from FLComp doesn't match FLQuants 
      # (can't check this in setValidity due to callNextMethod resulting in an invalid a4aStkParams object when range is supplied)
      if (abs(as.numeric(dimnames(.Object@m)$year[1]) - .Object@range["minyear"]) > 1e-9 ||
          abs(as.numeric(dimnames(.Object@m)$year[dim(.Object@m)[2]]) - .Object@range["maxyear"]) > 1e-9) {
            stop("range does not match supplied m and wt dimensions")
      }
      if (!missing(vcov)) .Object@vcov <- vcov
      if (!missing(units)).Object@units <- units
      .Object
    })

setValidity("a4aStkModel",
  function(object) {
    # the following is not throwing errors somehow...
    #if (length(object) == 0) {
    #  "There must be 3 submodels named f, n1, and sr"
    #} else
    # check dimensions of m and wt are the same
    if (!identical(unlist(dimnames(object@m)[2:5]),
                   unlist(dimnames(object@wt)[2:5]))) {
      "m and wt elements must share dimensions 2 to 5"
    } else {
      TRUE
    }
})


#
#  accessor methods
#


#' @rdname a4aStkModel-class
setMethod("m", "a4aStkModel", function(object) object@m)

#' @rdname a4aStkModel-class
setMethod("wt", "a4aStkModel", function(object) object@wt)

#' @rdname a4aStkModel-class
setMethod("vcov", "a4aStkModel", function(object) object@vcov)


#
#  show methods
#

#' @rdname a4aStkModel-class
setMethod("show", "a4aStkModel",
  function(object)
  {
    cat("a4a stock model for:", object@name, "\n")
    show(submodels(object))
 })


#
# Other methods
#

#' @rdname a4aStkParams-class
#' @param obj the object to be subset
#' @param it iteration to be extracted
setMethod("iter", "a4aStkParams", function(obj, it){
  obj@vcov <- obj@vcov[,,it, drop=FALSE]
  obj@m <- iter(obj@m, it)
  obj@wt <- iter(obj@wt, it)
  obj
})

