#' @title Stock submodels class
#' @docType class
#' @name stk_submodels
#' @rdname stk_submodels-class
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

setClass(
  "stk_submodels",
  contains = c("FLComp", "submodels"),
  slots =
    c(
      "m" = "FLQuant",
      "stock.wt" = "FLQuant",
      "mat" = "FLQuant",
      "harvest.spwn" = "FLQuant",
      "m.spwn" = "FLQuant"
    )
)

#' @rdname a4aStkModel-class
#' @aliases a4aStkModel a4aStkModel-methods
#' @template Accessors
#' @template Constructors
setGeneric("stk_submodels")

setMethod(
  "initialize",
  "stk_submodels",
  function(.Object, ..., m, wt, mat) {
      # fill up FLComp slots first
      .Object <- callNextMethod(.Object, ...)
      if (length(.Object) == 0) {
        .Object@.Data <-
          list(
            f = submodel(name = "f"),
            n1 = submodel(name = "n1"),
            sr = submodel(name = "sr")
          )
      }
      if (length(.Object) < 3) {
        stop("There must be 3 submodels named f, n1, and sr")
      }
      # if missing set dimensions of of m and wt based on range slot
      if (missing(m) || missing(wt)) {
        range <- range(.Object)
        df <-
          expand.grid(
            age = range["min"]:range["max"],
            year = range["minyear"]:range["maxyear"],
            data = NA,
            stringsAsFactors = FALSE
          )
        flq <- as.FLQuant(df)
      }
      .Object@m <- if (missing(m)) flq else m
      .Object@wt <- if (missing(wt)) flq else wt
      # throw error if range from FLComp doesn't match FLQuants
      # (can't check this in setValidity due to callNextMethod resulting in an invalid a4aStkParams object when range is supplied)
      if (
        abs(as.numeric(dimnames(.Object@m)$year[1]) -
          .Object@range["minyear"]) > 1e-9 ||
        abs(as.numeric(dimnames(.Object@m)$year[dim(.Object@m)[2]]) -
          .Object@range["maxyear"]) > 1e-9
        ) {
            stop("range does not match supplied m and wt dimensions")
      }
      if (!missing(vcov)) .Object@vcov <- vcov
      if (!missing(units)).Object@units <- units
      .Object
    })

setValidity(
  "stk_submodels",
  function(object) {
    # the following is not throwing errors somehow...
    #if (length(object) == 0) {
    #  "There must be 3 submodels named f, n1, and sr"
    #} else
    # check dimensions of m and wt are the same as f submodel
    if (
      !identical(
        unlist(dimnames(object@m)[2:5]),
        unlist(dimnames(object@wt)[2:5])
        )) {
      "m and wt elements must share dimensions 2 to 5"
    } else {
      TRUE
    }
  }
)


#
#  accessor methods
#


#' @rdname stk_submodels-class
setMethod("m", "stk_submodels", function(object) object@m)

#' @rdname stk_submodels-class
setMethod("wt", "stk_submodels", function(object) object@wt)


#
#  show methods
#

#' @rdname stk_submodels-class
setMethod(
  "show",
  "stk_submodels",
  function(object)
  {
    cat("a4a stock model for:", object@name, "\n")
    show(submodels(object))
  }
)


#
# Other methods
#

#' @rdname a4aStkParams-class
#' @param obj the object to be subset
#' @param it iteration to be extracted
setMethod(
  "iter",
  "stk_submodels",
  function(obj, it){
    obj@vcov <- obj@vcov[,,it, drop=FALSE]
    obj@m <- iter(obj@m, it)
    obj@wt <- iter(obj@wt, it)
    obj
  }
)
