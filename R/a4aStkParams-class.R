#' @title Stock parameters class
#' @docType class
#' @name a4aStkParams
#' @rdname a4aStkParams-class
#' @template ClassDescription
#' @section Slot:
#' \describe{
#'	\item{\code{fMod}}{F submodel \code{formula}}
#'	\item{\code{n1Mod}}{first year N \code{formula}}
#'	\item{\code{srMod}}{stock-recruitment submodel \code{formula}}
#'	\item{\code{params}}{\code{FLPar} with parameters}
#'	\item{\code{vcov}}{\code{array} with variance-covariance}
#'	\item{\code{centering}}{centering values \code{numeric}}
#'	\item{\code{distr}}{statistical distribution \code{character}}
#'	\item{\code{m}}{natural mortality \code{FLQuant}}
#'	\item{\code{units}}{data units \code{character}}
#' }
#' @aliases a4aStkParams-class

setClass("a4aStkParams",
  contains = "FLComp",
  slots =
    c(
      fMod      = "formula",
      n1Mod     = "formula",
      srMod     = "formula",
      params    = "FLPar",
      vcov      = "array",
      centering = "numeric",
      distr     = "character",
      m         = "FLQuant",
      wt         = "FLQuant",
      units     = "character"
    )
)

setValidity("a4aStkParams",
  function(object) {
    # check dimensions of m and wt are the same
    if (!identical(unlist(dimnames(object@m)[2:5]),
                   unlist(dimnames(object@wt)[2:5]))) {
      "m and wt elements must share dimensions 2 to 5"
    } else {
      TRUE
    }
})

setMethod("initialize", "a4aStkParams",
    function(.Object,
              fMod      = ~ 1,
              n1Mod     = ~ 1,
              srMod     = ~ 1,
              params    = FLPar(),
              vcov      = array(),
              centering = 0,
              distr     = "lnorm",
              m         = FLQuant(),
              wt        = FLQuant(),
              units     = "NA",
              ...) {
      # initialize FLComp slots
      .Object <- callNextMethod(.Object, ...)
      # initialize remaining slots
      .Object@fMod  <- fMod
      .Object@n1Mod <- n1Mod
      .Object@srMod <- srMod
      .Object@params <- params
      .Object@vcov <- vcov
      .Object@centering <- centering
      .Object@distr <- distr
      # if missing set dimensions of of m and wt based on range
      if (missing(m) || missing(wt))
        flq <- FLQuant(
                  matrix(NA,
                     nrow = .Object@range["max"] - .Object@range["min"] + 1,
                     ncol = .Object@range["maxyear"] - .Object@range["minyear"] + 1),
                     dimnames = list(age = .Object@range["min"]:.Object@range["max"],
                                     year = .Object@range["minyear"]:.Object@range["maxyear"])
                  )
      .Object@m <- if (missing(m)) flq else m
      .Object@wt <- if (missing(wt)) flq else wt
      # throw error if range from FLComp doesn't match FLQuants 
      # (can't check this in setValidity due to callNextMethod resulting in an invalid a4aStkParams object when range is supplied)
      if (abs(as.numeric(dimnames(.Object@m)$year[1]) - .Object@range["minyear"]) > 1e-9 ||
          abs(as.numeric(dimnames(.Object@m)$year[dim(.Object@m)[2]]) - .Object@range["maxyear"]) > 1e-9) {
            stop("range does not match supplied m and wt dimensions")
      }
      .Object@units <- units
      .Object
    })






#' @rdname a4aStkParams-class
#' @template bothargs
#' @aliases a4aStkParams a4aStkParams-methods
#' @template Accessors
#' @template Constructors
setGeneric("a4aStkParams", function(object, ...) standardGeneric("a4aStkParams"))
#' @rdname a4aStkParams-class
setMethod("a4aStkParams", signature(object="missing"),
  function(...) {
    # empty
    if(missing(...)){
        new("a4aStkParams")
    # or not
    } else {
      args <- list(...)
      args$Class <- 'a4aStkParams'
      do.call("new", args)
      }
  }
)


#' @rdname a4aStkParams-class
setMethod("m", signature(object="a4aStkParams"), function(object) object @ m)

#' @rdname a4aStkParams-class
setMethod("wt", signature(object="a4aStkParams"), function(object) object @ wt)

#' @rdname a4aStkParams-class
#' @aliases fMod fMod-methods
setGeneric("fMod", function(object, ...) standardGeneric("fMod"))
#' @rdname a4aStkParams-class
setMethod("fMod", "a4aStkParams", function(object) object@fMod)

#' @rdname a4aStkParams-class
#' @param value the new object
#' @aliases fMod<- fMod<--methods
setGeneric("fMod<-", function(object,value) standardGeneric("fMod<-"))
#' @rdname a4aStkParams-class
setReplaceMethod("fMod", signature("a4aStkParams","formula"), function(object, value){
    if(all.equal(is(value), is(object@fMod))) object@fMod <- value
    object
})

#' @rdname a4aStkParams-class
#' @aliases n1Mod n1Mod-methods
setGeneric("n1Mod", function(object, ...) standardGeneric("n1Mod"))
#' @rdname a4aStkParams-class
setMethod("n1Mod", "a4aStkParams", function(object) object@n1Mod)

#' @rdname a4aStkParams-class
#' @aliases n1Mod<- n1Mod<--methods
setGeneric("n1Mod<-", function(object,value) standardGeneric("n1Mod<-"))
#' @rdname a4aStkParams-class
setReplaceMethod("n1Mod", signature("a4aStkParams","formula"), function(object, value){
    if(all.equal(is(value), is(object@n1Mod))) object@n1Mod <- value
    object
})

#' @rdname a4aStkParams-class
#' @aliases srMod rMod-methods
setGeneric("srMod", function(object, ...) standardGeneric("srMod"))
#' @rdname a4aStkParams-class
setMethod("srMod", "a4aStkParams", function(object) object@srMod)

#' @rdname a4aStkParams-class
#' @aliases srMod<- srMod<--methods
setGeneric("srMod<-", function(object,value) standardGeneric("srMod<-"))
#' @rdname a4aStkParams-class
setReplaceMethod("srMod", signature("a4aStkParams","formula"), function(object, value){
    if(all.equal(is(value), is(object@srMod))) object@srMod <- value
    object
})

#' @rdname a4aStkParams-class
setMethod("params", "a4aStkParams", function(object) object@params)

#' @rdname a4aStkParams-class
setReplaceMethod("params", signature("a4aStkParams","FLPar"), function(object, value){
    if(all.equal(is(value), is(object@params))) object@params <- value
    object
})

#' @rdname a4aStkParams-class
setMethod("distr", "a4aStkParams", function(object) object@distr)

#' @rdname a4aStkParams-class
setReplaceMethod("distr", signature("a4aStkParams","character"), function(object, value){
    if(all.equal(is(value), is(object@distr))) object@distr <- value
    object
})

#' @rdname a4aStkParams-class
setMethod("vcov", "a4aStkParams", function(object) object@vcov)

#' @rdname a4aStkParams-class
setReplaceMethod("vcov", signature("a4aStkParams","array"), function(object, value){
    if(all.equal(is(value), is(object@vcov))) object@vcov <- value
    object
})

#' @rdname a4aStkParams-class
#' @param obj the object to be subset
#' @param it iteration to be extracted
setMethod("iter", "a4aStkParams", function(obj, it){
	obj@vcov <- obj@vcov[,,it, drop=FALSE]
	obj@params <- iter(obj@params, it)
	obj@m <- iter(obj@m, it)
	obj@wt <- iter(obj@wt, it)
	obj@centering <- obj@centering[it]
	obj
})































