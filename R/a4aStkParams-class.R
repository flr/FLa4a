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
  representation = 
    representation(
      "FLComp",
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
    ),
  prototype = 
    prototype(
      name      = character(0),
      desc      = character(0),
      range     = c(min=0, max=0, plusgroup=0, minyear=0, maxyear=0),
      fMod      = ~1,
      n1Mod     = ~1,
      srMod     = ~1,
      params    = FLPar(),
      vcov      = array(),
      centering = 0,
      distr     = "lnorm",
      m         = FLQuant(),
      wt         = FLQuant(),
      units     = "NA"
    )
)

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
#' @aliases rMod rMod-methods
setGeneric("rMod", function(object, ...) standardGeneric("rMod"))
#' @rdname a4aStkParams-class
setMethod("rMod", "a4aStkParams", function(object) object@rMod)

#' @rdname a4aStkParams-class
#' @aliases rMod<- rMod<--methods
setGeneric("rMod<-", function(object,value) standardGeneric("rMod<-"))
#' @rdname a4aStkParams-class
setReplaceMethod("rMod", signature("a4aStkParams","formula"), function(object, value){
    if(all.equal(is(value), is(object@rMod))) object@rMod <- value
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































