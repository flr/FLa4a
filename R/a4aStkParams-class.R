#' @title Stock parameters class
#' @docType class
#' @name a4aStkParams
#' @rdname a4aStkParams-class
#' @template ClassDescription
#' @section Slot: 
#' \describe{
#'	\item{\code{fMod}}{F submodel \code{formula}}
#'	\item{\code{n1Mod}}{First year N \code{formula}}
#'	\item{\code{srMod}}{Stock-recruitment submodel \code{formula}}
#'	\item{\code{params}}{\code{FLPar} with parameters}
#'	\item{\code{vcov}}{\code{array} with variance-covariance}
#'	\item{\code{centering}}{Centering values \code{numeric}}
#'	\item{\code{distr}}{Statistical distribution \code{character}}
#'	\item{\code{m}}{Natural mortality \code{FLQuant}}
#'	\item{\code{units}}{Data units \code{character}}
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
      units     = "NA"
    )
)

#' @rdname a4aStkParams-class
#' @aliases a4aStkParams a4aStkParams-methods a4aStkParams,missing-method
#' @template Accessors
#' @template Constructors
setGeneric("a4aStkParams", function(object, ...) standardGeneric("a4aStkParams"))
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


#' @rdname a4aGr-class
#' @aliases m,a4aStkParams-method
setMethod("m", signature(object="a4aStkParams"), function(object) object @ m)

#' @rdname a4aStkParams-class
#' @aliases fMod fMod-methods fMod,a4aStkParams-method
setGeneric("fMod", function(object, ...) standardGeneric("fMod"))
setMethod("fMod", "a4aStkParams", function(object) object@fMod)

#' @rdname a4aStkParams-class
#' @aliases fMod<- fMod<--methods fMod<-,a4aStkParams,formula-method
setGeneric("fMod<-", function(object,value) standardGeneric("fMod<-"))
setReplaceMethod("fMod", signature("a4aStkParams","formula"), function(object, value){
    if(all.equal(is(value), is(object@fMod))) object@fMod <- value
    object
})

#' @rdname a4aStkParams-class
#' @aliases n1Mod n1Mod-methods n1Mod,a4aStkParams-method
setGeneric("n1Mod", function(object, ...) standardGeneric("n1Mod"))
setMethod("n1Mod", "a4aStkParams", function(object) object@n1Mod)

#' @rdname a4aStkParams-class
#' @aliases n1Mod<- n1Mod<--methods n1Mod<-,a4aStkParams,formula-method
setGeneric("n1Mod<-", function(object,value) standardGeneric("n1Mod<-"))
setReplaceMethod("n1Mod", signature("a4aStkParams","formula"), function(object, value){
    if(all.equal(is(value), is(object@n1Mod))) object@n1Mod <- value
    object
})

#' @rdname a4aStkParams-class
#' @aliases rMod rMod-methods rMod,a4aStkParams-method
setGeneric("rMod", function(object, ...) standardGeneric("rMod"))
setMethod("rMod", "a4aStkParams", function(object) object@rMod)

#' @rdname a4aStkParams-class
#' @aliases rMod<- rMod<--methods rMod<-,a4aStkParams,formula-method
setGeneric("rMod<-", function(object,value) standardGeneric("rMod<-"))
setReplaceMethod("rMod", signature("a4aStkParams","formula"), function(object, value){
    if(all.equal(is(value), is(object@rMod))) object@rMod <- value
    object
})

#' @rdname a4aStkParams-class
#' @aliases params,a4aStkParams-method
setMethod("params", "a4aStkParams", function(object) object@params)

#' @name a4aStkParams
#' @rdname a4aStkParams-class
#' @aliases params<-,a4aStkParams,FLPar-method
setReplaceMethod("params", signature("a4aStkParams","FLPar"), function(object, value){
    if(all.equal(is(value), is(object@params))) object@params <- value
    object
})

#' @rdname a4aStkParams-class
#' @aliases distr distr-methods distr,a4aStkParams-method
setGeneric("distr", function(object, ...) standardGeneric("distr"))
setMethod("distr", "a4aStkParams", function(object) object@distr)

#' @rdname a4aStkParams-class
#' @aliases distr<- distr<--methods distr<-,a4aStkParams,character-method
setGeneric("distr<-", function(object, value) standardGeneric("distr<-"))
setReplaceMethod("distr", signature("a4aStkParams","character"), function(object, value){
    if(all.equal(is(value), is(object@distr))) object@distr <- value
    object
})

#' @rdname a4aStkParams-class
#' @aliases vcov,a4aStkParams-method
setMethod("vcov", "a4aStkParams", function(object) object@vcov)

#' @name vcov<-
#' @rdname a4aStkParams-class
#' @aliases vcov<-,a4aStkParams-method
setReplaceMethod("vcov", signature("a4aStkParams","array"), function(object, value){
    if(all.equal(is(value), is(object@vcov))) object@vcov <- value
    object
})

