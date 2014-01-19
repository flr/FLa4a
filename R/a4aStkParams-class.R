
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


setMethod("m", signature(object="a4aStkParams"),
  function(object) object @ m)


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


# accessors
setGeneric("fMod", function(object, ...) standardGeneric("fMod"))
setMethod("fMod", "a4aStkParams", function(object) object@fMod)

setGeneric("fMod<-", function(object,value) standardGeneric("fMod<-"))
setReplaceMethod("fMod", "a4aStkParams", function(object, value){
    if(all.equal(is(value), is(object@fMod))) object@fMod <- value
    object
})

setGeneric("n1Mod", function(object, ...) standardGeneric("n1Mod"))
setMethod("n1Mod", "a4aStkParams", function(object) object@n1Mod)

setGeneric("n1Mod<-", function(object,value) standardGeneric("n1Mod<-"))
setReplaceMethod("n1Mod", "a4aStkParams", function(object, value){
    if(all.equal(is(value), is(object@n1Mod))) object@n1Mod <- value
    object
})

setGeneric("rMod", function(object, ...) standardGeneric("rMod"))
setMethod("rMod", "a4aStkParams", function(object) object@rMod)

setGeneric("rMod<-", function(object,value) standardGeneric("rMod<-"))
setReplaceMethod("rMod", "a4aStkParams", function(object, value){
    if(all.equal(is(value), is(object@rMod))) object@rMod <- value
    object
})


setMethod("params", "a4aStkParams", function(object) object@params)

setGeneric("params<-", function(object, value) standardGeneric("params<-"))
setReplaceMethod("params", "a4aStkParams", function(object, value){
    if(all.equal(is(value), is(object@params))) object@params <- value
    object
})

setGeneric("distr", function(object, ...) standardGeneric("distr"))
setMethod("distr", "a4aStkParams", function(object) object@distr)

setGeneric("distr<-", function(object, value) standardGeneric("distr<-"))
setReplaceMethod("distr", "a4aStkParams", function(object, value){
    if(all.equal(is(value), is(object@distr))) object@distr <- value
    object
})

setMethod("vcov", "a4aStkParams", function(object) object@vcov)

setGeneric("vcov<-", function(object, value) standardGeneric("vcov<-"))
setReplaceMethod("vcov", "a4aStkParams", function(object, value){
    if(all.equal(is(value), is(object@vcov))) object@vcov <- value
    object
})

