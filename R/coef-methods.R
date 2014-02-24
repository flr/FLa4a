#==================================================================== 
#    coef  methods
#==================================================================== 

#' Methods to extract and replace the model coefficients
#' @name coef
#' @docType methods
#' @rdname coef-methods

setGeneric("coef", function(object, ...) standardGeneric("coef"))

#' @rdname coef-methods
#' @aliases coef,a4aFitSA-method
setMethod("coef", signature(object = "a4aFitSA"),
  function(object) {
	  coef(pars(object))
  })


#' @rdname coef-methods
#' @aliases coef,SCAPars-method
setMethod("coef", signature(object = "SCAPars"),
  function(object) {
    list(
      stkmodel = coef(stkmodel(object)),
      qmodel   = coef(qmodel(object)),
      vmodel   = coef(vmodel(object))
    )
  })

#' @rdname coef-methods
#' @aliases coef,a4aStkParams-method
setMethod("coef", signature(object = "a4aStkParams"),
  function(object) {
      object @ params
  })

#' @rdname coef-methods
#' @aliases coef,submodels-method
setMethod("coef", signature(object = "submodels"),
  function(object) {
      lapply(object, coef)
  })


#' @rdname coef-methods
#' @aliases coef,submodel-method
setMethod("coef", signature(object = "submodel"),
  function(object) {
      object @ params
  })

#==================================================================== 
#    coef<-  methods
#==================================================================== 

setGeneric("coef<-", function(object, ..., value) standardGeneric("coef<-"))

#' @rdname coef-methods
#' @aliases coef<-,a4aFitSA,numeric-method
setMethod("coef<-", signature(object = "a4aFitSA", value = "numeric"),
  function(object, ..., value) {
    coef(object @ pars) <- value
    object
  })


#' @rdname coef-methods
#' @aliases coef<-,SCAPars,numeric-method
setMethod("coef<-", signature(object = "SCAPars", value = "numeric"),
  function(object, ..., value) {
    v <- coef(object)
    old <- unlist(v)
    new <- rep_len(unlist(value), length = length(old))
    
    coef(object @ stkmodel) <- new[grep("stkmodel", names(old))]
    coef(object @ qmodel) <- new[grep("qmodel.", names(old))]
    coef(object @ vmodel) <- new[grep("vmodel.", names(old))]

    object
  })


#' @rdname coef-methods
#' @aliases coef<-,a4aStkParams,numeric-method
setMethod("coef<-", signature(object = "a4aStkParams", value = "numeric"),
  function(object, ..., value) {    
    object @ params[] <- value
    object
  })

#' @rdname coef-methods
#' @aliases coef<-,submodels,numeric-method
setMethod("coef<-", signature(object = "submodels", value = "numeric"),
  function(object, ..., value) {
    v <- coef(object)
    old <- unlist(v)
    new <- rep_len(unlist(value), length = length(old))
    
    for (i in seq_along(object)) {
      object[[i]] @ params[] <- new[grep(object[[i]] @ name, names(old))]  
    }
    object
  })


#' @rdname coef-methods
#' @aliases coef<-,submodel,numeric-method
setMethod("coef<-", signature(object = "submodel", value = "numeric"),
  function(object, ..., value) {
      object @ params[] <- value
      object
  })

