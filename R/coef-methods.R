#==================================================================== 
#    coef  methods
#==================================================================== 

#' @title coefficients extract and replacement
#' @name coef
#' @docType methods
#' @rdname coef-methods
#' @template bothargs
#' @aliases coef,a4aFitSA-methods
#' @description Methods to extract and replace the model coefficients.

setGeneric("coef", function(object, ...) useAsDefault=stats::coef)

#' @rdname coef-methods
setMethod("coef", signature(object = "a4aFitSA"),
  function(object) {
	  coef(pars(object))
  })


#' @rdname coef-methods
setMethod("coef", signature(object = "SCAPars"),
  function(object) {
    list(
      stkmodel = coef(stkmodel(object)),
      qmodel   = coef(qmodel(object)),
      vmodel   = coef(vmodel(object))
    )
  })

#' @rdname coef-methods
setMethod("coef", signature(object = "a4aStkParams"),
  function(object) {
      object @ params
  })

#' @rdname coef-methods
setMethod("coef", signature(object = "submodels"),
  function(object) {
      lapply(object, coef)
  })


#' @rdname coef-methods
setMethod("coef", signature(object = "submodel"),
  function(object) {
      object @ params
  })

#==================================================================== 
#    coef<-  methods
#==================================================================== 

#' @rdname coef-methods
#' @param value the new object
#' @aliases coef<-,a4aFitSA-methods
setGeneric("coef<-", function(object, ..., value) standardGeneric("coef<-"))

#' @rdname coef-methods
setMethod("coef<-", signature(object = "a4aFitSA", value = "numeric"),
  function(object, ..., value) {
    coef(object @ pars) <- value
    object
  })

#' @rdname coef-methods
setMethod("coef<-", signature(object = "SCAPars", value = "numeric"),
  function(object, ..., value) {
    v <- coef(object)
    old <- unlist(v)
    new <- rep_len(unlist(value), length.out = length(old))
    
    coef(object @ stkmodel) <- new[grep("stkmodel", names(old))]
    coef(object @ qmodel) <- new[grep("qmodel.", names(old))]
    coef(object @ vmodel) <- new[grep("vmodel.", names(old))]

    object
  })


#' @rdname coef-methods
setMethod("coef<-", signature(object = "a4aStkParams", value = "numeric"),
  function(object, ..., value) {    
    object @ params[] <- value
    object
  })

#' @rdname coef-methods
setMethod("coef<-", signature(object = "submodels", value = "numeric"),
  function(object, ..., value) {
    v <- coef(object)
    old <- unlist(v)
    new <- rep_len(unlist(value), length.out = length(old))
    
    for (i in seq_along(object)) {
      object[[i]] @ params[] <- new[grep(object[[i]] @ name, names(old))]  
    }
    object
  })


#' @rdname coef-methods
setMethod("coef<-", signature(object = "submodel", value = "numeric"),
  function(object, ..., value) {
      object @ params[] <- value
      object
  })

