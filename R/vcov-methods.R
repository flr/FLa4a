#==================================================================== 
#    vcov  methods
#==================================================================== 

#' Methods to extract and replace the variance-covariance matrix
#' @name vcov
#' @rdname vcov-methods
#' @aliases vcov,a4aFitSA-method
setMethod("vcov", signature(object = "a4aFitSA"),
  function(object) {
    vcov(pars(object))
  })


#' @rdname vcov-methods
#' @aliases vcov,SCAPars-method
setMethod("vcov", signature(object = "SCAPars"),
  function(object) {
    list(
      stkmodel = vcov(stkmodel(object)),
      qmodel   = vcov(qmodel(object)),
      vmodel   = vcov(vmodel(object))
    )
  })

#' @rdname vcov-methods
#' @aliases vcov,a4aStkParams-method
setMethod("vcov", signature(object = "a4aStkParams"),
  function(object) {
      object @ vcov
  })

#' @rdname vcov-methods
#' @aliases vcov,submodels-method
setMethod("vcov", signature(object = "submodels"),
  function(object) {
      lapply(object, vcov)
  })


#' @rdname vcov-methods
#' @aliases coef,submodel-method
setMethod("vcov", signature(object = "submodel"),
  function(object) {
      object @ vcov
  })

#==================================================================== 
#    vcov<-  methods
#==================================================================== 

#' @rdname vcov-methods
#' @aliases vcov<-,a4aFitSA,numeric-method
setMethod("vcov<-", signature(object = "a4aFitSA", value = "numeric"),
  function(object, ..., value) {
    vcov(object @ pars) <- value
    object
  })


#' @rdname vcov-methods
#' @aliases vcov<-,SCAPars,numeric-method
setMethod("vcov<-", signature(object = "SCAPars", value = "numeric"),
  function(object, ..., value) {
    v <- vcov(object)
    old <- unlist(v)
    new <- rep_len(unlist(value), length = length(old))
    
    vcov(object @ stkmodel) <- new[grep("stkmodel", names(old))]
    vcov(object @ qmodel) <- new[grep("qmodel.", names(old))]
    vcov(object @ vmodel) <- new[grep("vmodel.", names(old))]

    object
  })


#' @rdname vcov-methods
#' @aliases vcov<-,a4aStkParams,numeric-method
setMethod("vcov<-", signature(object = "a4aStkParams", value = "numeric"),
  function(object, ..., value) {    
    object @ vcov[] <- value
    object
  })

#' @rdname vcov-methods
#' @aliases vcov<-,submodels,numeric-method
setMethod("vcov<-", signature(object = "submodels", value = "numeric"),
  function(object, ..., value) {
    v <- vcov(object)
    old <- unlist(v)
    new <- rep_len(unlist(value), length = length(old))
    
    for (i in seq_along(object)) {
      object[[i]] @ vcov[] <- new[grep(object[[i]] @ name, names(old))]  
    }
    object
  })

#' @rdname vcov-methods
#' @aliases vcov<-,submodel,numeric-method
setMethod("vcov<-", signature(object = "submodel", value = "numeric"),
  function(object, ..., value) {
      object @ vcov[] <- value
      object
  })

