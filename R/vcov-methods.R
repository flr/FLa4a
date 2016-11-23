#==================================================================== 
#    vcov  methods
#==================================================================== 

#' @title Variance-covariance matrix
#' @description Methods to extract and replace the variance-covariance matrix.
#' @template bothargs
#' @rdname vcov-methods
setMethod("vcov", signature(object = "a4aFitSA"),
  function(object) {
    vcov(pars(object))
  })


#' @rdname vcov-methods
setMethod("vcov", signature(object = "SCAPars"),
  function(object) {
    list(
      stkmodel = vcov(stkmodel(object)),
      qmodel   = vcov(qmodel(object)),
      vmodel   = vcov(vmodel(object))
    )
  })

#' @rdname vcov-methods
setMethod("vcov", signature(object = "submodels"),
  function(object) {
      lapply(object, vcov)
  })


#' @rdname vcov-methods
setMethod("vcov", signature(object = "submodel"),
  function(object) {
      object @ vcov
  })

#==================================================================== 
#    vcov<-  methods
#==================================================================== 

#' @rdname vcov-methods
#' @param value the new object
setMethod("vcov<-", signature(object = "a4aFitSA", value = "numeric"),
  function(object, ..., value) {
    vcov(object @ pars) <- value
    object
  })


#' @rdname vcov-methods
setMethod("vcov<-", signature(object = "SCAPars", value = "numeric"),
  function(object, ..., value) {
    v <- vcov(object)
    old <- unlist(v)
    new <- rep_len(unlist(value), length.out = length(old))
    
    vcov(object @ stkmodel) <- new[grep("stkmodel", names(old))]
    vcov(object @ qmodel) <- new[grep("qmodel.", names(old))]
    vcov(object @ vmodel) <- new[grep("vmodel.", names(old))]

    object
  })


#' @rdname vcov-methods
setMethod("vcov<-", signature(object = "a4aStkParams", value = "numeric"),
  function(object, ..., value) {    
    object @ vcov[] <- value
    object
  })

#' @rdname vcov-methods
setMethod("vcov<-", signature(object = "submodels", value = "numeric"),
  function(object, ..., value) {
    v <- vcov(object)
    old <- unlist(v)
    new <- rep_len(unlist(value), length.out = length(old))
    
    for (i in seq_along(object)) {
      object[[i]] @ vcov[] <- new[grep(object[[i]] @ name, names(old))]  
    }
    object
  })

#' @rdname vcov-methods
setMethod("vcov<-", signature(object = "submodel", value = "numeric"),
  function(object, ..., value) {
      object @ vcov[] <- value
      object
  })

