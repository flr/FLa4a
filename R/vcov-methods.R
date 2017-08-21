#==================================================================== 
#    vcov  methods
#==================================================================== 

#' @title Variance-covariance matrix
#' @name vcov
#' @aliases vcov,a4aFitSA-method
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
    nmodels <- length(object)
    corblocks <- corBlocks(object)
    
    cormat <-
      do.call(rbind, 
        lapply(1:nmodels, function(i) {
          do.call(cbind,
            lapply(1:nmodels, function(j) {
              if (i == j) {
                # return correlation matrix
                cov2cor(vcov(object[[i]]))
              } else {
                out <- corblocks[[paste(names(object)[sort(c(i, j))], collapse = ".")]]
                if (i > j) {
                  t(out)
                } else {
                  out
                }
              }
            }))
          })
      )
    vardiag <- unlist(lapply(object, function(x) diag(vcov(x))))
    V <- diag(sqrt(vardiag)) %*% cormat %*% diag(sqrt(vardiag))
    # add names in?
    npar <- sapply(object, function(x) length(coef(x)))
    parnames <- lapply(object, function(x) dimnames(coef(x))$params)
    colnames(V) <- rownames(V) <- paste(rep(names(object), npar), unlist(parnames), sep = ":")
    V
  })


#' @rdname vcov-methods
setMethod("vcov", signature(object = "submodel"),
  function(object) {
      object@vcov
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
setMethod("vcov<-", signature(object = "submodel", value = "numeric"),
  function(object, ..., value) {
      object @ vcov[] <- value
      object
  })

#' @rdname vcov-methods
setMethod("vcov<-", signature(object = "submodel", value = "matrix"),
  function(object, ..., value) {
      object @ vcov[] <- value
      object
  })



