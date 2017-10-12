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
    if (nmodels == 1) {
      return(vcov(object[[1]]))
    }
    corblocks <- corBlocks(object)
    niters <- unique(sapply(object, function(x) dim(vcov(x))[3]))
    V.lst <- lapply(1:niters,
                    function(it) {
                      cormat <-
                        do.call(rbind,
                                lapply(1:nmodels, function(i) {
                                  do.call(cbind,
                                          lapply(1:nmodels, function(j) {
                                            if (i == j) {
                                              # return correlation matrix
                                              out <- vcov(object[[i]])
                                              cov2cor(dropMatrixIter(out, iter = it))
                                            } else {
                                              out <- corblocks[[paste(names(object)[sort(c(i, j))], collapse = ".")]]
                                              out <- dropMatrixIter(out, it)
                                              if (i > j) {
                                                t(out)
                                              } else {
                                                out
                                              }
                                            }
                                          })
                                   )
                                })
                        )
                      vardiag <- unlist(lapply(object, function(x) diag(as.matrix(vcov(x)[,,it]))))
                      diag(sqrt(vardiag)) %*% cormat %*% diag(sqrt(vardiag))
                    })
    npar <- sapply(object, function(x) dim(coef(x))[1])
    V <- array(unlist(V.lst), dim = c(sum(npar), sum(npar), niters))
    # add names in?
    parnames <- lapply(object, function(x) dimnames(coef(x))$params)
    dimnames(V) <- list(paste(rep(names(object), npar), unlist(parnames), sep = ":"),
                        paste(rep(names(object), npar), unlist(parnames), sep = ":"),
                        1:niters)
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

#' @rdname vcov-methods
setMethod("vcov<-", signature(object = "submodel", value = "array"),
  function(object, ..., value) {
      object @ vcov[] <- value
      object
  })


