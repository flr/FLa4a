
#' @title Submodels class
#' @docType class
#' @name submodels
#' @rdname submodels-class
#' @template ClassDescription
#' @note This class is similar to other 'plural' calsses in \code{FLR}. It is a list constrained to having all elements of the same class, in this case \code{submodel}. Otherwise it works exacly as any other list.
#' @aliases submodels-class
submodels <-
  setClass("submodels",
    contains = "FLComps",
    slots = c("corBlocks" = "list"))

#' @rdname submodels-class
#' @template Constructors
#' @template bothargs
#' @aliases submodels submodels-methods
setGeneric("submodels")

setMethod("initialize", "submodels",
  function(.Object,
           ...,
           corBlocks,
           names) {
      .Object <- callNextMethod(.Object, ...)
      if (!missing(names)) {
        # need to apply new() recursively to maintain a valid object
        asListObject <- as(.Object, "list")
        names(asListObject) <- names
        for (i in seq_along(asListObject)) name(asListObject[[i]]) <- names[i]
        .Object <- new("submodels", asListObject, corBlocks = corBlocks)
      }
      # this is needed to avoid attempted evaluation of names argument
      # when calling the names function in following if statment
      names <- ""
      if (any(is.na(names(.Object)) | names(.Object) == "")) {
        names(.Object) <- unname(sapply(.Object, name))
      }
      # finally check for corrupt submodels and apply a simple naming scheme
      if (any(names(.Object) == "")) {
        names <- names(.Object)
        names[names == ""] <- letters[1:sum(names == "")]
        .Object <- new("submodels", as(.Object, "list"), corBlocks = corBlocks, names = make.unique(names))
      }
      if (!missing(corBlocks)) {
        .Object@corBlocks <- corBlocks
      } else {
        # generate from submodel dimensions
        nmodels <- length(.Object)
        # catch case of empty and single objects
        if (nmodels > 1) {
          npar <- sapply(.Object, function(x) length(coef(x)))
          parnames <- lapply(.Object, function(x) dimnames(coef(x))$params)
          modelpairs <- combn(seq(nmodels), 2)
          .Object@corBlocks <-
            lapply(seq(ncol(modelpairs)),
                   function(i)
                     array(0,
                            dim = c(npar[modelpairs[1,i]],
                                    npar[modelpairs[2,i]],
                                    iter = 1),
                            dimnames = c(parnames[modelpairs[,i]], "1")))
          names(.Object@corBlocks) <- apply(modelpairs, 2, function(x) paste(names(.Object)[x], collapse = "."))
        }
      }
      .Object
})


setValidity("submodels",
  function(object) {
    # All items are submodel-class
    if(!all(sapply(object, is, 'submodel'))) {
      "Components must be submodel"
    } else {
      TRUE
    }
})



#
#  accessor methods
#

#' @rdname submodels-class
#' @aliases corBlocks corBlock-methods
setGeneric("corBlocks", function(object, ...) standardGeneric("corBlocks"))
#' @rdname submodels-class
setMethod("corBlocks", "submodels", function(object) object@corBlocks)


#' @rdname submodels-class
setMethod("params", "submodels", function(object) lapply(object, coef))

#' @rdname submodels-class
setMethod("sMod", "submodels", function(object) lapply(object, sMod))

#' @rdname submodels-class
setMethod("vcov", "submodels", function(object) lapply(object, vcov))

#' @rdname submodels-class
setMethod("formula", "submodels", function(x) lapply(x, formula))


#
#  assignment methods
#

#' @rdname submodels-class
#' @param value value the new object
#' @aliases corBlocks<-
setGeneric("corBlocks<-", function(object, ..., value) standardGeneric("corBlocks<-"))

#' @rdname vcov-methods
#' @param value value the new object
setMethod("corBlocks<-", signature(object = "submodels", value = "list"),
  function(object, ..., value) {
    object@corBlocks[] <- value
    object
  })

# method.skeleton("$<-", signature(object = "submodels", value = "submodel"),  file = stdout())

setMethod("$<-",
  signature(x = "submodels", value = "submodel"),
  function(x, name, value) {
    x[[name]] <- value
    x
  })

setMethod("[[<-",
  c("submodels", "character", "missing"),
  function (x, i, j, ..., value)
  {
    lst <- as(x, "list")
    names(lst) <- names(x)
    lst[[i]] <- value
    new("submodels", lst, corBlocks = x@corBlocks)
  }
)

setMethod("[[<-",
  c("submodels", "numeric", "missing"),
  function (x, i, j, ..., value)
  {
    lst <- as(x, "list")
    names(lst) <- names(x)
    lst[[i]] <- value
    new("submodels", lst, corBlocks = x@corBlocks)
  }
)


#
#  show methods
#

setMethod("show", "submodels",
  function(object)
  {
    cat("submodels:\n")
    if (length(object) == 0) {
      cat("empty object\n")
    } else {
      fmt <- paste0("\t %", max(nchar(sapply(object, name))), "smodel: ")
      for (i in object) {
        cat(sprintf(fmt, name(i))); print(formula(i), showEnv = FALSE)
      }
    }
 })



#
# Coersion methods
#

# method.skeleton("coerce", "submodels",  file = stdout())

setMethod("coerce",
  signature(from = "submodels", to = "submodel"),
  function (from, to, strict = TRUE)
  {
    stop("No method described yet.")
  }
)

#
# Other methods
#


# method.skeleton("propagate", "submodels",  file = stdout())
#' @rdname submodels-class
setMethod("propagate",
  signature(object = "submodels"),
  function (object, iter, fill.iter = TRUE)
  {
    if (length(object) == 0) {
      stop("propagate can only extend a submodels object that has at least one submodel")
    }

    lst <- as(object, "list")
    lst <- lapply(lst, propagate, iter = iter, fill.iter = fill.iter)

    # propagate corblocks
    corBlocks <-
      lapply(object@corBlocks, function(x)
      {
        dob <- dim(x)

        if (iter != dob[3]) {
          # CHECK no iters in object
          if (dob[3] > 1) stop("propagate can only extend objects with no iters")

          out <- array(NA, dim = c(dob[1:2], iter), dimnames = c(dimnames(x)[1:2], list(1:iter)))
          if (fill.iter) {
            out[] <- as.vector(x)
          } else {
            out[,,1] <- as.vector(x)
          }
        } else {
          out <- x
        }
        out
      })

    new("submodels", lst, corBlocks = corBlocks)
  }
)


#' @rdname submodels-class
#' @param obj the object to be subset
#' @param it iteration to be extracted
setMethod("iter", "submodels", function(obj, it){
  out <- submodels(lapply(obj, iter, it))
  names(out) <- names(obj)
  out
})


