
#' @title Submodels class
#' @docType class
#' @name submodels
#' @rdname submodels-class
#' @template ClassDescription
#' @note This class is similar to other 'plural' calsses in \code{FLR}. It is a list constrained to having all elements of the same class, in this case \code{submodel}. Otherwise it works exacly as any other list.
#' @aliases submodels-class
submodels <- setClass("submodels", contains = "FLComps")

#' @rdname submodels-class
#' @template Constructors
#' @template bothargs
#' @aliases submodels submodels-methods
setGeneric("submodels")

#' @rdname submodels-class
setMethod("submodels", "missing",
  function(...) {
    # empty
    if(missing(...)){
      new("submodels")
    # or not
    } else {
      args <- list(...)
      object <- args[!names(args)%in%c('names', 'desc', 'lock')]
      args <- args[!names(args)%in%names(object)]
      do.call('submodels',  c(list(object=object), args))
    }
  }
)


#' @rdname submodels-class
setMethod("submodels", "list",
  function(object, ...) {
    
    args <- list(...)
    
    # names in args, ... 
    if("names" %in% names(args)) {
      names <- args[['names']]
    } else {
    # ... or in object,
      if(!is.null(names(object))) {
        names <- names(object)
    # ... or in elements, ...
      } else {
        names <- unlist(lapply(object, name))
        # ... or 1:n
        idx <- names == "NA" | names == ""
        if(any(idx))
          names[idx] <- as.character(length(names))[idx]
      }
    }

    # desc & lock
    args <- c(list(Class="submodels", .Data=object, names=names),
      args[!names(args) %in% 'names'])

    return(
      do.call('new', args)
      )
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
setMethod("params", "submodels", function(object) lapply(object, params))

#' @rdname submodels-class
setMethod("sMod", "submodels", function(object) lapply(object, sMod))

#' @rdname submodels-class
setMethod("vcov", "submodels", function(object) lapply(object, vcov))

#' @rdname submodels-class
setMethod("formula", "submodels", function(x) lapply(x, formula))



#
#  show methods
#

setMethod("show", "submodels",
  function(object)
  {
    cat("submodels:\n")
    fmt <- paste0("\t %", max(nchar(sapply(object, name))), "smodel: ")
    for (i in object) {
      cat(sprintf(fmt, name(i))); print(formula(i), showEnv = FALSE)
    }
 })


#
# Other methods
#

#' @rdname submodels-class
#' @param obj the object to be subset
#' @param it iteration to be extracted 
setMethod("iter", "submodels", function(obj, it){
  out <- submodels(lapply(obj, iter, it))
  names(out) <- names(obj)
  out
})