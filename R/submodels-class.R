
#' @title Submodels class
#' @docType class
#' @name submodels
#' @rdname submodels-class
#' @template ClassDescription
#' @note This class is similar to other 'plural' calsses in \code{FLR}. It is a list constrained to having all elements of the same class, in this case \code{submodel}. Otherwise it works exacly as any other list.
#' @aliases submodels-class
setClass("submodels", contains = "FLComps")

setValidity("submodels", 
  function(object) {
    # All items are submodel-class
    if(!all(sapply(object, is, 'submodel'))) {
      "Components must be submodel"
    } else {
      TRUE
    }
})




#' @rdname submodels-class
#' @template Constructors
#' @template bothargs
#' @aliases submodels submodels-methods
setGeneric("submodels", function(object, ...)
	standardGeneric("submodels"))
#' @rdname submodels-class
setMethod("submodels", signature(object="missing"),
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
setMethod("submodels", signature(object="submodel"), function(object, ...) {
    lst <- c(object, list(...))
    submodels(lst)
})

#' @rdname submodels-class
setMethod("submodels", signature(object="list"),
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
      args[!names(args)%in%'names'])

    return(
      do.call('new', args)
      )

}) # }}}

#' @rdname submodels-class
#' @param obj the object to be subset
#' @param it iteration to be extracted 
setMethod("iter", "submodels", function(obj, it){
	out <- submodels(lapply(obj, iter, it))
	names(out) <- names(obj)
	out
})

