vsmods <- function(object){
	
  # All items are FLStock
  if(!all(unlist(lapply(object, is, 'submodel'))))
      return("Components must be submodel")	
	
	return(TRUE)
}

#' @title Submodels class
#' @docType class
#' @name submodel
#' @rdname submodels-class
#' @template ClassDescription
#' @note This class is similar to other 'plural' calsses in \code{FLR}. It's a list constraint by having all elements of the same class, in this case \code{submodel}. Otherwise it works exacly as any other list.
#' @aliases submodels-class
setClass("submodels", contains="FLComps",
	validity=vsmods
)

#' @rdname submodels-class
#' @template Constructors
#' @aliases submodels submodels-methods submodels,missing-method

setGeneric("submodels", function(object, ...)
	standardGeneric("submodels"))

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
#' @aliases submodels,submodel-method
setMethod("submodels", signature(object="submodel"), function(object, ...) {
    lst <- c(object, list(...))
    submodels(lst)
})

#' @rdname submodels-class
#' @aliases submodels,list-method
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

