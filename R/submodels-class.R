vsmods <- function(object){
	
  # All items are FLStock
  if(!all(unlist(lapply(object, is, 'submodel'))))
      return("Components must be submodel")	
	
	return(TRUE)
}

# submodel-class - «Short one line description»
# submodel-class


#' a4aFitSA extends \code{"a4aFit"} class.
#'
#' Some details about this class and my plans for it in the body.
#'
#' \describe{
#'    \item{myslot1}{A logical keeping track of something.}
#'
#'    \item{myslot2}{An integer specifying something else.}
#' 
#'    \item{myslot3}{A data.frame holding some data.}
#'  }
#' @name submodels-class
#' @rdname submodels-class
#' @exportClass submodels
setClass("submodels", contains="FLComps",
	validity=vsmods
)


# constructor

#' @export
setGeneric("submodels", function(object, ...)
	standardGeneric("submodels"))

#' @export
setMethod("submodels", signature(object="submodel"), function(object, ...) {
    lst <- c(object, list(...))
    submodels(lst)
})

#' @export
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

#' @export
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

