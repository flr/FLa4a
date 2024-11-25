#====================================================================
# plural class for a4aFitSA (used for model averaging)
#====================================================================

#' @rdname a4aFitSA-class
#' @aliases a4aFitSAs-class

setClass("a4aFitSAs",
  contains="a4aFits"
)

setValidity("a4aFitSAs",
  function(object) {
    if(!all(sapply(object, is, 'a4aFitSA'))) {
      "Components must be a4aFitSA"
    } else {
      TRUE
    }
})

#' @rdname a4aFitSA-class
#' @aliases a4aFitSAs a4aFitSAs-methods
setGeneric("a4aFitSAs", function(object, ...) standardGeneric("a4aFitSAs"))

#' @rdname a4aFitSA-class
setMethod("a4aFitSAs", signature(object="list"),
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
    args <- c(list(Class="a4aFitSAs", .Data=object, names=names),
      args[!names(args)%in%'names'])

    return(
      do.call('new', args)
      )

})

#' @rdname a4aFitSA-class
setMethod("a4aFitSAs", signature(object="a4aFitSA"), function(object, ...) {
    lst <- c(object, list(...))
    a4aFitSAs(lst)
})

#' @rdname a4aFitSA-class
setMethod("a4aFitSAs", signature(object="missing"),
  function(...) {
    # empty
    if(missing(...)){
      new("a4aFitSAs")
    # or not
    } else {
      args <- list(...)
      object <- args[!names(args)%in%c('names', 'desc', 'lock')]
      args <- args[!names(args)%in%names(object)]
      do.call('a4aFitSAs',  c(list(object=object), args))
    }
  }
)


