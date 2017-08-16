#' @title S4 class \code{a4aFitSA}
#'
#' @description The \code{a4aFitSA} class extends \code{a4aFit} to store information about the parameters of the model.
#'
#' @section Slots:
#' \describe{
#'    \item{call}{The function call}
#'    \item{clock}{Information on call duration}
#'    \item{fitSumm}{Fit summary}
#'    \item{stock.n}{Estimates of stock numbers-at-age}
#'    \item{harvest}{Estimates of fishing mortality at age}
#'    \item{catch.n}{Estimates of catch numbers-at-age}
#'    \item{index}{Estimates of survey or CPUE indices-at-age}
#'    \item{pars}{an object of class \code{SCAPars} with information about model parameters}
#' }
#'
#' @template Accessors
#' @template Constructors
#' @docType class
#' @name a4aFitSA-class
#' @rdname a4aFitSA-class
#' @aliases a4aFitSA-class
#' @examples
#' data(ple4)
#' data(ple4.index)
#'
#' obj <- sca(stock=ple4, indices=FLIndices(ple4.index), fit="assessment")
#' obj
#'
#' slotNames(obj)
#' clock(obj)
#' fitSumm(obj)
#'
#' flq <- stock.n(obj)
#' is(flq)
#' flq <- index(obj)
#' is(flq)
#'
#' logLik(obj)
#' AIC(obj)
#' BIC(obj)
#'
#' is(pars(obj))

a4aFitSA <-
  setClass("a4aFitSA",
           contains = "a4aFit",
           slots = c(pars = "SCAPars"))

setValidity("a4aFitSA",
  function(object) {
    # no validation at present
    TRUE
})

setMethod("initialize", "a4aFitSA",
    function(.Object, ..., pars) {
      if (!missing(pars)) .Object@pars <- pars
      .Object <- callNextMethod(.Object, ...)
      .Object
})

#' @rdname a4aFitSA-class
setMethod("show", "a4aFitSA",
  function(object)
  {
    show(a4aFit(object))

    cat("\nSubmodels:\n")
    submodels(object)
 })


#' @rdname a4aFitSA-class
setMethod("submodels", signature(object = "a4aFitSA"),
  function(object)
  {
    cat("\t fmodel: "); print( fmodel(pars(object)), showEnv = FALSE)
    cat("\tsrmodel: "); print(srmodel(pars(object)), showEnv = FALSE)
    cat("\tn1model: "); print(n1model(pars(object)), showEnv = FALSE)

    # something to format the qmodel and vmodel formulas
    printFormulaList <- function(mods) {
      if (length(mods) == 0) return(invisible(NA))
      mnames <- names(mods)
      maxname <- max(sapply(mnames, nchar))
      for (i in seq_along(mods)) {
        cat("\t   ", mnames[i], ": ", rep(" ", maxname - nchar(mnames[i])), sep = "")
        print(mods[[i]], showEnv = FALSE)
      }
    }

    cat("\t qmodel:\n")
    printFormulaList(sMod(qmodel(pars(object))))
    cat("\t vmodel:\n")
    printFormulaList(sMod(vmodel(pars(object))))

 })


#' @rdname a4aFitSA-class
#' @template bothargs
#' @aliases a4aFitSA a4aFitSA-methods
setGeneric("a4aFitSA")

#' @rdname a4aFitSA-class
setMethod("a4aFit", "a4aFitSA",
  function(object, ...) {
    as(object, "a4aFit")
  }
)


#' @rdname a4aFitSA-class
#' @aliases pars pars-methods
setGeneric("pars", function(object) standardGeneric("pars"))
#' @rdname a4aFitSA-class
setMethod("pars", "a4aFitSA", function(object) object@pars)

#' @rdname a4aFitSA-class
setMethod("m", signature(object="a4aFitSA"), function(object) m(pars(object)))

#' @rdname a4aFitSA-class
setMethod("wt", signature(object="a4aFitSA"), function(object) wt(pars(object)))

#' @rdname a4aFitSA-class
#' @param obj the object to be subset
#' @param it iteration to be extracted
setMethod("iter", "a4aFitSA", function(obj, it){
	obj@fitSumm <- obj@fitSumm[,it, drop=FALSE]
	obj@harvest <- iter(obj@harvest, it)
	obj@stock.n <- iter(obj@stock.n, it)
	obj@catch.n <- iter(obj@catch.n, it)
	obj@index <- iter(obj@index, it)
	obj@pars <- iter(obj@pars, 1)
	obj
})


#====================================================================
# plural class for a4aFitSA (used for model averaging)
#====================================================================

#' @rdname a4aFitSA-class
#' @aliases a4aFitSAs-class

setClass("a4aFitSAs",
	contains="FLComps",
	validity=function(object){
		if(!all(unlist(lapply(object, is, 'a4aFitSA'))))
			return("Components must be a4aFitSA")
		return(TRUE)}
)

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



