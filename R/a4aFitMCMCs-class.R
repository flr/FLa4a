#====================================================================
# plural class for a4aFitMCMC (used for model averaging)
#====================================================================

#' @rdname a4aFitMCMC-class
#' @aliases a4aFitMCMCs-class

setClass("a4aFitMCMCs",
  contains="a4aFitSAs"
)

setValidity("a4aFitMCMCs",
  function(object) {
    if(!all(sapply(object, is, 'a4aFitMCMC'))) {
      "Components must be a4aFitMCMC"
    } else {
      TRUE
    }
})


#' @rdname a4aFitMCMC-class
#' @aliases a4aFitMCMCs a4aFitMCMCs-methods
setGeneric("a4aFitMCMCs", function(object, ...) standardGeneric("a4aFitMCMCs"))

#' @rdname a4aFitMCMC-class
setMethod("a4aFitMCMCs", signature(object="list"),
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
    args <- c(list(Class="a4aFitMCMCs", .Data=object, names=names),
      args[!names(args)%in%'names'])

    return(
      do.call('new', args)
      )

})

#' @rdname a4aFitMCMC-class
setMethod("a4aFitMCMCs", signature(object="a4aFitMCMC"), function(object, ...) {
    lst <- c(object, list(...))
    a4aFitMCMCs(lst)
})

#' @rdname a4aFitMCMC-class
setMethod("a4aFitMCMCs", signature(object="missing"),
  function(...) {
    # empty
    if(missing(...)){
      new("a4aFitMCMCs")
    # or not
    } else {
      args <- list(...)
      object <- args[!names(args)%in%c('names', 'desc', 'lock')]
      args <- args[!names(args)%in%names(object)]
      do.call('a4aFitMCMCs',  c(list(object=object), args))
    }
  }
)

#' @title Plot of metrics of multiple fits
#' @name plot metrics of multiple fits
#' @docType methods
#' @rdname plot-mmcfits
#' @aliases plot,a4aFitMCMCs, missing-method
#' @description Method to plot fitting statistics of multiple fits, useful to compare fits.
#' @param x an \code{a4aFitMCMCs} object with multiple fits
#' @param y ignored
#' @param ... additional argument list that might never be used
#' @return a \code{plot} with fitting statistics
#' @examples
#' data(ple4)
#' data(ple4.index)
#' qmods <- list(list(~s(age, k=6)))
#' fmods = list()
#' for(i in 1:6) {
#'   fmods[[paste0(i)]] <-
#'     as.formula(
#'       paste0("~te(age, year, k = c(6,", i+14,"), bs = 'tp') + s(age, k = 6)")
#'     )
#' }
#' myFits <-
#'   scax(
#'     FLStocks(ple4), list(FLIndices(ple4.index)),
#'     fmodel = fmods, qmodel=qmods, fit="MCMC", mcmc=SCAMCMC()
#'   )
#' plot(myFits)

setMethod("plot", c("a4aFitMCMCs", "missing"), function(x, y=missing, ...){
	args <- list()
	par(mar=c(5, 4, 4, 4) + 0.1)
	accrate = lapply(x,function(x) fitSumm(x)['accrate',])
	df <- data.frame(unlist(accrate))
	df$fit <- as.numeric(gsub("fit", "",names(accrate)))
	names(df) <- c("accrate","fit")
	plot(df$fit, df$accrate, type = "b", col = "blue", ylab = "accrate", xlab = "fit", main="Analysis of fit metrics", ylim=c(0.1,0.5))
	abline(h=0.3, col = "blue",lty = 2)
})
