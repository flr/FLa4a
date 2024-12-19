#====================================================================
# plural class for a4aFit (used for model averaging)
#====================================================================

#' @rdname a4aFit-class
#' @aliases a4aFits-class

setClass("a4aFits",
  contains="FLComps"
)

setValidity("a4aFits",
  function(object) {
    if(!all(sapply(object, is, 'a4aFit'))) {
      "Components must be a4aFit"
    } else {
      TRUE
    }
})


#' @rdname a4aFit-class
#' @aliases a4aFits a4aFits-methods
setGeneric("a4aFits", function(object, ...) standardGeneric("a4aFits"))

#' @rdname a4aFit-class
setMethod("a4aFits", signature(object="list"),
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
    args <- c(list(Class="a4aFits", .Data=object, names=names),
      args[!names(args)%in%'names'])

    return(
      do.call('new', args)
      )

})

#' @rdname a4aFit-class
setMethod("a4aFits", signature(object="a4aFit"), function(object, ...) {
    lst <- c(object, list(...))
    a4aFits(lst)
})

#' @rdname a4aFit-class
setMethod("a4aFits", signature(object="missing"),
  function(...) {
    # empty
    if(missing(...)){
      new("a4aFits")
    # or not
    } else {
      args <- list(...)
      object <- args[!names(args)%in%c('names', 'desc', 'lock')]
      args <- args[!names(args)%in%names(object)]
      do.call('a4aFits',  c(list(object=object), args))
    }
  }
)


#' @title Plot of metrics of multiple fits
#' @name plot metrics of multiple fits
#' @docType methods
#' @rdname plot-mfits
#' @aliases plot,a4aFits, missing-method
#' @description Method to plot fitting statistics of multiple fits, useful to compare fits.
#' @param x an \code{a4aFits} object with multiple fits
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
#'   scas(
#'     FLStocks(ple4), list(FLIndices(ple4.index)),
#'     fmodel = fmods, qmodel=qmods, fit="MP"
#'   )
#' plot(myFits)

setMethod("plot", c("a4aFits", "missing"), function(x, y=missing, ...){
	args <- list()
	par(mar=c(5, 4, 4, 4) + 0.1)
	gcv = lapply(x,function(x) fitSumm(x)['gcv',])
	bic = lapply(x, function(x) BIC(x))
	df <- data.frame(unlist(gcv), unlist(bic))
	df$fit <- as.numeric(gsub("fit", "",names(gcv)))
	names(df) <- c("GCV","BIC","fit")
	df <- df[complete.cases(df),]
	plot(df$fit, df$GCV, type = "b", col = "blue", ylab = "GCV", xlab = "fit", main="Analysis of fit metrics")
	par(new = TRUE)
	plot(df$fit, df$BIC, type = "b", col = "red", axes = FALSE, xlab = "", ylab = "")
	axis(4)
	mtext("BIC", side=4, line=3)
	abline(v=df[min(df$GCV)==df$GCV,]$fit, col = "blue",lty = 2)
	abline(v=df[min(df$BIC)==df$BIC,]$fit, col = "red",lty = 2)
	legend("topleft", legend = c("GCV", "BIC"), col = c("blue", "red"), lty = 1, bg="white")
})
