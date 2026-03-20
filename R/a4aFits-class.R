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
	gcv = lapply(myFits,function(x) fitSumm(x)['gcv',])
	bic = lapply(myFits, function(x) BIC(x))
	df <- data.frame(unlist(gcv), unlist(bic))
	df$fit <- as.numeric(gsub("fit", "",names(gcv)))
	names(df) <- c("GCV","BIC","fit")
	col_gcv <- "#0072B2"  
	col_bic <- "#D55E00"  
	par(mar = c(5, 5, 4, 5)) 
	gcv_lims <- range(pretty(df$GCV))
	bic_lims <- range(pretty(df$BIC))
	plot(df$fit, df$GCV, type = "n", axes = FALSE, 
	     xlab = "", ylab = "", main = "Analysis of Fit Metrics",
	     ylim = gcv_lims) 
	grid(nx = NA, ny = NULL, col = "gray90", lty = 1, lwd = 1)
	lines(df$fit, df$GCV, col = col_gcv, lwd = 2.5)
	points(df$fit, df$GCV, col = col_gcv, pch = 16, cex = 1.2)
	axis(2, col = col_gcv, col.axis = col_gcv, las = 1, lwd = 1.5)
	mtext("GCV", side = 2, line = 3.5, col = col_gcv, font = 2)
	par(new = TRUE)
	plot(df$fit, df$BIC, type = "n", axes = FALSE, xlab = "", ylab = "",
	     ylim = bic_lims) 
	lines(df$fit, df$BIC, col = col_bic, lwd = 2.5)
	points(df$fit, df$BIC, col = col_bic, pch = 16, cex = 1.2)
	axis(4, col = col_bic, col.axis = col_bic, las = 1, lwd = 1.5)
	mtext("BIC", side = 4, line = 3.5, col = col_bic, font = 2)
	axis(1, at = df$fit, col = "gray30", col.axis = "gray30", lwd = 1.5)
	mtext("Fit Index", side = 1, line = 3, font = 2, col = "gray30")
	abline(v = df$fit[which.min(df$GCV)], col = col_gcv, lty = 2, lwd = 2)
	abline(v = df$fit[which.min(df$BIC)], col = col_bic, lty = 2, lwd = 2)
	legend("bottomleft", legend = c("GCV", "BIC"), 
	       col = c(col_gcv, col_bic), lwd = 2.5, pch = 16, 
	       bty = "n", inset = 0.02)
})
