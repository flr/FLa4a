# a4aFitResiduals-class - «Short one line description»
# a4aFitResiduals-class

#' a4aFitResiduals extends \code{"FLQuant"} class.
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
#' @name a4aFitResiduals-class
#' @rdname a4aFitResiduals-class
#' @exportClass a4aFitResiduals
setClass("a4aFitResiduals", contain="FLQuants")

# constructor

#' Computes log residuals of catches and indices 
#'
#' @param a4aFit object
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitResiduals-methods
#'
#' @name a4aFit
#' @docType methods
#' @rdname a4aFit-methods
#' @aliases a4aFit,a4aFit-method
#' @examples
#' library(FLa4a)
#' data(ple4)
#' data(ple4.index)
#' fmodel <- ~factor(age) + factor(year)
#' qmodel <- list(~factor(age))
#' fit. <- a4a(stock=ple4, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(ple4.index), fit ="assessment", wkdir="test")
#' flqs <- residuals(fit., ple4, FLIndices(idx=ple4.index))

setMethod("residuals", signature(object="a4aFit"), function(object, stock, indices, ...) {
	args <- list(...)
	# object holder
	lst <- list()
	length(lst) <- length(indices) + 1	
	# catch
	lst[[1]] <- stdlogres(catch.n(stock), catch.n(object))
	# indices
	idx <- index(fit.)
	for(i in 1:length(indices)){
		lst[[i+1]] <- stdlogres(index(indices[[i]]), idx[[i]])
	}
	# out
	names(lst) <- c("catch.n", names(indices))
	new("a4aFitResiduals", FLQuants(lst))
  }
)


# std log residuals

#' Standardized log residuals 
#' @name stdlogres
#' @docType methods
#' @rdname stdlogres-methods
#' @aliases stdlogres,stdlogres-method
#' Method to get K values
#'
#' @param obs a \code{FLQuant} object with the observations
#' @param fit a \code{FLQuant} object with the fitted value
#' @return a \code{FLQuant} with stardardized log residuals
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
#' @examples
#' library(FLa4a)
#' data(ple4)
#' data(ple4.index)
#' fmodel <- ~factor(age) + factor(year)
#' qmodel <- list(~factor(age))
#' fit. <- a4a(stock=ple4, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(ple4.index), fit ="assessment", wkdir="test")
#' flqs <- residuals(fit., ple4, FLIndices(idx=ple4.index))

setGeneric("stdlogres", function(obs, fit, ...) standardGeneric("stdlogres"))

setMethod("stdlogres", c("FLQuant","FLQuant"), function(obs, fit, ...){

	flq <- log(obs/fit)	
	res <- apply(flq, 2:6, scale, center=FALSE)
	dimnames(res) <- dimnames(flq)
	FLQuant(res)
	
}) 


#' plot of standardized log residuals 
#' @name plot
#' @docType methods
#' @rdname plot-methods
#' @aliases plot,plot-method
#' Method to plot scatterplot of standardized residuals
#'
#' @param x a \code{a4aFitResiduals} object with the standardized residuals
#' @param ... Additional argument list that might not ever be used.
#' @return a \code{plot} with stardardized log residuals
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
#' @examples
#' library(FLa4a)
#' data(ple4)
#' data(ple4.index)
#' fmodel <- ~factor(age) + factor(year)
#' qmodel <- list(~factor(age))
#' fit. <- a4a(stock=ple4, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(ple4.index), fit ="assessment", wkdir="test")
#' stdlogres(catch.n(ple4), catch.n(fit.))
#' plot(flqs)

setMethod("plot", c("a4aFitResiduals", "missing"), function(x, y=missing, ...){

	xyplot(data~year|factor(age)*qname, type=c("p", "smooth"), groups=qname, data=x, cex=0.3, lwd=2, ylab="standardized residuals", xlab="", panel=function(x,y,...){panel.abline(h=0, col.line="gray80");panel.xyplot(x,y,...)}, par.settings=list(superpose.symbol=list(col="gray50", pch=19, cex=0.2), superpose.line=list(col=1, lty=1, lwd=2), strip.background=list(col="gray90"), strip.border=list(col="gray90"), box.rectangle=list(col="gray90")), main="log residuals of catch and abundance indices", ...)

})

#' qqplot of standardized log residuals 
#' @name plot
#' @docType methods
#' @rdname qqmath-methods
#' @aliases qqmath,qqmath-method
#' Method to plot qqplots of standardized residuals
#'
#' @param x a \code{a4aFitResiduals} object with the standardized residuals
#' @param ... Additional argument list that might not ever be used.
#' @return a \code{qqplot} with stardardized log residuals
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
#' @examples
#' library(FLa4a)
#' data(ple4)
#' data(ple4.index)
#' fmodel <- ~factor(age) + factor(year)
#' qmodel <- list(~factor(age))
#' fit. <- a4a(stock=ple4, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(ple4.index), fit ="assessment", wkdir="test")
#' stdlogres(catch.n(ple4), catch.n(fit.))
#' qqmath(flqs)

setGeneric("qqmath", function(x, data, ...) standardGeneric("qqmath"))

setMethod("qqmath", c("a4aFitResiduals", "missing"), function(x, data=missing, ...){

	qqmath(~data|factor(age)*qname, data=as.data.frame(x), ylab="standardized residuals", xlab="", prepanel=prepanel.qqmathline, panel = function(x, ...){panel.qqmathline(x, col="gray50"); panel.qqmath(x, ...)}, col=1, pch=19, cex=0.2, par.settings=list(strip.background=list(col="gray90"), strip.border=list(col="gray90"), box.rectangle=list(col="gray90")), main="quantile-quantile plot of log residuals of catch and abundance indices", ...)

})

#' bubbles plot of standardized log residuals 
#' @name plot
#' @docType methods
#' @rdname bubbles-methods
#' @aliases bubbles,bubbles-method
#' Method to plot bubbles of standardized residuals
#'
#' @param x a \code{a4aFitResiduals} object with the standardized residuals
#' @param ... Additional argument list that might not ever be used.
#' @return a \code{bubbles} plot with stardardized log residuals
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
#' @examples
#' library(FLa4a)
#' data(ple4)
#' data(ple4.index)
#' fmodel <- ~factor(age) + factor(year)
#' qmodel <- list(~factor(age))
#' fit. <- a4a(stock=ple4, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(ple4.index), fit ="assessment", wkdir="test")
#' stdlogres(catch.n(ple4), catch.n(fit.))
#' bubbles(flqs)

setMethod("bubbles", c("a4aFitResiduals", "missing"), function(x, data=missing, ...){

	bubbles(age~year|qname, data=x, par.settings=list(strip.background=list(col="gray90"), strip.border=list(col="gray90"), box.rectangle=list(col="gray90")), main="log residuals of catch and abundance indices", ...)

})


