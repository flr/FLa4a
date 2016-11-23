#' @title S4 class \code{a4aFitResiduals}
#' @description The \code{a4aFitResiduals} class extends \code{FLQuants} to store residuals of the a4a stock assessment fit. By default, these should be log residuals of catches and indices.
#' @docType class
#' @name a4aFitResiduals-class
#' @rdname a4aFitResiduals-class
#' @aliases a4aFitResiduals-class

setClass("a4aFitResiduals", contain="FLQuants")

#' @rdname a4aFitResiduals-class
#' @aliases a4aFitResiduals a4aFitResiduals-methods residuals,a4aFit-method
#' @template bothargs
#' @param stock \code{FLStock} object used to fit the model
#' @param indices \code{FLIndices} object used to fit the model
#' @examples
#' data(ple4)
#' data(ple4.index)
#' obj <- sca(ple4, FLIndices(ple4.index))
#' flqs <- residuals(obj, ple4, FLIndices(idx=ple4.index))
setMethod("residuals", signature(object="a4aFit"), function(object, stock, indices, ...) {
	args <- list(...)
	# object holder
	lst <- list()
	length(lst) <- length(indices) + 2	
	# catch
	lst[[1]] <- stdlogres(catch.n(stock), catch.n(object))
	# indices
	idx <- index(object)
	for(i in 1:length(indices)){
		lst[[i+1]] <- stdlogres(index(indices[[i]]), idx[[i]])
	}
	lst[[length(lst)]] <- stdlogres(catch(stock), computeCatch(stock + object))
	# out
	names(lst) <- c("catch.n", names(indices), "catch")
	new("a4aFitResiduals", FLQuants(lst))
  }
)

#' @title Standardized log residuals 
#' @description Method to compute the standardized residuals on the log scale for index- and catch-at-age residuals in the a4a stock assessment framework.
#' @name stdlogres
#' @docType methods
#' @rdname stdlogres-methods
#' @aliases stdlogres stdlogres-methods stdlogres,FLQuant,FLQuant-method
#' @param obs an \code{FLQuant} object with the observations
#' @param fit an \code{FLQuant} object with the fitted value
#' @template dots
#' @return an \code{FLQuant} with stardardized log residuals
#' @examples
#' data(ple4)
#' data(ple4.index)
#' obj <- sca(ple4, FLIndices(ple4.index))
#' flqs <- residuals(obj, ple4, FLIndices(idx=ple4.index))
#'
#' stdlogres(catch.n(ple4), catch.n(obj))
#' # which is the same as the following (because residuals() uses stdlogres):
#' flqs$catch.n
#' # check:
#' stdlogres(catch.n(ple4),catch.n(obj)) - flqs$catch.n
setGeneric("stdlogres", function(obs, fit, ...) standardGeneric("stdlogres"))

#' @rdname stdlogres-methods
setMethod("stdlogres", c("FLQuant","FLQuant"), function(obs, fit, ...){
	flq <- log(obs/fit)	
	#res <- apply(flq, c(1,3:6), scale, center=FALSE)
	res <- flq/sqrt(yearVars(flq)[,rep(1,ncol(flq))])
	dimnames(res) <- dimnames(flq)
	as(res, "FLQuant")
}) 

#' @title Plot of standardized log residuals
#' @name plot
#' @docType methods
#' @rdname plot-methods
#' @aliases plot,a4aFitResiduals,missing-method
#' @description Method to produce scatterplots of standardized residuals
#' @param x an \code{a4aFitResiduals} object with the standardized residuals
#' @param y ignored
#' @param ... additional argument list that might never be used
#' @return a \code{plot} with stardardized log residuals
#' @examples
#' data(ple4)
#' data(ple4.index)
#' obj <- sca(ple4, FLIndices(ple4.index))
#' flqs <- residuals(obj, ple4, FLIndices(idx=ple4.index))
#' plot(flqs)

setMethod("plot", c("a4aFitResiduals", "missing"), function(x, y=missing, ...){
	args <- list()
	args$data <- as.data.frame(x)
	args$x <- as.formula("data~year|factor(age)*qname")
	args$type=c("p", "smooth")
	args$groups <- quote(qname)
	args$cex=0.3
	args$lwd=2
	args$ylab="standardized residuals"
	args$xlab=""
	args$panel=function(x,y,...){
		panel.abline(h=0, col.line="gray80")
		panel.xyplot(x,y,...)
		}
	args$par.settings=list(
		superpose.symbol=list(col="gray50", pch=19, cex=0.2), 
		superpose.line=list(col=1, lty=1, lwd=2), 
		strip.background=list(col="gray90"), 
		strip.border=list(col="black"), 
		box.rectangle=list(col="gray90"))
	args$main="log residuals of catch and abundance indices by age"
	if(is(latticeExtra::useOuterStrips, "function")) latticeExtra::useOuterStrips(do.call("xyplot", args)) else do.call("xyplot", args)
})

#' @title qqplot of standardized log residuals
#' @name qqmath
#' @docType methods
#' @rdname qqmath-methods
#' @aliases qqmath,a4aFitResiduals,missing-method
#' @description Method to produce qqplots of standardized residuals
#' @param x an \code{a4aFitResiduals} object with the standardized residuals
#' @param data ignored
#' @param ... additional argument list that might never be used
#' @return a \code{qqplot} with stardardized log residuals
#' @examples
#' data(ple4)
#' data(ple4.index)
#' obj <- sca(ple4, FLIndices(ple4.index))
#' flqs <- residuals(obj, ple4, FLIndices(idx=ple4.index))
#' qqmath(flqs)

setGeneric("qqmath", function(x, data, ...) standardGeneric("qqmath"))
#' @rdname qqmath-methods
setMethod("qqmath", c("a4aFitResiduals", "missing"), function(x, data=missing, ...){
	qqmath(~data|factor(age)*qname, data=as.data.frame(x), ylab="standardized residuals", xlab="", prepanel=prepanel.qqmathline, panel = function(x, ...){panel.qqmathline(x, col="gray50"); panel.qqmath(x, ...)}, col=1, pch=19, cex=0.2, par.settings=list(strip.background=list(col="gray90"), strip.border=list(col="gray90"), box.rectangle=list(col="gray90")), main="quantile-quantile plot of log residuals of catch and abundance indices", ...)
})

#' @title Bubbles plot of standardized log residuals
#' @name bubbles
#' @docType methods
#' @rdname bubbles-methods
#' @aliases bubbles,a4aFitResiduals,missing-method
#' @description Method to produce bubble plots of standardized residuals
#' @param x an \code{a4aFitResiduals} object with the standardized residuals
#' @param data ignored
#' @param ... additional argument list that might never be used
#' @return a \code{bubbles} plot with stardardized log residuals
#' @examples
#' data(ple4)
#' data(ple4.index)
#' obj <- sca(ple4, FLIndices(ple4.index))
#' flqs <- residuals(obj, ple4, FLIndices(idx=ple4.index))
#' bubbles(flqs)
setMethod("bubbles", c("a4aFitResiduals", "missing"), function(x, data=missing, ...){
	bubbles(age~year|qname, data=x, par.settings=list(strip.background=list(col="gray90"), strip.border=list(col="gray90"), box.rectangle=list(col="gray90")), main="log residuals of catch and abundance indices", ...)
})


