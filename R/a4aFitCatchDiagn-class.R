
#' @title S4 class \code{a4aFitCatchDiagn}
#' @description The \code{a4aFitCatchDiagn} class extends \code{FLQuants} to store information to run diagnostics on aggregated catch estimated  by the a4a stock assessment fit.
#' @docType class
#' @name a4aFitCatchDiagn-class
#' @rdname a4aFitCatchDiagn-class
#' @aliases a4aFitCatchDiagn-class

setClass("a4aFitCatchDiagn", contain="FLQuants")
# VALIDITY CHECK NAMES


#' @rdname a4aFitCatchDiagn-class
#' @aliases a4aFitCatchDiagn a4aFitCatchDiagn-methods
#' @template bothargs
#' @param stock \code{FLStock} object used to fit the model
#' @param indices \code{FLIndices} object used to fit the model
#' @examples
#' data(ple4)
#' data(ple4.index)
#' fit <- sca(ple4, ple4.index)
#' flqs <- computeCatchDiagnostics(fit, ple4)
setGeneric("computeCatchDiagnostics", function(object, ...) standardGeneric("computeCatchDiagnostics"))

#' @rdname a4aFit-class
#' @aliases clock clock-methods
setMethod("computeCatchDiagnostics", signature(object="a4aFit"), function(object, stock, ...) {
	args <- list(...)
	it <- 500
	pred <- predict(object)
	if(dims(object)$iter==1) fits <- simulate(object, it) else it <- dims(object)$iter
	stk_new <- stock + object
	cth_est <- catch(stk_new)
	cth_obs <- catch(stock)
	# sim with observation error
	flqoe <- quantSums(rlnorm(it, log(catch.n(object)), pred$vmodel$catch)*catch.wt(stock))
	# sim with estimation error
	flqee <- catch(stock + fits)
	# both
	flqe <- quantSums(rlnorm(it, log(catch.n(object)), sqrt(pred$vmodel$catch^2 + iterVars(log(catch.n(fits)))))*catch.wt(stock))
	# standardized residuals
	resstd <- stdlogres(cth_obs, cth_est)
	# pearson residuals
	sdlog <- sqrt(iterVars(log(flqoe)))
	resprs <- stdlogres(cth_obs, cth_est, sdlog=sdlog)
	# deviances
	devs <- log(cth_obs/cth_est)
	# out
	flqs <- FLQuants(list(obs = cth_obs, est = cth_est, oe=flqoe, ee=flqee, oee=flqe, resstd=resstd, resprs=resprs, resraw=devs))
	new("a4aFitCatchDiagn", flqs)
  }
)

#' @title Plot of aggregated catch standardized log residuals
#' @name plot of catch residuals
#' @docType methods
#' @rdname plot-catch
#' @aliases plot,a4aFitCatchDiagn,missing-method
#' @description Method to produce scatterplots of aggregated catch residuals
#' @param x an \code{a4aFit} object with the model fit
#' @param y the \code{FLStock} object used to fit the model
#' @param ... additional argument list that might never be used
#' @return a \code{plot} with stardardized log residuals
#' @examples
#' data(ple4)
#' data(ple4.index)
#' fit <- sca(ple4, ple4.index)
#' flqs <- computeCatchDiagnostics(fit, ple4)
#' plot(flqs)

setMethod("plot", c("a4aFitCatchDiagn", "missing"), function(x, y=missing, probs=c(0.1, 0.9), type="all", ...){
	args <- list()

	#----------------------------------------------------------------
	# lattice stores the code for the plot not the plot itself, which means
	# one has to have datasets for each plot and formulas must match
	#----------------------------------------------------------------

	#----------------------------------------------------------------
	# absolute catches
	#----------------------------------------------------------------
	
	# build datasets
	ci <- probs[2] - probs[1]
	probs <- c(probs[1], 0.5, probs[2])
	d1 <- as.data.frame(quantile(x$oe, probs))
	d2 <- as.data.frame(quantile(x$ee, probs))
	d3 <- as.data.frame(quantile(x$oee, probs))
	obs <- x$obs
	d1$iter <- factor(d1$iter, labels=c("l","m","u"))
	d2$iter <- factor(d2$iter, labels=c("l","m","u"))
	d3$iter <- factor(d3$iter, labels=c("l","m","u"))

	# these are absolute catches which should be in a similar scale
	# ylimits across catch plots
	v1 <- c(d1$data, d2$data, d3$data)
	mn1 <- min(v1)*.95
	mx1 <- max(v1)*1.05

	# arguments for xyplot
	pset <- list(regions=list(col="gray95"), axis.line = list(col = "gray75"))
	pfun <- function(x,y,subscripts,groups,...){
		panel.polygon(c(x[groups=="l"],rev(x[groups=="u"])), c(y[groups=="l"],rev(y[groups=="u"])), col="gray85", border=0)
		panel.grid(col="gray95")
		panel.xyplot(x[groups=="m"], y[groups=="m"], lty=2, col=1, lwd=1.5, ...)
		panel.xyplot(dimnames(obs)[[2]], c(obs), type="l", col=1, lwd=1.5)
	}

	# oe
	p1 <- xyplot(data~year, groups=iter, data=d1, type="l", par.settings=pset, panel=pfun, xlab="", ylab="", main="Observation error", ylim=c(mn1, mx1))

	# ee
	p2 <- xyplot(data~year, groups=iter, data=d2, type="l", par.settings=pset, panel=pfun, xlab="", ylab="", main="Estimation error", ylim=c(mn1, mx1)) 

	# oee
	p3 <- xyplot(data~year, groups=iter, data=d3, type="l", par.settings=pset, panel=pfun, xlab="", ylab="", main="Prediction error", ylim=c(mn1, mx1)) 

	#----------------------------------------------------------------
	# catch residuals
	#----------------------------------------------------------------

	# these are residuals, which should be around a N(0,1) scale, and
	# have negative values
	# ylimits across residual plots
	v2 <- max(abs(unlist(x[c("resprs","resstd","resraw")])))*1.05
	mn2 <- -v2
	mx2 <- v2

	# build dataset
	y1 <- x1 <- x2 <- as.numeric(dimnames(x$obs)[[2]])
	y1[] <- 0
	df0 <- data.frame(x1=x1, x2=x2, y1=y1, resprs=c(x$resprs), resstd=c(x$resstd), resraw=c(x$resraw))

	# arguments for xyplot
	args <- list(...)
	args$x = y1~x1
	args$data = df0
	args$type="p"
	args$cex=0.5
	args$pch=19
	args$col=1
	args$ylim=c(mn2, mx2)
	args$ylab=""
	args$xlab=""
	args$par.settings=list(axis.line = list(col = "gray75"))

	# pearson
	args$panel=function(x1, y1, ...){
		panel.grid(col="gray95")
		panel.segments(df0$x1, df0$y1, df0$x2, df0$resprs, ...)
		panel.xyplot(df0$x1, df0$resprs, ...)
		}
	args$main="Pearson residuals"
	p4 <- do.call("xyplot", args)

	# standardized
	args$panel=function(x1, y1, ...){
		panel.grid(col="gray95")
		panel.segments(df0$x1, df0$y1, df0$x2, df0$resstd, ...)
		panel.xyplot(df0$x1, df0$resstd, ...)
		}
	args$main="Standardized residuals"
	p5 <- do.call("xyplot", args)

	# deviances
	args$panel=function(x1, y1, ...){
		panel.grid(col="gray95")
		panel.segments(df0$x1, df0$y1, df0$x2, df0$resraw, ...)
		panel.xyplot(df0$x1, df0$resraw, ...)
		}
	args$main="Raw residuals (deviances)"
	p6 <- do.call("xyplot", args)
	
	#----------------------------------------------------------------
	# build plot
	#----------------------------------------------------------------

	subtext <- paste("(shaded area = CI", ci*100, "%, dashed line = median, solid line = observed \n", sep="")

	if(type=="prediction"){
		p3$sub <- list(label=subtext, cex=1)
		p3$main <- list(label="Prediction error" , cex=1.5)
		print(p3)
	}
	if(type=="all"){
		grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3, as.table=FALSE, top=textGrob("Aggregated catch diagnostics \n", gp=gpar(fontface = "bold", cex=1.5)), bottom=textGrob(subtext, gp=gpar(cex=1)))
	}
})

