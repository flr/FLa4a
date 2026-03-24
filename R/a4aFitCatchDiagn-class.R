
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

# setMethod("plot", c("a4aFitCatchDiagn", "missing"), function(x, y=missing, probs=c(0.1, 0.9), type="all", ...){
# 	args <- list()
#
# 	#----------------------------------------------------------------
# 	# lattice stores the code for the plot not the plot itself, which means
# 	# one has to have datasets for each plot and formulas must match
# 	#----------------------------------------------------------------
#
# 	#----------------------------------------------------------------
# 	# absolute catches
# 	#----------------------------------------------------------------
#
# 	# build datasets
# 	ci <- probs[2] - probs[1]
# 	probs <- c(probs[1], 0.5, probs[2])
# 	d1 <- as.data.frame(quantile(x$oe, probs))
# 	d2 <- as.data.frame(quantile(x$ee, probs))
# 	d3 <- as.data.frame(quantile(x$oee, probs))
# 	obs <- x$obs
# 	d1$iter <- factor(d1$iter, labels=c("l","m","u"))
# 	d2$iter <- factor(d2$iter, labels=c("l","m","u"))
# 	d3$iter <- factor(d3$iter, labels=c("l","m","u"))
#
# 	# these are absolute catches which should be in a similar scale
# 	# ylimits across catch plots
# 	v1 <- c(d1$data, d2$data, d3$data)
# 	mn1 <- min(v1)*.95
# 	mx1 <- max(v1)*1.05
#
# 	# arguments for xyplot
# 	pset <- list(regions=list(col="gray95"), axis.line = list(col = "gray75"))
# 	pfun <- function(x,y,subscripts,groups,...){
# 		panel.polygon(c(x[groups=="l"],rev(x[groups=="u"])), c(y[groups=="l"],rev(y[groups=="u"])), col="gray85", border=0)
# 		panel.grid(col="gray95")
# 		panel.xyplot(x[groups=="m"], y[groups=="m"], lty=2, col=1, lwd=1.5, ...)
# 		panel.xyplot(dimnames(obs)[[2]], c(obs), type="l", col=1, lwd=1.5)
# 	}
#
# 	# oe
# 	p1 <- xyplot(data~year, groups=iter, data=d1, type="l", par.settings=pset, panel=pfun, xlab="", ylab="", main="Observation error", ylim=c(mn1, mx1))
#
# 	# ee
# 	p2 <- xyplot(data~year, groups=iter, data=d2, type="l", par.settings=pset, panel=pfun, xlab="", ylab="", main="Estimation error", ylim=c(mn1, mx1))
#
# 	# oee
# 	p3 <- xyplot(data~year, groups=iter, data=d3, type="l", par.settings=pset, panel=pfun, xlab="", ylab="", main="Prediction error", ylim=c(mn1, mx1))
#
# 	#----------------------------------------------------------------
# 	# catch residuals
# 	#----------------------------------------------------------------
#
# 	# these are residuals, which should be around a N(0,1) scale, and
# 	# have negative values
# 	# ylimits across residual plots
# 	v2 <- max(abs(unlist(x[c("resprs","resstd","resraw")])))*1.05
# 	mn2 <- -v2
# 	mx2 <- v2
#
# 	# build dataset
# 	y1 <- x1 <- x2 <- as.numeric(dimnames(x$obs)[[2]])
# 	y1[] <- 0
# 	df0 <- data.frame(x1=x1, x2=x2, y1=y1, resprs=c(x$resprs), resstd=c(x$resstd), resraw=c(x$resraw))
#
# 	# arguments for xyplot
# 	args <- list(...)
# 	args$x = y1~x1
# 	args$data = df0
# 	args$type="p"
# 	args$cex=0.5
# 	args$pch=19
# 	args$col=1
# 	args$ylim=c(mn2, mx2)
# 	args$ylab=""
# 	args$xlab=""
# 	args$par.settings=list(axis.line = list(col = "gray75"))
#
# 	# pearson
# 	args$panel=function(x1, y1, ...){
# 		panel.grid(col="gray95")
# 		panel.segments(df0$x1, df0$y1, df0$x2, df0$resprs, ...)
# 		panel.xyplot(df0$x1, df0$resprs, ...)
# 		}
# 	args$main="Pearson residuals"
# 	p4 <- do.call("xyplot", args)
#
# 	# standardized
# 	args$panel=function(x1, y1, ...){
# 		panel.grid(col="gray95")
# 		panel.segments(df0$x1, df0$y1, df0$x2, df0$resstd, ...)
# 		panel.xyplot(df0$x1, df0$resstd, ...)
# 		}
# 	args$main="Standardized residuals"
# 	p5 <- do.call("xyplot", args)
#
# 	# deviances
# 	args$panel=function(x1, y1, ...){
# 		panel.grid(col="gray95")
# 		panel.segments(df0$x1, df0$y1, df0$x2, df0$resraw, ...)
# 		panel.xyplot(df0$x1, df0$resraw, ...)
# 		}
# 	args$main="Raw residuals (deviances)"
# 	p6 <- do.call("xyplot", args)
#
# 	#----------------------------------------------------------------
# 	# build plot
# 	#----------------------------------------------------------------
#
# 	subtext <- paste("(shaded area = CI", ci*100, "%, dashed line = median, solid line = observed \n", sep="")
#
# 	if(type=="prediction"){
# 		p3$sub <- list(label=subtext, cex=1)
# 		p3$main <- list(label="Prediction error" , cex=1.5)
# 		print(p3)
# 	}
# 	if(type=="all"){
# 		grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3, as.table=FALSE, top=textGrob("Aggregated catch diagnostics \n", gp=gpar(fontface = "bold", cex=1.5)), bottom=textGrob(subtext, gp=gpar(cex=1)))
# 	}
# })

setMethod("plot", c("a4aFitCatchDiagn", "missing"), function(x, y=missing, probs=c(0.1, 0.9), type="all", cex=1, ...){

  #----------------------------------------------------------------
  # 1. SETUP & DATA EXTRACTION
  #----------------------------------------------------------------
  ci <- probs[2] - probs[1]
  probs <- c(probs[1], 0.5, probs[2])

  d1 <- as.data.frame(quantile(x$oe, probs))
  d2 <- as.data.frame(quantile(x$ee, probs))
  d3 <- as.data.frame(quantile(x$oee, probs))

  d1$iter <- factor(d1$iter, labels=c("l","m","u"))
  d2$iter <- factor(d2$iter, labels=c("l","m","u"))
  d3$iter <- factor(d3$iter, labels=c("l","m","u"))

  res_prs <- c(x$resprs)
  res_std <- c(x$resstd)
  res_raw <- c(x$resraw)

  years <- as.numeric(as.character(unique(d1$year)))
  obs <- c(x$obs)
  years_obs <- as.numeric(dimnames(x$obs)[[2]])

  # --- THE FIX: Get the exact min and max years ---
  x_lims <- range(c(years, years_obs))

  # Create nice internal ticks (e.g., every 5 years) but drop ones outside the data
  base_ticks <- pretty(x_lims)
  base_ticks <- base_ticks[base_ticks >= x_lims[1] & base_ticks <= x_lims[2]]

  # Force the axis to draw from the exact start year to the exact end year
  x_ticks <- unique(c(x_lims[1], base_ticks, x_lims[2]))

  #----------------------------------------------------------------
  # 2. AESTHETICS & HELPER FUNCTIONS
  #----------------------------------------------------------------
  col_fill <- "gray90"
  col_line <- "black"
  col_obs  <- "black"
  col_res  <- "black"

  plot_ribbon <- function(df, main_title) {
    l <- df$data[df$iter == "l"]
    m <- df$data[df$iter == "m"]
    u <- df$data[df$iter == "u"]

    panel_min <- min(c(l, obs), na.rm = TRUE) * 0.95
    panel_max <- max(c(u, obs), na.rm = TRUE) * 1.05

    # Add xlim = x_lims here
    plot(years, m, type = "n", xlim = x_lims, ylim = c(panel_min, panel_max),
         xlab = "", ylab = "", main = main_title,
         axes = FALSE, cex.main = cex*1.5, font.main = 1)

    grid(nx = NA, ny = NULL, col = "gray90", lty = 1)
    axis(1, at = x_ticks, col = "gray50", col.axis = "gray30")
    axis(2, las = 1, col = "gray50", col.axis = "gray30")

    polygon(c(years, rev(years)), c(l, rev(u)), col = col_fill, border = NA)
    lines(years, m, lty = 2, col = col_line, lwd = cex*1.75)
    lines(years_obs, obs, lty = 1, col = col_obs, lwd = cex*1.5)
    # points(years_obs, obs, pch = 16, col = col_obs, cex = 1.2)
  }

  plot_lollipop <- function(res_data, main_title) {
    panel_max <- max(abs(res_data), na.rm = TRUE) * 1.05
    panel_ylim <- c(-panel_max, panel_max)

    # Add xlim = x_lims here too
    plot(years_obs, res_data, type = "n", xlim = x_lims, ylim = panel_ylim,
         xlab = "", ylab = "", main = main_title,
         axes = FALSE, cex.main = cex*1.5, font.main = 1)

    grid(nx = NA, ny = NULL, col = "gray90", lty = 1)
    axis(1, at = x_ticks, col = "gray50", col.axis = "gray30")
    axis(1, at = x_ticks, col = "gray50", col.axis = "gray30")
    axis(2, las = 1, col = "gray50", col.axis = "gray30")

    abline(h = 0, col = "gray50", lwd = cex)
    segments(years_obs, 0, years_obs, res_data, col = col_res, lwd = cex*1.2)
    points(years_obs, res_data, pch = 16, col = col_res, cex = cex*1.2)
  }

  #----------------------------------------------------------------
  # 3. BUILD LAYOUT & PLOTS
  #----------------------------------------------------------------
  subtext <- sprintf("(shaded area = CI %s%%, dashed line = median, solid line = observed)", ci*100)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if (type == "prediction") {
    par(mar = c(5, 4, 4, 2) + 0.1, oma = c(2, 0, 2, 0))
    plot_ribbon(d3, "Prediction error")
    mtext(subtext, side = 1, line = 4, cex = cex, col = "gray30")

  } else if (type == "all") {
    par(mfcol = c(3, 2),
        mar = c(2, 4, 3, 1),
        oma = c(4, 2, 4, 2))

    plot_ribbon(d1, "Observation error")
    plot_ribbon(d2, "Estimation error")
    plot_ribbon(d3, "Prediction error")

    plot_lollipop(res_prs, "Pearson residuals")
    plot_lollipop(res_std, "Standardized residuals")
    plot_lollipop(res_raw, "Raw residuals (deviances)")

    mtext("Aggregated catch diagnostics", side = 3, outer = TRUE,
          font = 2, cex = cex*1.3, line = 1)
    mtext(subtext, side = 1, outer = TRUE,
          cex = cex*1.2, col = "gray30", line = 1)
  }
})
