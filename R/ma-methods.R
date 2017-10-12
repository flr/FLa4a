#' @title Model averaging (experimental)
#' @name ma
#' @rdname ma-methods
#' @description Method to average across a set of models. This is still experimental. Use with care.
#' @param object an \code{a4aFits} object with the fits to be averaged across
#' @param stock a \code{stock} object with the original data used for fitting
#' @param FUN a \code{function} to compute the weights, which must return a named vector with weights. Note the weights will be normalized to sum 1 by \code{ma()}
#' @param nsim a \code{numeric} with the number of simulations to be drawn
#' @template dots
#' @return an \code{FLStock} object with iterations defined by \code{nsim}
#' @aliases ma ma-methods
#' @examples
#' data(ple4)
#' data(ple4.indices)
#' fmod <- ~ factor(age) + s(year, k=20)
#' qmod <- list(~ s(age, k = 4), ~ s(age, k = 4), ~ s(age, k = 3))
#' f1 <- sca(ple4, ple4.indices, fmodel=fmod, qmodel=qmod, fit = "assessment")
#' qmod <- list(~ s(age, k = 4)+year, ~ s(age, k = 4), ~ s(age, k = 3))
#' f2 <- sca(ple4, ple4.indices, fmodel=fmod, qmodel=qmod, fit = "assessment")
#' # AIC weighting
#' aicwt <- function(object){
#'  ICs <- -1 * sapply(object, AIC)
#'  exp( 0.5 * (ICs - max(ICs)))
#' }
#' stock.sim <- ma(a4aFitSAs(list(f1=f1, f2=f2)), ple4, aicwt, nsim = 100)
#' # equal weighting
#' eqwt <- function(object){
#'  v <- rep(1, length(object))
#'  names(v) <- names(object)
#'  v
#' }
#' stock.sim <- ma(a4aFitSAs(list(f1=f1, f2=f2)), ple4, eqwt, nsim = 100)

setGeneric("ma", function(object, ...) standardGeneric("ma"))
#' @rdname ma-methods
setMethod("ma", "a4aFitSAs", function(object, stock, FUN, nsim = 1000){
	if(sum(unlist(lapply(object, function(x){dim(x@stock.n)[6]}))>1)>0) stop("\nNot working with stock assessment uncertainty yet, sorry !")
	FUN <- match.fun(FUN)
	warning("This method is experimental, use at your own risk ! \n")
	# calculate weights
#	if(FUN %in% c("AIC", "BIC")){
#		ICs <- -1 * sapply(object, FUN)
#		eICs <- exp( 0.5 * (ICs - max(ICs)))
#	} else {
#	}
#	eICs <- sapply(object, FUN)	
	eICs <- do.call(FUN, list(object=object))	
	weights <- eICs / sum(eICs)

	wt.table <- data.frame("weight (perc)" = round(weights * 100, 3))
	rownames(wt.table) <- names(object)
  
	message("model weights are \n\t", paste(capture.output(wt.table), collapse = "\n\t"))

	# now sample from each model 1000 times and randomly select those to combine
	sim <- sample(seq_along(object), nsim, prob = weights, replace = TRUE)

	stock.sim <- propagate(stock, nsim)
	for (i in seq_along(object)) {
		if (sum(sim == i)) stock.sim[,,,,,sim == i] <- stock.sim[,,,,,sim == i] + object[[i]]
	}
	# DONE !!
	stock.sim
})

