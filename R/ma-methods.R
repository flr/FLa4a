#' @title Model averaging (experimental)
#' @name ma
#' @rdname ma-methods
#' @description Method to average across a set of models. This is still experimental. Use with care.
#' @param object an \code{a4aFits} object with the fits to be averaged across
#' @param stock a \code{stock} object with the original data used for fitting
#' @param fun a \code{function} that will be used to extract the values for weighting; for now it must be "AIC", "BIC" or "LogLik"
#' @param nsim a \code{numeric} with the number of simulations to be drawn
#' @return an \code{FLStock} object with iterations defined by \code{nsim}
#' @aliases ma ma-methods ma,a4aFitSAs-method
#' @examples
#' data(ple4)
#' data(ple4.indices)
#' f1 <- sca(ple4, ple4.indices, fmodel=~ factor(age) + s(year, k=20), qmodel=list(~ s(age, k = 4), ~ s(age, k = 4), ~ s(age, k = 3)), fit = "assessment")
#' f2 <- sca(ple4, ple4.indices, fmodel=~ factor(age) + s(year, k=20), qmodel=list(~ s(age, k = 4)+year, ~ s(age, k = 4), ~ s(age, k = 3)), fit = "assessment")
#' stock.sim <- ma(a4aFitSAs(list(f1=f1, f2=f2)), ple4, AIC, nsim = 100)
#' stks <- FLStocks(f1=ple4+f1, f2=ple4+f2, ma=stock.sim)
#' flqs <- lapply(stks, ssb)
#' flqs <- lapply(flqs, iterMedians)
#' xyplot(data~year, groups=qname, data=flqs, type="l")
#' plot(stks)

setGeneric("ma", function(object, ...) standardGeneric("ma"))
setMethod("ma", "a4aFitSAs", function(object, stock, FUN = AIC, nsim = 1000){
  FUN <- match.fun(FUN)
  # calculate weights
  ICs <- -1 * sapply(object, FUN)
  eICs <- exp( 0.5 * (ICs - max(ICs)))
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

