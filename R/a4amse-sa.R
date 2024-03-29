

#' Call sca inside the mp function
#'
#' This function provides an interface to sca() to be used inside the mp()
#' function of the mse package.
#'
#' @param stk The FLStock input object.
#' @param idx The FLIndices input object.
#' @param genArgs The mse arguments used by mp().
#' @param update Should the fmodel be updated with the default?
#' @param dfm data points fraction to be used to set the spline ks.
#' @param ... Any other arguments to sca()
#'
#' @return A list containing the estimated stock (stk, of class FLStock), and the tracking FLQuant, including convergence flags.
#'
#' @name sca.sa
#' @rdname sca.sa
#' @aliases sca.sa
#' @keywords classes

sca.sa <- function(stk, idx, update=TRUE, dfm=c(0.75, 0.75), args, tracking, ...){
	
  args0 <- list(...)
	
  if(update)
    args0$fmodel <- defaultFmod(stk, dfm=dfm)

	stk <- replaceZeros(stk)
	idx <- replaceZeros(idx)
	
  args0$stock <- stk
	args0$indices <- idx
	
  if(is.null(args0$fit))
    args0$fit <- 'MP'

  fit <- do.call('sca', args0)
	stk <- stk + fit

	track(tracking, "conv.est", ac(range(stk)["maxyear"] + 1)) <- fit@fitSumm["maxgrad",]

	list(stk = stk, tracking = tracking)
}

#' Call a separable SA inside the mp function
#'
#' This function provides an interface to a call to a separable model based on
#' sca() to be used inside the mp() function of the mse package.
#'
#' @param stk The FLStock input object.
#' @param idx The FLIndices input object.
#' @param genArgs The mse arguments used by mp().
#' @param update Should the fmodel be updated with the default?
#' @param dfm data points fraction to be used to set the spline ks.
#' @param ... Any other arguments to sca()
#'
#' @return A list containing the estimated stock (stk, of class FLStock), and the tracking FLQuant, including convergence flags.
#'
#' @name sep.sa
#' @rdname sep.sa
#' @aliases sep.sa
#' @keywords classes

sep.sa <- function(stk, idx, args, update=TRUE, dfm=c(0.75, 0.75), ...){
	args0 <- list(...)
	# set model
	if(update){
		dis <- dims(stk)
		KY=floor(dfm[1] * dis$year)
		KA=ceiling(dfm[2] * dis$age)
		if (KA >= 3) {
			KA <- min(max(3, KA), 6)
			KB <- min(max(3, KA), 10)
			fmodel <- formula(paste("~s(year, k =", KY,") + s(age, k=", KB, ")"))
		} else {
			fmodel <- formula(paste("~ age + s(year, k = ", KY,")"))
		}
		args0$fmodel <- fmodel
	}
	
	# replace zeros
	stk <- replaceZeros(stk)
	idx <- replaceZeros(idx)

	# add to args	
	args0$stock <- stk
	args0$indices <- idx
	if(is.null(args0$fit)) args0$fit <- 'MP'
	tracking <- args0$tracking
	args0$tracking <- NULL
	fit <- do.call('sca', args0)
	stk <- stk + fit
	
  track(tracking, "conv.est", ac(range(stk)["maxyear"] + 1)) <- fit@fitSumm["maxgrad",]
	
	list(stk = stk, tracking = tracking)
}
