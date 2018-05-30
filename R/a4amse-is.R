# Management implementation function: k()

#' Evaluate the chosen Management implementation function
#'
#' Evaluate the chosen HCR function using the current stock perception and a control.
#' For example, TAC control or effort control.
#' Returns a ctrl object for projecting the OM.
#' @param method Name of the chosen implementation function.
#' @param stk The perceived stock.
#' @param ay The current year. The management control (e.g. TAC or effort) will be set in ay+1.
#' @param EFF Effort array (only used if effort management is being used).
#' @param EFF0 Tracking array.
#' @param control The control object for the chosen management implementation function. A list of parameters.

#k <- function(...){
#  args <- list(...)
#  method <- args$method
#  args$method <- NULL
#  # Check inputs
#  if(!is(args$stk,"FLS")) stop("stk argument must be an FLStock")
#  # dispatch
#  out <- do.call(method, args)
#  # check outputs
#  if(!is(out$ctrl, "fwdControl")) stop("The HCR must return an object of class fwdControl")	
#  # return
#  out  
#}

#' TAC implementation function
#'
#' Performs a short term forecast (STF) to hit the target F in year ay+1.
#' The resulting catch in year ay+1 is the TAC, i.e. the TAC that will result in Fbar = Ftarget.
#' The STF uses geometric mean recruitment. Fbar in the intermediate years (i.e. years between last data year and ay+1) are set as the mean of the last nsqy years.
#' The control argument is a list of parameters:
#' nsqy - number of years to average over to get Fbar for STF
#' delta_tac_min - constraint on the TAC
#' delta_tac_max - constraint on the TAC
#' @param stk The perceived FLStock.
#' @param imp_control A list with the elements: nsqy, delta_tac_min, delta_tac_max
#' @param ay The year for which the target F is set, based on the SSB in year (ay - control$ssb_lag).
#' @param EFF0 The tracking array
#' @param EFF Not used by this function but may be used by the other implementation functions

tac.is <- function(stk, ctrl, ay, nsqy=3, delta_tac_max=NA, delta_tac_min=NA, tracking){
	# reference value
	if(ay==iy) refCatch <- tracking["OM.catch", ac(ay-1)] else refCatch <- tracking["Implementation", ac(ay-1)]
	# Year range of perceived stock
	yrs <- as.numeric(dimnames(stock.n(stk))$year)
	last_data_yr <- yrs[length(yrs)]
	# Status quo years
	sqy <- (last_data_yr-nsqy+1):last_data_yr
	# Get the Fbar for the intermediate years
	fsq0 <- yearMeans(fbar(stk)[,ac(sqy)])
	# Number of intermediate years (between last data year and ay+1)
	ninter_yrs <- ay - last_data_yr
	# Control object for the STF
	ctrl <- getCtrl(c(rep(fsq0, times=ninter_yrs), ctrl@trgtArray[,"val",]), "f", (last_data_yr+1):(ay+1), dim(stock.n(stk))[6])
	# Number of projection years
	nproj_yrs <- (ay+1) - last_data_yr
	stkTmp <- stf(stk, nproj_yrs, wts.nyears=nsqy)
	# Set geomean sr relationship
	gmean_rec <- c(exp(yearMeans(log(rec(stk)[,ac(sqy)]))))
	# Project!
	stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=list(model="mean", params = FLPar(gmean_rec,iter=it)))
	# Get TAC for following year that results from hitting the F in STF
	TAC <- catch(stkTmp)[,ac(ay+1)]
	# catch stabelizers
	upper_limit <- refCatch * delta_tac_max
	lower_limit <- refCatch * delta_tac_min
	TAC <- pmin(c(upper_limit), c(TAC), na.rm=TRUE)
	TAC <- pmax(c(lower_limit), c(TAC), na.rm=TRUE)
	# new control file
	ctrl <- getCtrl(c(TAC), "catch", ay+1, it)
	list(ctrl = ctrl, tracking = tracking)
}


#' effort implementation function
#'
#' @param stk The perceived FLStock.
#' @param imp_control A list with the elements: nsqy, delta_tac_min, delta_tac_max
#' @param ay The year for which the target F is set, based on the SSB in year (ay - control$ssb_lag).

eff.is <- function(stk, ctrl, ay, tracking){
	# reference value
	if(ay==iy) fay <- tracking["OM.f",ac(ay-1)] else fay <- tracking["Implementation",ac(ay-1)]*tracking["Fperc",ac(ay)]	
	# target to reach defined by HCR (in f units)
	trgt <- ctrl@trgtArray[,"val",]
	# multiplier
	mult <- trgt/fay
	# new control file, in relative terms
	ctrl <- getCtrl(mult, "f", ay+1, it, rel.year=ay)
	list(ctrl = ctrl, tracking = tracking)
}


