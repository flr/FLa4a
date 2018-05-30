# Harvest Control Rule function: h()
#' Evaluate the chosen HCR function
#'
#' Evaluate the chosen HCR function using the current stock perception and a control
#' @param method Name of the chosen HCR function.
#' @param stk The perceived stock.
#' @param ay The current year. The management control (e.g. F or effort) will be set in ay+1.
#' @param EFF Effort array (if effort management is being used).
#' @param EFF0 Tracking array.
#' @param control The control object for the chosen HCR function. A list of parameters.
#h <- function(...){
#	args <- list(...)
#	method <- args$method
#	args$method <- NULL
#	# Check inputs
#	if(!is(args$stk,"FLStock")) stop("stk argument must be an FLStock")
#	# dispatch
#	ctrl <- do.call(method, args)
#	# check outputs
#	if(!is(ctrl, "fwdControl")) stop("The HCR must return and object of class fwdControl")	
#	# return
#	ctrl  
#}

#' The typical HCR used by ICES
#'
#' The typical HCR used by ICES which sets a target F based on the SSB based on 4 parameters: blim, bsafe, fmin and ftrg.
#' F increases linearly between SSB = blim and SSB = bsafe, from F = fmin to F = ftrg.
#' If:
#' B < Blim, F = Fbycatch;
#' B > trigger, F = Fmsy;
#' B > Blim & B < trigger, F linear between Fbycatch and Fmsy;
#' F = ftrg is the maximum F, F = fmin is the minimum F.
#' F is set in year ay, based on SSB in year ay - ssb_lag.
#' The control argument is a list of parameters used by the HCR.
#' @param stk The perceived FLStock.
#' @param control A list with the elements fmin, ftrg, blim, bsafe and ssb_lag, all of which are numeric.
#' @param ay The year for which the target F is set, based on the SSB in year (ay - control$ssb_lag).
ices.hcr <- function(stk, fmin, ftrg, blim, bsafe, ssb_lag=1, ay, tracking){
	# rule
  ssb <- ssb(stk)[, ac(ay-ssb_lag)]
  fout <- FLQuant(fmin, dimnames=list(iter=dimnames(ssb)$iter))
  fout[ssb >= bsafe] <- ftrg
  inbetween <- (ssb < bsafe) & (ssb > blim)
  gradient <- (ftrg - fmin) / (bsafe - blim)
  fout[inbetween] <- (ssb[inbetween] - blim) * gradient + fmin
	# create control file
	ctrl <- getCtrl(c(fout), "f", ay+1, dim(fout)[6])
	# return
	list(ctrl=ctrl, tracking=tracking)
}

#' A fixed target f
#'
#' No matter what get F= Ftarget
#' The control argument is a list of parameters used by the HCR.
#' @param stk The perceived FLStock.
#' @param control A list with the element ftrg (numeric).
fxdFtrg.hcr <- function(stk, ftrg, ay, tracking){
	# rule 
	if(!is(ftrg, "FLQuant")) ftrg <- FLQuant(ftrg, dimnames=list(iter=dimnames(stk@catch)$iter))
	# create control file
	ctrl <- getCtrl(c(ftrg), "f", ay+1, dim(ftrg)[6])
	# return
	list(ctrl=ctrl, tracking=tracking)
}


movFtrg.hcr <- function(stk, hcrpars, ay, tracking){
	# rule 
	if(!is(hcrpars, "FLQuant")) hcrpars <- FLQuant(hcrpars, dimnames=list(iter=dimnames(stk@catch)$iter))
	# create control file
	ctrl <- getCtrl(c(hcrpars), "f", ay+1, dim(hcrpars)[6])
	# return
	list(ctrl=ctrl, tracking=tracking)
}


# Test:
# library(FLCore)
# data(ple4)
# Test ICES HCR - single iter
#ple4t <- ple4
#control <- list(fmin=0.1, ftrg=1.0, blim=100000,bsafe=200000, ssb_lag=1)
#ay <- 2008
#stock.wt(ple4t) <- 1.0
#stock.n(ple4t)[,ac(ay-control$ssb_lag)] <- 0.0
#ssb <- seq(from=0, to = 300000, length=20)
#fout <- rep(NA, length(ssb))
#for (i in 1:length(ssb)){
#    stock.n(ple4t)[10,ac(ay-control$ssb_lag)] <- ssb[i]
#    #c(ssb(ple4t)[,ac(ay-control$ssb_lag)])
#    fout[i] <- ices_hcr(ple4t, control, ay)
#}
#plot(ssb, fout)

# Test ICES HCR - multiple iters
#niter <- 10
#ple4p <- propagate(ple4, niter)
#ay <- 2008
#stock.wt(ple4p)[] <- 1.0
#stock.n(ple4p)[,ac(ay-control$ssb_lag)] <- 0.0
#stock.wt(ple4p) <- stock.wt(ple4p) * rlnorm(prod(dim(stock.wt(ple4p))), mean=0, sd=0.1)
#control <- list(fmin=0.1, ftrg=1.0, blim=100000,bsafe=200000, ssb_lag=1)
#ssb <- seq(from=0, to = 300000, length=20)
#fout <- array(NA, dim=c(length(ssb), niter))
#ssbout <- array(NA, dim=c(length(ssb), niter))
#for (i in 1:length(ssb)){
#    stock.n(ple4p)[10,ac(ay-control$ssb_lag)] <- ssb[i]
#    fout[i,] <- ices_hcr(ple4p, control, ay)
#    ssbout[i,] <- c(ssb(ple4p)[,ac(ay-control$ssb_lag)])
#}
#plot(ssbout[,1], fout[,1])
#for (i in 2:niter){
#    points(ssbout[,i], fout[,i])
#}

# Test h()
#library(FLash)
#control <- list(fmin=0.1, ftrg=1.0, blim=100000,bsafe=200000, ssb_lag=1)
#niter <- 10
#ple4p <- propagate(ple4, niter)
#ay <- 2000
#stock.n(ple4p) <- stock.n(ple4p) * rlnorm(prod(dim(stock.n(ple4p))), mean=0, sd=0.1)
#EFF <- FLQuant(NA, dimnames=list(EFF="all", year=dimnames(stock.n(ple4p))$year, iter=dimnames(stock.n(ple4p))$iter))
#eff_dmns <- dimnames(EFF)
#names(eff_dmns)[1] <- "metric"
#eff_dmns[["metric"]] <- c("Fperc", "Fhcr", "Implementation", "IEM", "Hyper")
#EFF0 <- FLQuant(NA, dimnames=eff_dmns)
## from funs.R
#getCtrl <- function(values, quantity, years, it){
#	dnms <- list(iter=1:it, year=years, c("min", "val", "max"))
#	arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
#	arr0[,,"val"] <- unlist(values)
#	arr0 <- aperm(arr0, c(2,3,1))
#	ctrl <- fwdControl(data.frame(year=years, quantity=quantity, val=NA))
#	ctrl@trgtArray <- arr0
#	ctrl
#}
#test <- h(method = "ices_hcr", stk=ple4p, ay=ay, EFF0=EFF0, control=control)


