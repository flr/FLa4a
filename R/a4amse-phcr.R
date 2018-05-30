# HCR parametrization
#x <- function(...){
#	args <- list(...)
#	method <- args$method
#	args$method <- NULL
#	# Check inputs
#	if(!is(args$stk,"FLS")) stop("stk argument must be an FLStock")
#	# dispatch
#	out <- do.call(method, args)
#	# check outputs
#	#if(!is(refpts, "FLPar")) stop("The HCR parametrization must return and object of class FLPar")	
#	# return
#	out  
#}

movFtrg.phcr <- function(stk, frp="f0.1", model="missing", ay, iy, hcrpars, interval, tracking){
	if(ay==iy | (ay-iy)%%interval==0){
		if(!missing(model)){
			sr0 <- fmle(as.FLSR(stk, model=model))
			hcrpars <- c(refpts(FLBRP:::brp(FLBRP(stk, sr0)))[frp,"harvest"])
		} else {
			hcrpars <- c(refpts(brp(FLBRP(stk)))[frp,"harvest"])
		}
	}
	list(hcrpars=hcrpars, tracking=tracking)	
}

