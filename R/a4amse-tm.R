#example use
#w(stk0, scn, ay, vectorY, fy)   

#input: stk0=observation FLStock
#SCN: scenario
#ay: assessment Year
#vectorY: years to be filled in
#fy: last year
#output: new stk0 object with harvest(stk) changed

#w <- function(...){
#	args <- list(...)
#	method <- args$method
#	args$method <- NULL
#	# Check inputs
#	if(!is(args$stk,"FLS")) stop("stk argument must be an FLStock")
#	# dispatch
#	out <- do.call(method, args)
#	# check outputs
#	if(!is(out$snew, "FLQuant")) stop("The technical measures must return and object of class FLQuant")	
#	# return
#	out  
#}

mpa.tm <- function(stk, sqy, sel.objective, tracking){             
	sold <- snew <- yearMeans(harvest(stk)[,sqy])
	snew[] <- predict(sel.objective, x=as.numeric(dimnames(snew)[[1]]))
	v <- range(stk, "minfbar"):range(stk,"maxfbar")
	snew <- snew * mean(sold[v])/mean(snew[v])
	list(flq=snew, tracking=tracking)
}

