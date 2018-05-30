#--------------------------------------------------------------------
# fleet dynamics, mostly changes in f-at-age
#j <- function(...){
#	args <- list(...)
#	method <- args$method
#	args$method <- NULL
#	# Check inputs
#	if(!is(args$ctrl, "fwdControl")) stop("ctrl must be of class fwdControl")
#	# dispatch
#	out <- do.call(method, args)
#	# check outputs
#	if(!is(out$ctrl, "fwdControl")) stop("The HCR must return and object of class fwdControl")	
#	# return
#	out  
#}

hyperstability.fb <- function(ctrl, beta=1, maxF=2, alpha=maxF^(1-beta), tracking){
	# Only operates on F targets - so nothing happens to TAC
	# This function creates a control file to be later used in the fwd()
	# function where two optional relations are established between
	# fishing effort and fishing mortality
	# Beta is in this MSE either 1 for a 1:1 linear relationship between
	# F and effort, if beta = 0.7, the relation is not linear and it can
	# mimick a hyperstability scenario.
	# alpha = maxF^(1-beta) # linear meets curve at maxF
	ctrl@trgtArray[ctrl@target[,"quantity"]=="f",,] <- alpha * ctrl@trgtArray[ctrl@target[,"quantity"]=="f",,]^beta
	list(ctrl=ctrl, tracking=tracking)
	
}

