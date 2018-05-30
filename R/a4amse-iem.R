# in: ctrl = fwd control
#	in: method = character for the wrapper function on the assessment.
# out: ctrl = fwd control


#l <- function(...){
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

noise.iem <- function(ctrl, fun="rlnorm", mean=0, sd=0.1, multiplicative=TRUE, tracking){
  # eh? number of non NA values in target. But if there are NAs then we need to
  # know their position for the *ctrl@trgtArray later
	iem <- list(mean = mean, sd = sd, n = sum(!is.na(ctrl@trgtArray)))
	if(multiplicative){
		ctrl@trgtArray <- do.call(fun, iem)*ctrl@trgtArray
	} else {
		ctrl@trgtArray <- do.call(fun, iem) + ctrl@trgtArray
	}
	list(ctrl=ctrl, tracking=tracking)
}

