#o <- function(...){
#	# in: stk = FLStock, idx = FLIndices
#	#	in: method = character for the wrapper function on the assessment.
#	# out: list with elements stk and idx sampled from OM
#	args <- list(...)
#	method <- args$method
#	args$method <- NULL
#	# checks 
#	if(!is(args$stk, "FLS")) stop("stk must be of class FLStock")
#	# dispatch
#	out <- do.call(method, args)
#	if(!is(out$stk, "FLS")) stop("stk must be of class FLStock")
#	if(!is(out$idx, "FLIndices")) stop("idx must be of class FLIndices")
#	out
#}

sampling.oem <- function(stk, deviances, observations, vy0, ay, tracking, oe=c("both","index","catch")){
	# needs more work to remove the index OE, for now index OE is mandatory 

	dataYears <- vy0
	assessmentYear <- ac(ay)
	# dataYears is a position vector, not the years themselves

	# carry on stock information in the observations for "short-cut" approach
	stock.n(observations$stk)[,assessmentYear] <- stock.n(stk)[,assessmentYear]	
	
	# catch.n
	if(oe %in% c("both","catch")){
		catch.n(observations$stk)[,max(dataYears)] <- catch.n(stk)[,max(dataYears)]*deviances$stk$catch.n[,max(dataYears)]
		catch(observations$stk)[,max(dataYears)] <- computeCatch(observations$stk[,max(dataYears)])
		stk0 <- observations$stk[,dataYears]
	}

	# indices
	if(oe %in% c("both","index")){
		for (idx_count in 1:length(observations$idx)){
			index(observations$idx[[idx_count]])[,max(dataYears)] <- stock.n(stk)[,max(dataYears)]*deviances$idx[[idx_count]][,max(dataYears)]
		}
		idx0 <- lapply(observations$idx, function(x) x[,dataYears])
	}

	# return
	list(stk=stk0, idx=idx0, observations=observations, tracking=tracking)
}

perfectInfo.oem <- function(stk, deviances, observations, vy0, ay, tracking){
	dataYears <- vy0
	assessmentYear <- ac(ay)
	stk0 <- stk[,dataYears]
	idx0 <- FLIndices(a=FLIndex(index=stock.n(stk)[,dataYears]*0.01))
	range(idx0[[1]])[c("startf","endf")] <- c(0,0)
	list(stk=stk0, idx=idx0, deviances, observations, tracking=tracking)
}



