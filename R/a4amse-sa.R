#f <- function(...){
#	args <- list(...)
#	method <- args$method
#	args$method <- NULL
#	# in: stk = FLStock, idx = FLIndices
#	#	method = character for the wrapper function on the assessment.
#	# out: list with elements 'stk' updated with the assessment and 'convergence' diagnostics
#	# checks 
#	if(!is(args$stk, "FLS")) stop("stk must be of class FLStock")
#	if(!is(args$idx, "FLIndices")) stop("idx must be of class FLIndices")
#	# dispatch
#	out <- do.call(method, args)
#	if(!is(out$stk, "FLS")) stop("stk must be of class FLStock")
##	stk <- out$stk
##	tracking <- out$tracking
##	return(list(stk = stk, tracking = tracking))
#	out
#}

sca.sa <- function(stk, idx, genArgs, update=TRUE, dfm=c(0.75, 0.75), ...){
	args <- list(...)
	if(update) args$fmodel <- defaultFmod(stk, dfm=dfm)
	args$stock <- stk
	args$indices <- idx
	if(is.null(args$fit)) args$fit <- 'MP'
	tracking <- args$tracking
	args$tracking <- NULL
	fit <- do.call('sca', args)
	stk <- stk + fit
	tracking["conv.est",ac(range(stk)["maxyear"] + 1)] <- fit@fitSumm["maxgrad",]
	list(stk = stk, tracking = tracking)
}

#xsa.sa <- function(stk, idx, ...){
#	args <- list(...)
#	args$stock <- stk
#	args$indices <- idx
#	if(is.null(args$control)) args$control <- FLXSA.control()
#	tracking <- args$tracking
#	args$tracking <- NULL
#	fit <- do.call('FLXSA', args)
#	stk <- stk + fit
#	tracking["convergence",ac(range(stk)["maxyear"]+1)] <- fit@control@maxit
#	list(stk = stk, tracking = tracking)
#}

sep.sa <- function(stk, idx, genArgs, update=TRUE, dfm=c(0.75, 0.75), ...){
	args <- list(...)
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
		args$fmodel <- fmodel
	}
	args$stock <- stk
	args$indices <- idx
	if(is.null(args$fit)) args$fit <- 'MP'
	tracking <- args$tracking
	args$tracking <- NULL
	fit <- do.call('sca', args)
	stk <- stk + fit
	tracking["conv.est",ac(range(stk)["maxyear"]+1)] <- fit@fitSumm["maxgrad",]
	list(stk = stk, tracking = tracking)
}

# ma.wrapper <- function(stk, idx, xsa.control=list(), sca.control=list()){
#   
#   # XSA
#   if(is.null(xsa.control)){
#     print("Using default xsa control settings")
#     xsa.control  <- FLXSA.control()
#   }
#   # Fit XSA
#   fit <- FLXSA(stk, idx, xsa.control)
#   # convergence diagnostic (quick and dirty)
#   maxit <- fit@control@maxit
#   # Update stk0
#   stk.xsa <- stk+fit
# 
#   # SCA
#   fit <- sca(stk, idx)
#   stk.sca1 <- stk+fit  
# 
#   # SCA
#   if(!"fmod" %in% names(sca.control)){
#     print("Using default fmod sca settings")
#     fmod <- getFmod(stk, dfm=c(2/3, 1/2))
#   } else {
#     fmod <- control$fmod
#   }
#   if(!"qmod" %in% names(sca.control)){
#     print("Using default qmod sca settings")
#     qmod <- getQmod(idx)
#   } else {
#     qmod <- control$qmod
#   }
#   fit <- sca(stk, idx, fmodel=fmod, qmodel=qmod)
#   stk.sca <- stk+fit  
# 
#   # SCA separable
#   if(!"fmod" %in% names(sca.control)){
#     print("Using default fmod sca settings")
#     fmod <- getFmod(stk, dfm=c(2/3, 1/2), model='separable')
#   } else {
#     fmod <- control$fmod
#   }
#   if(!"qmod" %in% names(sca.control)){
#     print("Using default qmod sca settings")
#     qmod <- getQmod(idx)
#   } else {
#     qmod <- control$qmod
#   }
#   fit <- sca(stk, idx, fmodel=fmod, qmodel=qmod)
#   stk.sep <- stk+fit  
#   
#   catch.n(stk) <- (catch.n(stk.sca1)+catch.n(stk.sca)+catch.n(stk.xsa)+catch.n(stk.sep))/4
#   stock.n(stk) <- (stock.n(stk.sca1)+stock.n(stk.sca)+stock.n(stk.xsa)+stock.n(stk.sep))/4
#   harvest(stk) <- (harvest(stk.sca1)+harvest(stk.sca)+harvest(stk.xsa)+harvest(stk.sep))/4
#   list(stk=stk, convergence=NA)
# }
# 
# 
# ivw.wrapper <- function(stk, idx, xsa.control=list(), sca.control=list()){
#   
# 
#   # SCA
#   if(!"fmod" %in% names(sca.control)){
#     fmod <- getFmod(stk, dfm=c(2/3, 1/2))
#   } else {
#     fmod <- control$fmod
#   }
#   if(!"qmod" %in% names(sca.control)){
#     qmod <- getQmod(idx)
#   } else {
#     qmod <- control$qmod
#   }
#   fit <- sca(stk, idx, fmodel=fmod, qmodel=qmod)
#   stk.sca <- stk+fit  
# 	v1 <- sum((residuals(fit, stk, idx)$catch.n)^2)
# 
#   # SCA separable
#   if(!"fmod" %in% names(sca.control)){
#     print("Using default fmod sca settings")
#     fmod <- getFmod(stk, dfm=c(2/3, 1/2), model='separable')
#   } else {
#     fmod <- control$fmod
#   }
#   if(!"qmod" %in% names(sca.control)){
#     print("Using default qmod sca settings")
#     qmod <- getQmod(idx)
#   } else {
#     qmod <- control$qmod
#   }
#   fit <- sca(stk, idx, fmodel=fmod, qmodel=qmod)
#   stk.sep <- stk+fit  
# 	v2 <- sum((residuals(fit, stk, idx)$catch.n)^2)
#   
#   w1 <- (1/v1)/(1/v1+1/v2)
#   
#   catch.n(stk) <- w1*catch.n(stk.sca)+(1-w1)*catch.n(stk.sep)
#   stock.n(stk) <- w1*stock.n(stk.sca)+(1-w1)*stock.n(stk.sep)
#   harvest(stk) <- w1*harvest(stk.sca)+(1-w1)*harvest(stk.sep)
#   list(stk=stk, convergence=NA)
# }
# 
# 
# eqw.wrapper <- function(stk, idx, xsa.control=list(), sca.control=list()){
# 
#   # SCA
#   if(!"fmod" %in% names(sca.control)){
#     fmod <- getFmod(stk, dfm=c(2/3, 1/2))
#   } else {
#     fmod <- control$fmod
#   }
#   if(!"qmod" %in% names(sca.control)){
#     qmod <- getQmod(idx)
#   } else {
#     qmod <- control$qmod
#   }
#   fit <- sca(stk, idx, fmodel=fmod, qmodel=qmod)
#   stk.sca <- stk+fit  
# 
#   # SCA separable
#   if(!"fmod" %in% names(sca.control)){
#     print("Using default fmod sca settings")
#     fmod <- getFmod(stk, dfm=c(2/3, 1/2), model='separable')
#   } else {
#     fmod <- control$fmod
#   }
#   if(!"qmod" %in% names(sca.control)){
#     print("Using default qmod sca settings")
#     qmod <- getQmod(idx)
#   } else {
#     qmod <- control$qmod
#   }
#   fit <- sca(stk, idx, fmodel=fmod, qmodel=qmod)
#   stk.sep <- stk+fit  
# 
#   w1 <- 0.5
#   
#   catch.n(stk) <- w1*catch.n(stk.sca)+(1-w1)*catch.n(stk.sep)
#   stock.n(stk) <- w1*stock.n(stk.sca)+(1-w1)*stock.n(stk.sep)
#   harvest(stk) <- w1*harvest(stk.sca)+(1-w1)*harvest(stk.sep)
#   list(stk=stk, convergence=NA)
# }
