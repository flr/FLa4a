#==================================================================== 
#    simulate  methods
#==================================================================== 

#' @title Simulation methods for SCA
#' @description Simulation methods for a4a stock assessment fits.
#' @name simulate
#' @rdname simulate-methods
#' @examples
#' data(ple4)
#' data(ple4.index)
#' fmodel <- ~factor(age) + factor(year)
#' qmodel <- list(~factor(age))
#' fit1 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#' fit1
#' summary(fit1)
#' stock.n(fit1)

setGeneric("simulate", useAsDefault = stats::simulate)

#' @rdname simulate-methods
#' @template object
#' @param nsim number of iterations
#' @param seed \code{numeric} with random number seed 
#' @param empirical logical, shall the empirical method in MASS be used
#' @template dots
setMethod("simulate", signature(object = "a4aFitSA"),
  function(object, nsim = 1, seed = NULL, empirical=TRUE) {
    out <- object
    out @ pars <- simulate(pars(object), nsim = nsim, seed = seed, empirical=empirical)

    # now get catch.n, stock.n, harvest and index
    preds <- predict(out)
    out @ harvest <- preds $ stkmodel $  harvest
    out @ stock.n <- out @ catch.n <- out @ harvest
    out @ stock.n[1,] <- preds $ stkmodel $  rec
    out @ stock.n[-1,1] <- preds $ stkmodel $ ny1[-1,]

    # plusgroup?
    dms <- dims(object)
    plusgrp <- !is.na(dms $ plusgroup) && dms $ plusgroup >= dms $ max

	# fill stock.n (waste space save time)
	stkn <- stock.n(out)
    Zs <- harvest(out) + m(out)
    for (a in 2:dms $ age) {
      stkn[a,-1] <- stkn[a-1, 1:(dms $ year-1)] * exp( - Zs[a-1, 1:(dms $ year-1)] )
    }
    # if plus group
    if (plusgrp) {
      for (y in 1:(dms $ year-1)) 
        stkn[a,y+1,] <- stkn[a,y+1,] + stkn[a, y,] * exp( - Zs[a, y,] )
    } 

	out@stock.n <- stkn
 
    # calculate catch
    zfrac <- harvest(out) / Zs * (1 - exp(-Zs))
    out @ catch.n <- zfrac * stkn

    # work out indices
    out @ index <- idxs <- preds $ qmodel

    for (i in seq(idxs)) {
    	idx <- idxs[[i]]
    	dnms <- dimnames(idx)	
    	iages <- dnms$age
    	iyears <- dnms$year
		when <- mean(range(qmodel(pars(out))[[i]])[c("startf", "endf")])
		# is it a biomass idx ?
		bioidx <- FALSE
		if(attr(index(object)[[i]], "FLIndexBiomass")) bioidx <- TRUE 
     	# if biomass index use ages in range or all ages have to be accounted
     	# WARNING: spagheti code
	 	if(bioidx){
#		if(missing(stock)){
#		    warning("Can't simulate the biomass index. Please provide FLStock to get stock weights.")
#		    out @ index[[i]][] <- object@index[[i]][]
#	    } else {
			rng <- attr(index(object)[[i]], "range")
			if(is.na(rng["min"])) iages <- dimnames(stkn)[[1]] else iages <- ac(rng["min"]:rng["max"])
			stk <- stkn[iages]*exp(-Zs[iages] * when)
			stk <- mcf(list(e1=stk, e2=wt(out)[iages]))
			stk$e2[] <- stk$e2[,,,,,1]
			stk <- quantSums(do.call("*", stk))[,iyears]
			out @ index[[i]] <- stk * out @ index[[i]]
			attr(out@index[[i]], "FLIndexBiomass") <- attr(object@index[[i]], "FLIndexBiomass")
			attr(out@index[[i]], "range") <- attr(object@index[[i]], "range")

#	    }
	    # or else it's a age based index
		} else {
			stk <- (stkn * exp(-Zs * when))[iages, iyears]
			out @ index[[i]] <- stk * out @ index[[i]]
			attr(out@index[[i]], "FLIndexBiomass") <- attr(object@index[[i]], "FLIndexBiomass")
			attr(out@index[[i]], "range") <- attr(object@index[[i]], "range")
		}
    }
    out
})

#' @rdname simulate-methods
setMethod("simulate", signature(object = "SCAPars"),
  function(object, nsim = 1, seed=NULL, empirical=TRUE) {    
    out <- object
    out @ stkmodel <- simulate(object @ stkmodel, nsim = nsim, seed = seed, empirical=empirical)
    out @ qmodel <- simulate(object @ qmodel, nsim = nsim, seed = seed, empirical=empirical)
    out @ vmodel <- simulate(object @ vmodel, nsim = nsim, seed = seed, empirical=empirical)
    out
  })

#' @rdname simulate-methods
setMethod("simulate", signature(object = "a4aStkParams"),
  function(object, nsim = 1, seed=NULL, empirical=TRUE) {    

    # get parameter estimates
    b <- coef(object)
    
    # get parameter variance matrices
    V <- vcov(object)

	# iters and objects
    pitr <- dims(b)$iter
	vitr <- dim(V)[3]
	mitr <- max(c(nsim, pitr, vitr))
	it <- seq(mitr)

    # sanity checks - must be 1 or n
	if(!((nsim==1 | nsim==mitr) & (pitr==1 | pitr==mitr) & (vitr==1 | vitr==mitr)))
		stop("The number of iters and simulations must be 1 or n.")

	# if there are iters in pars or vcov simulation must be done by iter
	if(pitr!=1 | vitr != 1){
		# expand objects is needed
		if(pitr==mitr & vitr==1){
			V <- V[,,rep(1,mitr), drop=FALSE] 
			dimnames(V)$iters <- it
		} else if(pitr==1 & vitr==mitr){
			b <- propagate(b, mitr) 
		}
		# run simulations along iters
	    if(is.null(seed)){
		    parsim <- sapply(seq_along(it), function(i){ 
		    	mvrEmpT(1, c(b[,it[i]]), V[,,it[i]], empirical=empirical)
		    }) 
	    } else {
			set.seed(seed)
		    parsim <- sapply(seq_along(it), function(i){ 
		    	mvrEmpT(1, c(b[,it[i]]), V[,,it[i]], empirical=empirical)
		    }) 
 	   }
	} else {
		if(!is.null(seed)) set.seed(seed)
		parsim <- mvrEmpT(mitr, c(b), V[,,1,drop=T], empirical=empirical)
	}

	# add to object and return
	if(mitr > pitr) object@params <- propagate(object@params, mitr)
	object@params[] <- parsim
	return(object)

})

#' @rdname simulate-methods
setMethod("simulate", signature(object = "submodels"),
  function(object, nsim = 1, seed=NULL, empirical=TRUE) {
    out <- lapply(object, simulate, nsim = nsim, seed=seed, empirical=empirical)
    submodels(out)
  })

#' @rdname simulate-methods
setMethod("simulate", signature(object = "submodel"),
  function(object, nsim = 1, seed=NULL, empirical=TRUE) {
    # get parameter estimates
    b <- coef(object)
    
    # get parameter variance matrices
    V <- vcov(object)

	# iters and objects
    pitr <- dims(b)$iter
	vitr <- dim(V)[3]
	mitr <- max(c(nsim, pitr, vitr))
	it <- seq(mitr)

    # sanity checks - must be 1 or n
	if(!((nsim==1 | nsim==mitr) & (pitr==1 | pitr==mitr) & (vitr==1 | vitr==mitr)))
		stop("The number of iters and simulations must be 1 or n.")

	# if there are iters in pars or vcov simulation must be done by iter
	if(pitr!=1 | vitr != 1){
		# expand objects is needed
		if(pitr==mitr & vitr==1){
			V <- V[,,rep(1,mitr), drop=FALSE] 
			dimnames(V)$iters <- it
		} else if(pitr==1 & vitr==mitr){
			b <- propagate(b, mitr) 
		}
		# run simulations along iters
	    if(is.null(seed)){
		    parsim <- sapply(seq_along(it), function(i){ 
		    	mvrEmpT(1, c(b[,it[i]]), V[,,it[i]], empirical=empirical)
		    }) 
	    } else {
			set.seed(seed)
		    parsim <- sapply(seq_along(it), function(i){ 
		    	mvrEmpT(1, c(b[,it[i]]), V[,,it[i]], empirical=empirical)
		    }) 
 	   }
	} else {
		if(!is.null(seed)) set.seed(seed)
		parsim <- mvrEmpT(mitr, c(b), V[,,1,drop=T], empirical=empirical)
	}

	# add to object and return
	if(mitr > pitr) object@params <- propagate(object@params, mitr)
	object@params[] <- parsim
	return(object)
})

