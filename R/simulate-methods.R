#==================================================================== 
#    simulate  methods
#==================================================================== 

#' @title Simulation methods for SCA
#' @description Simulation methods for a4a stock assessment fits.
#' @name simulate
#' @rdname simulate-methods
#' @aliases simulate simulate-methods
#' @template Example-simulate
setGeneric("simulate", useAsDefault = simulate)

#' @rdname simulate-methods
#' @aliases simulate,a4aFitSA-method
setMethod("simulate", signature(object = "a4aFitSA"),
  function(object, nsim = 1, seed = NULL, empirical=FALSE) {
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
			if(is.na(rng["min"])) iages <- dimnames(stkn)[[1]] else iages <- rng["min"]:rng["max"]
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
#' @aliases simulate,SCAPars-method
setMethod("simulate", signature(object = "SCAPars"),
  function(object, nsim = 1, seed=NULL, empirical=TRUE) {    
    out <- object
    out @ stkmodel <- simulate(object @ stkmodel, nsim = nsim, seed = seed, empirical=empirical)
    out @ qmodel <- simulate(object @ qmodel, nsim = nsim, seed = seed, empirical=empirical)
    out @ vmodel <- simulate(object @ vmodel, nsim = nsim, seed = seed, empirical=empirical)
    out
  })

#' @rdname simulate-methods
#' @aliases simulate,a4aStkParams-method
setMethod("simulate", signature(object = "a4aStkParams"),
  function(object, nsim = 1, seed=NULL, empirical=TRUE) {    

    # sanity checks
browser()
	# iters and objects
    pitr <- dims(object@params)$iter
	vitr <- dim(object@vcov)[3]

	if(nsim > 1 | pitr > 1 | vitr > 1){
		mitr <- max(c(nsim, pitr, vitr))
		iter <- seq(mitr)
		# need to check if all are 1 or n
		if(pitr==1){
			object@params <- propagate(object@params, mitr) 
		}
		if(vitr==1){
			object@vcov <- object@vcov[,,rep(1,mitr), drop=FALSE] 
			dimnames(object@vcov)$iters <- iter
		}
	} else {
		iter <- 1
	}

    # get parameter estimates
    b <- coef(object)
    
    # get parameter variance matrices
    V <- vcov(object)

#	# code to deal with bug in mvrnorm when empirical = T
#	if(empirical & nsim < length(b)){
#		nsim0 <- nsim
#		nsim <- length(b)
#	}

	# simulate
    if(is.null(seed)){
	    parsim <- sapply(seq_along(iter), function(i){ t(mvrnorm(1, c(b[,iter[i]]), V[,,iter[i]], empirical=empirical))}) 
    } else {
		set.seed(seed)
	    parsim <- sapply(seq_along(iter), function(i){ t(mvrnorm(1, c(b[,iter[i]]), V[,,iter[i]], empirical=empirical))}) 
    }

#	if(empirical & nsim0 < length(b)){
#		nsim <- nsim0
#		parsim <- parsim[1:nsim,]
#	}

	# add to object and return
	object@params@.Data[] <- c(parsim)
	if(vitr==1) object@vcov <- object@vcov[,,1, drop=FALSE]
	return(object)

})

#' @rdname simulate-methods
#' @aliases simulate,submodels-method
setMethod("simulate", signature(object = "submodels"),
  function(object, nsim = 1, seed=NULL, empirical=TRUE) {
    out <- lapply(object, simulate, nsim = nsim, seed=seed, empirical=empirical)
    submodels(out)
  })

#' @rdname simulate-methods
#' @aliases simulate,submodel-method
setMethod("simulate", signature(object = "submodel"),
  function(object, nsim = 1, seed=NULL, empirical=TRUE) {

    # sanity checks

	# iters and objects
    pitr <- dims(object@params)$iter
	vitr <- dim(object@vcov)[3]

	if(nsim > 1 | pitr > 1 | vitr > 1){
		mitr <- max(c(nsim, pitr, vitr))
		iter <- seq(mitr)
		# need to check if all are 1 or n
		if(pitr==1){
			object@params <- propagate(object@params, mitr) 
		}
		if(vitr==1){
			object@vcov <- object@vcov[,,rep(1,mitr), drop=FALSE] 
			dimnames(object@vcov)$iters <- iter
		}
	} else {
		iter <- 1
	}

    # get parameter estimates
    b <- coef(object)
    
    # get parameter variance matrices
    V <- vcov(object)

#	# code to deal with bug in mvrnorm when empirical = T
#	if(empirical & nsim < length(b)){
#		nsim0 <- nsim
#		nsim <- length(b)
#	}

	# simulate
    if(is.null(seed)){
	    parsim <- sapply(seq_along(iter), function(i){ t(mvrnorm(1, c(b[,iter[i]]), V[,,iter[i]], empirical=empirical))}) 
    } else {
		set.seed(seed)
	    parsim <- sapply(seq_along(iter), function(i){ t(mvrnorm(1, c(b[,iter[i]]), V[,,iter[i]], empirical=empirical))}) 
    }

#	if(empirical & nsim0 < length(b)){
#		nsim <- nsim0
#		parsim <- parsim[1:nsim,]
#	}

	# add to object and return
	object@params@.Data[] <- c(parsim)
	if(vitr==1) object@vcov <- object@vcov[,,1, drop=FALSE]
	return(object)
})

