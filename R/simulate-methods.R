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
  function(object, nsim = 1, seed = NULL, stock=NULL) {

    out <- object
    out @ pars <- simulate(pars(object), nsim = nsim, seed = seed)

    # now get catch.n, stock.n, harvest and index
    preds <- predict(out)
    out @ harvest <- preds $ stkmodel $  harvest
    out @ stock.n <- out @ catch.n <- out @ harvest
    out @ stock.n[1,] <- preds $ stkmodel $  rec
    out @ stock.n[-1,1] <- preds $ stkmodel $ ny1[-1,]

    # plusgroup?
    dms <- dims(object)
    plusgrp <- !is.na(dms $ plusgroup) && dms $ plusgroup >= dms $ max
  
    # build stock
    Zs <- harvest(out) + m(out)
    for (a in 2:dms $ age) {
      out @ stock.n[a,-1] <- out @ stock.n[a-1, 1:(dms $ year-1)] * exp( - Zs[a-1, 1:(dms $ year-1)] )
    }
    # if plus group
    if (plusgrp) {
      for (y in 1:(dms $ year-1)) 
        out @ stock.n[a,y+1,] <- out @ stock.n[a,y+1,] + out @ stock.n[a, y,] * exp( - Zs[a, y,] )
    } 
 
    # calculate catch
    zfrac <- harvest(out) / Zs * (1 - exp(-Zs))
    out @ catch.n <- zfrac * out @ stock.n

    # work out indices
    out @ index <- preds $ qmodel
    for (i in seq(out @ index)) {
      idx <- index(out)[[i]]
      dnms <- dimnames(idx)	
      iages <- dnms$age
      iyears <- dnms$year
	  when <- mean(range(qmodel(pars(out))[[i]])[c("startf", "endf")])
	  # is it a biomass idx ?
	  bioidx <- FALSE
	  if(length(iages)==1) if(iages=="all") bioidx <- TRUE 
      # if biomass index all ages have to be accounted
      # WARNING: spagheti code
	  if(bioidx){
		if(missing(stock)){
		   warning("Can't simulate the biomass index. Please provide FLStock to get stock weights.")
		   out @ index[[i]][] <- object@index[[i]][]
		} else {
			stk <- object@stock.n*exp(-Zs * when)
			stk <- mcf(list(e1=stk, e2=stock.wt(stock)))
			stk$e2[] <- stk$e2[,,,,,1]
			stk <- quantSums(do.call("*", stk))[,iyears]
	        out @ index[[i]] <- stk * out @ index[[i]]
		}
	  # or else it's a age based index
	  } else {
	   stk <- (object @ stock.n * exp(-Zs * when))[iages, iyears]
       out @ index[[i]] <- stk * out @ index[[i]]
      }
    }

    out
})

#' @rdname simulate-methods
#' @aliases simulate,SCAPars-method
setMethod("simulate", signature(object = "SCAPars"),
  function(object, nsim = 1, seed=NULL) {    
    out <- object
    out @ stkmodel <- simulate(object @ stkmodel, nsim = nsim, seed = seed)
    out @ qmodel <- simulate(object @ qmodel, nsim = nsim, seed = seed)
    out @ vmodel <- simulate(object @ vmodel, nsim = nsim, seed = seed)
    out
  })

#' @rdname simulate-methods
#' @aliases simulate,a4aStkParams-method
setMethod("simulate", signature(object = "a4aStkParams"),
  function(object, nsim = 1, seed=NULL) {    

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

	# simulate
    if(is.null(seed)){
	    parsim <- sapply(seq_along(iter), function(i){ t(mvrnorm(1, c(b[,iter[i]]), V[,,iter[i]]))}) 
    } else {
		set.seed(seed)
	    parsim <- sapply(seq_along(iter), function(i){ t(mvrnorm(1, c(b[,iter[i]]), V[,,iter[i]]))}) 
    }

	# add to object and return
	object@params@.Data[] <- c(parsim)
	if(vitr==1) object@vcov <- object@vcov[,,1, drop=FALSE]
	return(object)

})

#' @rdname simulate-methods
#' @aliases simulate,submodels-method
setMethod("simulate", signature(object = "submodels"),
  function(object, nsim = 1, seed=NULL) {
    out <- lapply(object, simulate, nsim = nsim, seed=seed)
    submodels(out)
  })

#' @rdname simulate-methods
#' @aliases simulate,submodel-method
setMethod("simulate", signature(object = "submodel"),
  function(object, nsim = 1, seed=NULL) {

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

	# simulate
    if(is.null(seed)){
	    parsim <- sapply(seq_along(iter), function(i){ t(mvrnorm(1, c(b[,iter[i]]), V[,,iter[i]]))}) 
    } else {
		set.seed(seed)
	    parsim <- sapply(seq_along(iter), function(i){ t(mvrnorm(1, c(b[,iter[i]]), V[,,iter[i]]))}) 
    }

	# add to object and return
	object@params@.Data[] <- c(parsim)
	if(vitr==1) object@vcov <- object@vcov[,,1, drop=FALSE]
	return(object)
})

