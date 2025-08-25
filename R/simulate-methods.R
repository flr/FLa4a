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
#' fit1 <-  sca(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
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
  function(object, nsim = 1, seed = NULL, empirical=TRUE, obserror=FALSE) {
    out <- object
    out @ pars <- simulate(pars(object), nsim = nsim, seed = seed, empirical=empirical)
    # now get harvest, rec, ny1 and index
    preds <- predict(out)

    #--------------------------------------------------------
    # work out stock.n and catch.n

    out @ harvest <- preds $ stkmodel $  harvest
    out @ stock.n <- out @ catch.n <- out @ harvest
    out @ stock.n[1,] <- preds $ stkmodel $  rec
    out @ stock.n[-1,1] <- preds $ stkmodel $ ny1[-1,]
    flqs <- genStknCthn(harvest(out), m(out), stock.n(out)[,1], stock.n(out)[1], plusgrp=TRUE)
    out@catch.n <- flqs$catch.n
    out@stock.n <- flqs$stock.n
    if(obserror){
      out@catch.n <- exp(log(out@catch.n) + rnorm(nsim, 0, sqrt(preds$vmodel$catch)))
    }

    #--------------------------------------------------------
    # work out indices
    out @ index <- idxs <- preds $ qmodel
    Zs <- harvest(out) + m(out)
    stkn <- out@stock.n

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
			rng <- attr(index(object)[[i]], "range")
			if(is.na(rng["min"])) iages <- dimnames(stkn)[[1]] else iages <- ac(rng["min"]:rng["max"])
			stk <- stkn[iages]*exp(-Zs[iages] * when)
			stk <- mcf(list(e1=stk, e2=wt(out)[iages]))
			stk$e2[] <- stk$e2[,,,,,1]
			stk <- quantSums(do.call("*", stk))[,iyears]
			out @ index[[i]] <- stk * out @ index[[i]]
			attr(out@index[[i]], "FLIndexBiomass") <- attr(object@index[[i]], "FLIndexBiomass")
			attr(out@index[[i]], "range") <- attr(object@index[[i]], "range")
	    # or else it's a age based index
		} else {
			stk <- (stkn * exp(-Zs * when))[iages, iyears]
			out @ index[[i]] <- stk * out @ index[[i]]
            if(obserror){
              out @ index[[i]] <- rlnorm(nsim, log(out @ index[[i]]), sqrt(preds$vmodel[[i+1]]))
            }
            attr(out@index[[i]], "FLIndexBiomass") <- attr(object@index[[i]], "FLIndexBiomass")
			attr(out@index[[i]], "range") <- attr(object@index[[i]], "range")
		}
    }

    #--------------------------------------------------------
    # spit
    units(out) <- units(object)
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
			dimnames(V)[[3]] <- it
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
	if(mitr > pitr) coef(object) <- propagate(coef(object), mitr)
	coef(object) <- parsim
	return(object)

})

#' @rdname simulate-methods
setMethod("simulate", signature(object = "submodels"),
  function(object, nsim = 1, seed = NULL, empirical = TRUE) {
    # get the joined up coefficients
    blist <- lapply(object, coef)

    # get parameter variance matrices
    V <- vcov(object)

    # iters and objects
    pitr <- sapply(blist, function(x) dims(x)$iter)
    if (!all(pitr[1] == pitr)) stop("Not valid object, all submodels should have the same iterations!")
    pitr <- unique(pitr)
    if (length(dim(V)) == 2) dim(V) <- c(dim(V), 1)
    vitr <- dim(V)[3]
    mitr <- max(c(nsim, pitr, vitr))

    # sanity checks - must be 1 or nsim
    if (!((nsim == 1 || nsim == mitr) &
          (pitr == 1 || pitr == mitr) &
          (vitr == 1 || vitr == mitr))) {
      stop("The number of iters and simulations must be 1 or nsim.")
    }

    if (!is.null(seed)) set.seed(seed)
    # get the random generator from the dist slot
    distr <- unique(sapply(object, slot, "distr"))
    if (length(distr) > 1) stop("For now all submodels must have the same distribution")
    simfunc <- match.fun(paste0("mvr", distr))

    # combine blist into a single FLPar
    if (length(blist) == 1) {
      b <- blist[[1]]
    } else {
      b <- do.call(rbind, blist)
    }

    # if there are iters in pars or vcov simulation must be done by iter
    if (pitr != 1 || vitr != 1) {
      # expand objects as needed
      if (pitr == mitr & vitr == 1) {
        V <- V[,, rep(1,mitr), drop=FALSE]
        dimnames(V)[[3]] <- seq(mitr)
      } else
        if(pitr == 1 & vitr == mitr) {
          b <- propagate(b, mitr)
        }
      # run simulations along iters
      parsim <-
        sapply(seq(mitr),
               function(i) {
                 ## empirical must be set to FALSE for one iteration
                 ## it only makes sense to use empirical=TRUE when there are more samples than dimensions
                 simfunc(1, as(b[,i], "vector"), V[,,i], empirical = FALSE)
                })
    } else
    if (nsim > 1) {
      # for mvrnorm at least samples are in rows,
      # so transform is required for FLPar
      if (nsim < dim(V)[1]) {
        parsim <- t(simfunc(nsim, as(b, "vector"), V[,,1], empirical = FALSE))
      } else {
        parsim <- t(simfunc(nsim, as(b, "vector"), V[,,1], empirical = empirical))
      }
    } else {
      ## empirical must be set to FALSE for one iteration
      ## it only makes sense to use empirical=TRUE when there are more samples than dimensions
      parsim <- simfunc(1, as(b, "vector"), V[,,1], empirical = FALSE)
      dim(parsim) <- c(length(parsim), 1)
    }

    # add to object and return
    if (mitr > pitr) {
      for (i in seq_along(object)) {
        coef(object[[i]]) <- propagate(coef(object)[[i]], mitr)
      }
    }
    npar <- c(0, sapply(object, function(x) nrow(coef(x))))
    for (i in seq_along(object)) {
      coef(object[[i]]) <- as(parsim[sum(npar[1:i]) + 1:npar[i+1],], "vector")
    }
    object
  })

#' @rdname simulate-methods
setMethod("simulate", signature(object = "submodel"),
  function(object, nsim = 1, seed = NULL, empirical = TRUE) {
    # get parameter estimates
    b <- coef(object)

    # get parameter variance matrices
    V <- vcov(object)

    # iters and objects
    pitr <- dims(b)$iter
    if (length(dim(V)) == 2) dim(V) <- c(dim(V), 1)
    vitr <- dim(V)[3]
    mitr <- max(c(nsim, pitr, vitr))

    # sanity checks - must be 1 or nsim
    if (!((nsim == 1 || nsim == mitr) &
         (pitr == 1 || pitr == mitr) &
         (vitr == 1 || vitr == mitr))) {
      stop("The number of iters and simulations must be 1 or nsim.")
    }

    if (!is.null(seed)) set.seed(seed)
    # get the random generator from the dist slot
    simfunc <- match.fun(paste0("mvr", object@distr))

    # if there are iters in pars or vcov simulation must be done by iter
    if (pitr != 1 || vitr != 1) {
      # expand objects as needed
      if (pitr == mitr & vitr == 1) {
        V <- V[,, rep(1,mitr), drop=FALSE]
        dimnames(V)[[3]] <- seq(mitr)
      } else
        if(pitr == 1 & vitr == mitr) {
          b <- propagate(b, mitr)
        }
      # run simulations along iters
      parsim <-
        sapply(seq(mitr),
               function(i) {
                 ## empirical must be set to FALSE for one iteration
                 ## it only makes sense to use empirical=TRUE when there are multiple samples
                 simfunc(1, as(b[,i], "vector"), V[,,i], empirical = FALSE)
                })
    } else
    if (nsim > 1) {
      # for mvrnorm at least samples are in rows,
      # so transform is required for FLPar
      if (nsim < dim(V)[1]) {
        parsim <- t(simfunc(nsim, as(b, "vector"), V[,,1], empirical = FALSE))
      } else {
        parsim <- t(simfunc(nsim, as(b, "vector"), V[,,1], empirical = empirical))
      }
    } else {
      ## empirical must be set to FALSE for one iteration
      ## it only makes sense to use empirical=TRUE when there are multiple samples
      parsim <- simfunc(1, as(b, "vector"), V[,,1], empirical = FALSE)
    }

    # add to object and return
    if (mitr > pitr) coef(object) <- propagate(coef(object), mitr)
    coef(object) <- as(parsim, "vector")
    object
  })

