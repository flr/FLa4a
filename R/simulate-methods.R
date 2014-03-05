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
  function(object, nsim = 1, iter = NULL) {
    out <- object
    out @ pars <- simulate(pars(object), nsim = nsim, iter = iter)
    
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
      iages <- rownames(out @ index[[i]])
      iyears <- colnames(out @ index[[i]])
      when <- mean(range(qmodel(pars(out))[[i]])[c("startf", "endf")])
      out @ index[[i]] <- (out @ stock.n * exp(-Zs * when))[iages, iyears] * out @ index[[i]]
    }

    out
})

#' @rdname simulate-methods
#' @aliases simulate,SCAPars-method
setMethod("simulate", signature(object = "SCAPars"),
  function(object, nsim = 1, iter = NULL) {    
    out <- object
    
    out @ stkmodel <- simulate(object @ stkmodel, nsim = nsim, iter = iter)
    out @ qmodel <- simulate(object @ qmodel, nsim = nsim, iter = iter)
    out @ vmodel <- simulate(object @ vmodel, nsim = nsim, iter = iter)

    out
  })

#' @rdname simulate-methods
#' @aliases simulate,a4aStkParams-method
setMethod("simulate", signature(object = "a4aStkParams"),
  function(object, nsim = 1, iter = NULL) {    

    # sanity checks
    if (is.null(iter)) {
      if (nsim == 1) iter <- seq(dim(object @ params)[2])
      if (nsim > 1) iter <- 1
    } else
    {
      if (nsim == 1) {
        if (any(iter > dim(object @ params)[2]) | any(iter < 0)) {
          message("supplied values of iter are not sensible... simulating from all iters")
          iter <- seq_along(dim(object @ params)[2])
        }
    } else
      {
        if (length(iter) > 1) stop("if nsim > 1 iter must be of length 1")
        if (iter > dim(object @ params)[2] | iter < 0) {
          message("supplied values of iter are not sensible... simulating from iter = 1")
          iter <- 1
        }
      }
    }
    
    # get parameter estimates
    b <- coef(object)
    
    # get parameter variance matrices
    V <- vcov(object)

    # simulate some new params from the first iteration only!
    if (dim(V)[3] == 1) {
      itervar <- rep(1, length(iter))
    } else {
      itervar = iter
    }
    parsim <-
      sapply(seq_along(iter),
        function(i) 
          t(mvrnorm(nsim, c(b[,iter[i]]), V[,,itervar[i]])))

    
    # load simpars into a SCAPars object and return
    out <- object
    
    if (nsim == 1) { 
      out @ params <- object @ params
    } else {
      out @ params <- propagate(object @ params[,iter], nsim)
    }
    out @ params[] <- c(parsim)

    ####
    # note we set the variance matrices to zero
    # since having a variance no longer makes sense...
    ####
    vcov(out) <- 0

    return(out)
})

#' @rdname simulate-methods
#' @aliases simulate,submodels-method
setMethod("simulate", signature(object = "submodels"),
  function(object, nsim = 1, iter = NULL) {
    out <- lapply(object, simulate, nsim = nsim, iter = iter)
    submodels(out)
  })

#' @rdname simulate-methods
#' @aliases simulate,submodel-method
setMethod("simulate", signature(object = "submodel"),
  function(object, nsim = 1, iter = NULL) {

    # sanity checks
    if (is.null(iter)) {
      if (nsim == 1) iter <- seq(dim(object @ params)[2])
      if (nsim > 1) iter <- 1
    } else
    {
      if (nsim == 1) {
        if (any(iter > dim(object @ params)[2]) | any(iter < 0)) {
          message("supplied values of iter are not sensible... simulating from all iters")
          iter <- seq_along(dim(object @ params)[2])
        }
    } else
      {
        if (length(iter) > 1) stop("if nsim > 1 iter must be of length 1")
        if (iter > dim(object @ params)[2] | iter < 0) {
          message("supplied values of iter are not sensible... simulating from iter = 1")
          iter <- 1
        }
      }
    }
    
    # get parameter estimates
    b <- coef(object)
    
    # get parameter variance matrices
    V <- vcov(object)

    # simulate some new params from the first iteration of V only!
    if (dim(V)[3] == 1) {
      itervar <- rep(1, length(iter))
    } else {
      itervar = iter
    }
    parsim <- 
      sapply(seq_along(iter), 
            function(i) 
              t(mvrnorm(nsim, c(b[,iter[i]]), V[,,itervar[i]])))

    # load simpars into a submodel object and return
    out <- object
    
    if (nsim == 1) { 
      out @ params <- object @ params
    } else {
      out @ params <- propagate(object @ params[,iter], nsim)
    }
    out @ params[] <- c(parsim)

    ####
    # note we set the variance matrices to zero
    # since having a variance no longer makes sense...
    ####
    vcov(out) <- 0

    return(out)
})

