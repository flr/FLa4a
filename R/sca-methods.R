#' @title Statistical catch-at-age method
#' @name sca
#' @docType methods
#' @rdname sca
#' @description Statistical catch-at-age method of the a4a stock assessment framework.
#' @details [REQUIRES REVISION] This method is the advanced method for stock assessment, it gives the user access to a set of arguments that the \code{sca} method doesn't. In particular, the default for the \code{fit} argument is 'assessment'. For detailed information about using the \code{sca} read the vignette 'The a4a Stock Assessment Modelling Framework' (\code{vignette('sca')}).
#' @param stock an \code{FLStock} object containing catch and stock information
#' @param indices an \code{FLIndices} object containing survey indices
#' @param fmodel a formula object depicting the model for log fishing mortality at age
#' @param qmodel a list of formula objects depicting the models for log survey catchability at age
#' @param srmodel a formula object depicting the model for log recruitment
#' @param n1model a formula object depicting the model for the population in the first year of the time series
#' @param vmodel a list of formula objects depicting the model for the variance of fishing mortality and the indices
#' @param covar a list with covariates to be used by the submodels. The formula must have an element with the same name as the list element.
#' @param wkdir used to set a working directory for the admb optimiser; if wkdir is set, all admb files are saved to this folder, otherwise they are deleted.
#' @param verbose if true, admb fitting information is printed to the screen.
#' @param fit character with type of fit: 'MP' or 'assessment'; the former does not require the hessian to be computed, while the latter does.
#' @param center, logical defining if the data should be centered before fitting.
#' @param mcmc an \code{SCAMCMC} object with the arguments to run MCMC
#' @template dots
#' @return an \code{a4aFit} object if fit is "MP" or an \code{a4aFitSA} object if fit is "assessment"
#' @aliases sca sca-methods
#' @examples
#' data(ple4)
#' data(ple4.index)
#'
#' # fishing mortality by age and year (separable) AND catchability at age without year trend
#' fmodel <- ~factor(age) + factor(year)
#' qmodel <- list(~factor(age))
#' fit1 <-  sca(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#'
#' # fishing mortality as a smoother by age and year (but still separable) AND
#' # catchability at age without year trend
#' fmodel <- ~ s(age, k=4) + s(year, k=10)
#' qmodel <- list(~factor(age))
#' fit2 <-  sca(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#'
#' # fishing mortality as a smoother by age and year (but still separable) AND
#' # catchability as a smoother by age without year trend
#' fmodel <- ~ s(age, k=4) + s(year, k=10)
#' qmodel <- list(~s(age, k=4))
#' fit3 <-  sca(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#'
#' # fishing mortality as a smoother by age and year (but still separable) AND
#' # catchability as a smoother by age with year trend
#' fmodel <- ~ s(age, k=4) + s(year, k=10)
#' qmodel <- list(~s(age, k=4) + year)
#' fit4 <-  sca(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#'
#' # It's a statistical model
#' BIC(fit1, fit2, fit3, fit4)
#'
#' # fishing mortality as a smoother by age and year with interactions (i.e. non-separable) AND
#' # catchability as a smoother by age without year trend
#' fmodel <- ~ te(age, year, k=c(4, 10))
#' qmodel <- list(~s(age, k=4))
#' fit5 <-  sca(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#'
#' # fit3 + smoother in recruitment
#' fmodel <- ~ s(age, k=4) + s(year, k=20)
#' qmodel <- list(~s(age, k=4))
#' rmodel <- ~s(year, k=20)
#' fit6 <-  sca(fmodel=fmodel, qmodel=qmodel, srmodel=rmodel, ple4, FLIndices(ple4.index))
#'
#' # fit3 + bevholt
#' rmodel <- ~ bevholt(CV=0.05)
#' fit7 <-  sca(fmodel=fmodel, qmodel=qmodel, srmodel=rmodel, ple4, FLIndices(ple4.index))
setGeneric("sca", function(stock, indices, ...) standardGeneric("sca"))

#' @rdname sca
setMethod("sca", signature("FLStock", "FLIndex"),
  function(stock, indices, ...) {
    sca(stock, FLIndices(IND=indices), ...)
  }
)

#' @rdname sca
setMethod("sca", signature("FLStock", "FLIndices"),
	function(stock, indices, fmodel = missing, qmodel = missing, srmodel = missing, n1model = missing, vmodel = missing, covar = missing, wkdir = missing, verbose = FALSE, fit = "assessment", center = TRUE, mcmc = missing) {

  #-----------------------------------------------------------------
  # get fit type

  fit <- match.arg(fit, c("MP", "assessment", "MCMC"))

  #-----------------------------------------------------------------
  # set models if missing

  if(missing(fmodel)) fmodel <- defaultFmod(stock)
  if(missing(qmodel)) qmodel <- defaultQmod(indices)
  if(missing(n1model)) n1model <- defaultN1mod(stock)
  if(missing(vmodel)) vmodel <- defaultVmod(stock, indices)
  if(missing(srmodel)) srmodel <- defaultSRmod(stock)

  #-----------------------------------------------------------------
  # now to deal with iterations ...

  # create a df for dimension information:
  dms <- do.call(rbind.data.frame, c(list(catch = c(dims(stock), startf = NA, endf = NA)), lapply(indices, dims)))

  # average stock over seasons
  # NOTE: Do we have a warning msg about this ? YES !
  stock <- collapseSeasons(stock)

  # only allow 1 season for surveys
  if (any(dms$season[-1] > 1)) stop("only one season per survey - please split into seperate surveys.")

  # now do a fit for each combination of unit, area and iter...
  # if fit = MP then we return an a4aFit with the same dimensions as stock
  # if fit = assessment then we return a4aFitSA with same dimensions as stock....  \TODO only true with iters so far

  grid <- do.call(expand.grid, c(dimnames(catch.n(stock))[c(3,5)], list(iter = 1:max(dms$iter))))
  #if (!identical(sort(unique(dms$iter)), sort(unique(c(1L, max(dms$iter))))))
  if(length(unique(dms$iter[dms$iter>1]))>1)
  	stop("incosistent number of iterations in stock and indices")
  it <- max(dms$iter)
  if(fit=="MCMC" & it>1) stop("You can not run MCMC with iters on your data objects")
  if(fit=="MCMC" & it==1) it <- getN(mcmc)

  # set up objects
  # stk
  dms <- dimnames(stock.n(stock))
  dms$iter <- 1:it
  ini <- FLQuant(NA, dimnames=dms)
  out <- if (fit %in% c("MP", "sim")) a4aFit() else a4aFitSA()
  out@desc <- desc(stock)
  out@name <- name(stock)
  out@range <- range(stock)
  out@call <- match.call()
  out@harvest <- ini
  out@stock.n <- ini
  out@catch.n <- ini
  # idx
  ini <- lapply(indices, function(x){
  	dms <- dimnames(index(x))
  	dms$iter <- 1:it
	FLQuant(NA, dimnames=dms)
  })
  out@index <- FLQuants(ini)

  if (fit %in% c("assessment", "MCMC")) {
    out@pars@stkmodel@fMod <- fmodel
    out@pars@stkmodel@n1Mod <- n1model
    out@pars@stkmodel@srMod <- srmodel
    # and the same for indices
  }

  time.used <- matrix(NA, nrow = 4, ncol = nrow(grid))
  ifit <- if (fit == "sim") "assessment" else fit
  # note this niters are not the same as it. niters will come from data objects with iters to loop and
  # fit the model several times, while it is used to build the output objects, which in the case of MCMC
  # also need FLQuants with iterations.
  niters <- nrow(grid)
  for (i in seq(niters)) {
    istock <- stock[,, grid$unit[i], grid$area[i]]
	  istock <- iter(istock, min(grid$iter[i], dims(stock)$iter))

	  # check: do we need indices to have matching units, areas?
    #iindices <- lapply(indices, function(x) x[,, grid$unit[i], grid$area[i], , min(grid$iter[i], dims(x)$iter)])
    iindices <-
      lapply(indices, function(x) {
        idx <- x[,, grid$unit[i], grid$area[i]]
		    iter(idx, min(grid$iter[i], dims(x)$iter))
    })
    iindices <- FLIndices(iindices)

    # check: do we need indices to have matching units, areas?
    if (!missing(covar) & !missing(wkdir)) {
      icovar <- lapply(covar, function(x) x[,, grid$unit[i], grid$area[i], , min(grid$iter[i], dims(x)$iter)])
	    outi <- a4aInternal(fmodel = fmodel, qmodel = qmodel, srmodel = srmodel, n1model = n1model, vmodel = vmodel, stock = istock, indices = iindices, covar = icovar, wkdir = wkdir, verbose = verbose, fit = ifit, center = center, mcmc=mcmc)
    } else if(!missing(covar) & missing(wkdir)){
      icovar <- lapply(covar, function(x) x[,, grid$unit[i], grid$area[i], , min(grid$iter[i], dims(x)$iter)])
	    outi <- a4aInternal(fmodel = fmodel, qmodel = qmodel, srmodel = srmodel, n1model = n1model, vmodel = vmodel, stock = istock, indices = iindices, covar = icovar, verbose = verbose, fit = ifit, center = center, mcmc=mcmc)
    } else if(missing(covar) & !missing(wkdir)){
	    outi <- a4aInternal(fmodel = fmodel, qmodel = qmodel, srmodel = srmodel, n1model = n1model, vmodel = vmodel, stock = istock, indices = iindices, wkdir=wkdir, verbose = verbose, fit = ifit, center = center, mcmc=mcmc)
	  } else {
	    outi <- a4aInternal(fmodel = fmodel, qmodel = qmodel, srmodel = srmodel, n1model = n1model, vmodel = vmodel, stock = istock, indices = iindices, verbose = verbose, fit = ifit, center = center, mcmc=mcmc)
	  }
    if (i == 1) {
      tmpSumm <- outi@fitSumm
      out@fitSumm <- array(0, c(dim(tmpSumm), niters), c(dimnames(tmpSumm), list(iters = 1:niters)))
    }
    out@fitSumm[,i] <- outi@fitSumm
    if (fit == "MP") {
      # copy results
      out@harvest[,, grid$unit[i], grid$area[i], , grid$iter[i]] <- harvest(outi)
      out@stock.n[,, grid$unit[i], grid$area[i], , grid$iter[i]] <- stock.n(outi)
      out@catch.n[,, grid$unit[i], grid$area[i], , grid$iter[i]] <- catch.n(outi)
      # add indices
      for (j in 1:length(iindices)) {
        out@index[[j]][,, grid$unit[i], grid$area[i], , grid$iter[i]] <- index(outi)[[j]]
      }
    }

    if (fit == "assessment") {
      # store everything in a a4aFitSA object
      out@harvest[,, grid$unit[i], grid$area[i], , grid$iter[i]] <- harvest(outi)
      out@stock.n[,, grid$unit[i], grid$area[i], , grid$iter[i]] <- stock.n(outi)
      out@catch.n[,, grid$unit[i], grid$area[i], , grid$iter[i]] <- catch.n(outi)
      # add indices
      for (j in 1:length(iindices)) {
        out@index[[j]][,, grid$unit[i], grid$area[i], , grid$iter[i]] <- index(outi)[[j]]
      }

      # fill up models
      if (i == 1) {
        # on i == 1 do the initial propagation
        out@pars <- propagate(outi@pars, niters)
        for (j in seq_along(indices)) {
          # add stock centering to link qmodel back to stock size
          out@pars@qmodel[[j]]@centering[,1] <- outi@pars@qmodel[[j]]@centering - outi@pars@stkmodel@centering
        }
      } else {
        # fill in 2nd, 3rd iterations etc.
        # now the a4aFitSA bits
        out@pars@stkmodel@centering[,i] <- outi@pars@stkmodel@centering
        out@pars@stkmodel@coefficients[,i]   <- outi@pars@stkmodel@coefficients
        out@pars@stkmodel@vcov[,,i]    <- outi@pars@stkmodel@vcov
        out@pars@stkmodel@m[,,,,,i]    <- outi@pars@stkmodel@m
        out@pars@stkmodel@wt[,,,,,i]    <- outi@pars@stkmodel@wt
        out@pars@stkmodel@mat[,,,,,i]    <- outi@pars@stkmodel@mat
        # qmodel
        for (j in seq_along(indices)) {
          # add stock centering to link qmodel back to stock size
          out@pars@qmodel[[j]]@centering[,i] <- outi@pars@qmodel[[j]]@centering - outi@pars@stkmodel@centering
          out@pars@qmodel[[j]]@coefficients[,i] <- outi@pars@qmodel[[j]]@coefficients
          out@pars@qmodel[[j]]@vcov[,,i]  <- outi@pars@qmodel[[j]]@vcov
        }
        for (j in seq_along(outi@pars@qmodel@corBlocks)) {
          out@pars@qmodel@corBlocks[[j]][,,i]  <- outi@pars@qmodel@corBlocks[[j]]
        }

        # vmodel
        for (j in seq_along(out@pars@vmodel)) {
          out@pars@vmodel[[j]]@centering[,i] <- outi@pars@vmodel[[j]]@centering
          out@pars@vmodel[[j]]@coefficients[,i] <- outi@pars@vmodel[[j]]@coefficients
          out@pars@vmodel[[j]]@vcov[,,i]   <- outi@pars@vmodel[[j]]@vcov
        }
        for (j in seq_along(outi@pars@vmodel@corBlocks)) {
          out@pars@vmodel@corBlocks[[j]][,,i]  <- outi@pars@vmodel@corBlocks[[j]]
        }
      }

    }

    if (fit == "MCMC") {
      out <- a4aFitMCMC(out, mcmc=mcmc)
      # store everything
      out@harvest[,, grid$unit[i], grid$area[i], , ] <- harvest(outi)
      out@stock.n[,, grid$unit[i], grid$area[i], , ] <- stock.n(outi)
      out@catch.n[,, grid$unit[i], grid$area[i], , ] <- catch.n(outi)
      # add indices
      for (j in 1:length(iindices)) {
        out@index[[j]][,, grid$unit[i], grid$area[i], , ] <- index(outi)[[j]]
        units(out@index[[j]]) <- units(indices[[j]]@index)
      }

      # fill up models
      out@pars@stkmodel   <- outi@pars@stkmodel
	  # CHECK: NOT SURE WE NEED PROPAGATE HERE
      out@pars@stkmodel@m         <- propagate(outi@pars@stkmodel@m, it)
      out@pars@stkmodel@wt        <- propagate(outi@pars@stkmodel@wt, it)
      out@pars@stkmodel@units     <- units(catch.n(stock))
      # qmodel
      out@pars@qmodel             <- outi@pars@qmodel
      for (j in seq_along(indices)) {
          # add stock centering to link qmodel back to stock size
          out@pars@qmodel[[j]]@centering <- outi@pars@qmodel[[j]]@centering - outi@pars@stkmodel@centering
      }
      # vmodel
      out@pars@vmodel               <- outi@pars@vmodel
	}
    # keep timing info
    time.used[,i] <- outi@clock
  }

  units(out@harvest) <- "f"
  units(out@catch.n) <- units(stock@catch.n)
  units(out@stock.n) <- units(stock@catch.n)

	# tag biomass indices with attribute
	for(i in 1:length(indices)){
		attr(out@index[[i]], "FLIndexBiomass") <- is(indices[[i]], "FLIndexBiomass")
		attr(out@index[[i]], "range") <- range(indices[[i]])
	}
  # add in combined timings
  out@clock <- outi@clock # to get names
  out@clock[] <- rowSums(time.used)

  # return out
  out
})
