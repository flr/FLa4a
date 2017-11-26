globalVariables(c("obs", "year", "age", "fleet"))

#' @title a4a
#' @name a4a
#' @docType methods
#' @rdname a4a
#' @template dots
#' @description This was the old stock assessment framework user's interface. It was replaced by \code{sca} or the advanced version \code{a4aSCA}.
#' @aliases a4a
a4a <- function(...){
  .Deprecated("a4aSCA",
              msg = 'The a4a function has been renamed: please use the function a4aSCA instead.\nSee help("Deprecated")')
}

#' @title Collapse seasons
#' @name collapseSeasons
#' @docType methods
#' @rdname collapseSeasons
#' @description Method to collapse seasons of \code{FLStock} objects. M and catch-at-age are summed while mean weights at age, maturity at age and mortalities before spawning are averaged.
#' @param stock an FLStock object
#' @return a FLStock object
#' @aliases collapseSeasons
collapseSeasons <- function (stock) {

  if (dims(stock)$season == 1) return (stock) # do nothing

  out <- FLStock(catch.n      = seasonSums(catch.n(stock)),
                 # straight forward averages
                 m.spwn       = seasonMeans(m.spwn(stock)),
                 harvest.spwn = seasonMeans(harvest.spwn(stock)),
                 # weighted averages
                 catch.wt     = seasonSums(catch.n(stock) * catch.wt(stock)) / seasonSums(catch.n(stock)),
                 #landings.wt  = seasonSums(landings.n(stock) * landings.wt(stock)) / seasonSums(landings.n(stock)),
                 #discards.wt  = seasonSums(discards.n(stock) * discards.wt(stock)) / seasonSums(discards.n(stock)),
                 # should be weighted means but we do not have stock.n
                 # should be:
                 #stock.wt(res) <- seasonSums(stock.n(stock) * stock.wt(stock)) / seasonSums(stock.n(stock))
                 #mat(res)      <- seasonSums(stock.n(stock) * mat(stock)) / seasonSums(stock.n(stock))
                 # but we do:
                 stock.wt    = seasonMeans(stock.wt(stock)),
                 mat         = seasonMeans(mat(stock)),
                 m           = seasonSums(m(stock))
                )
  units(harvest(out)) <- units(harvest(stock))

  message("Note: Seasonal M's are summed: i.e. we assume that the values are M * season length.")
  message("Note: we do not use stock.n to weight the seasonal stock.wt or maturity.")

  out
}

#' @title Statistical catch-at-age method
#' @name sca
#' @docType methods
#' @rdname sca
#' @description User interface to the statistical catch-at-age method of the a4a stock assessment framework.
#' @details This method is the simple method for stock assessment. Some arguments are not accessible to make it simpler for the user. The advanced method is  \code{a4aSCA}. In particular, the default for the \code{fit} argument is 'MP'. For detailed information about using the \code{sca} read the vignette 'The a4a Stock Assessment Modelling Framework' (\code{vignette('sca')}).
#' @param stock an \code{FLStock} object
#' @param indices an \code{FLIndices} object
#' @param fmodel a formula object depicting the model for log fishing mortality at age
#' @param qmodel a list of formula objects depicting the models for log survey catchability at age
#' @param srmodel a formula object depicting the model for log recruitment
#' @param fit character with type of fit: 'MP' or 'assessment'; the former does not require the hessian to be computed, while the latter does.
#' @param mcmc an \code{SCAMCMC} object with the arguments to run MCMC
#' @param useADMB if FALSE (default) use TMB if TRUE use ADMB (note MCMC cannot be used with TMB)
#' @template dots
#' @return an \code{a4aFit} or \code{a4aFitSA} object with the fit results.
#' @aliases sca sca-methods
#' @examples
#' data(ple4)
#' data(ple4.index)
#'
#' # fit using the default submodels
#' fit1 <- sca(ple4, FLIndices(ple4.index))
#' plot(ple4 + fit1)
#'
#' # see default submodels (set through automated procedure)
#' sca(ple4, FLIndices(ple4.index), fit='assessment')
#'
#' # fishing mortality by age and year (separable)
#' fit2 <- sca(ple4, FLIndices(ple4.index), fmodel=~factor(age) + factor(year))
#' plot(ple4 + fit2)
#' wireframe(data~year*age, data=harvest(fit2), zlab="F")
#'
#' # fit2 + catcability as a smoother by age without year trend
#' fmod <- ~factor(age) + factor(year)
#' qmod <- list(~s(age, k=4))
#' fit3 <- sca(ple4, FLIndices(ple4.index), fmodel=fmod, qmodel=qmod)
#' plot(ple4 + fit3)
#'
#' # fit3 + srmodel as a smoother by year
#' fmod <- ~factor(age) + factor(year)
#' qmod <- list(~s(age, k=4))
#' srmod <- ~s(year, k=45)
#' fit4 <- sca(ple4, FLIndices(ple4.index), fmodel=fmod, qmodel=qmod, srmodel=srmod)
#' plot(ple4 + fit4)
#'
#' AIC(fit1, fit2, fit3, fit4)
#' BIC(fit1, fit2, fit3, fit4)
setGeneric("sca", function(stock, indices, ...) standardGeneric("sca"))

#' @rdname sca
setMethod("sca", signature("FLStock", "FLIndex"),
  function(stock, indices, ...) {
    sca(stock, FLIndices(IND=indices), ...)
  }
)

#' @rdname sca
setMethod("sca", signature("FLStock", "FLIndices"), function(stock, indices, fmodel, qmodel, srmodel = ~ factor(year), fit = "MP", mcmc=missing, useADMB = FALSE)
{
  stkdms <- dims(stock)
  if(stkdms$quant != "age"){
    stop("FLStock object must be age based") #Happy now, Rob?
  }
  if(missing(fmodel)){
	  # set up a robust model...  this needs some work...
	  ka <- stkdms$age
	  #ka <- if (ka < 3) ka else min(max(3, floor(.5 * ka)), 6)
	  ky <- floor(.5 * stkdms$year)
	  if (ka >= 3) {
		ka <- min(max(3, floor(.4 * ka)), 6)
		kb <- min(max(3, floor(2 * ka)), 10)
	    #fmodel <- formula(paste("~ te(age, year, k = c(", ka,",", ky,"), bs = 'tp')"))
	    fmodel <- formula(paste("~ te(age, year, k = c(", ka,",", ky,"), bs = 'tp') + s(age, k=", kb, ")"))
	  } else {
		fmodel <- formula(paste("~ age + s(year, k = ", ky,")"))
	  }
  }

  if(missing(qmodel)){
	  inddms <- lapply(indices, dims)
	  ka <- sapply(inddms, "[[", "age")
	  qmodel <- lapply(ka, function(i){
	  	# one age ~1 2 ages ~ age 3 or more s(age)
	  	if (i == 1) ~ 1 else if (i <= 3) ~ age else formula(paste("~ s(age, k = ", pmax(3, pmin(floor(.7 * ka), 7)), ")"))
	  })
  }

  if (stkdms$age > 10) {
    n1model <- ~ s(age, k = 10)
  } else {
    n1model <- ~ factor(age)
  }

  vmodel  <- lapply(seq(length(indices) + 1), function(i) ~ 1)
  vmodel[[1]] <- ~ s(age, k = 3)

  a4aSCA(fmodel = fmodel, qmodel = qmodel, srmodel = srmodel, n1model = n1model, vmodel = vmodel, stock = stock, indices = indices, fit = fit, useADMB = useADMB)

})

#setMethod("sca", signature("FLStock", "FLIndex"), function(stock, indices, fmodel  = "missing", qmodel  = "missing", srmodel = ~ factor(year), fit = "MP")
#{
#browser()
#	if(missing(fmodel)){
#		sca(stock=stock, indices=FLIndices(indices), qmodel=qmodel, srmodel=srmodel, fit=fit)
#	}
#	if(missing(qmodel)){
#		sca(stock=stock, indices=FLIndices(indices), fmodel=fmodel, srmodel=srmodel, fit=fit)
#	}
#	if(missing(qmodel) & missing(fmodel)){
#		sca(stock=stock, indices=FLIndices(indices), srmodel=srmodel, fit=fit)
#	} else {
#		sca(stock=stock, indices=FLIndices(indices), fmodel=fmodel, qmodel=qmodel, srmodel=srmodel, fit=fit)
#	}
#})

#' @title Statistical catch-at-age advanced method
#' @name a4aSCA
#' @docType methods
#' @rdname a4aSCA
#' @description Advanced user interface to the statistical catch-at-age method of the a4a stock assessment framework.
#' @details This method is the advanced method for stock assessment, it gives the user access to a set of arguments that the \code{sca} method doesn't. In particular, the default for the \code{fit} argument is 'assessment'. For detailed information about using the \code{sca} read the vignette 'The a4a Stock Assessment Modelling Framework' (\code{vignette('sca')}).
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
#' @param useADMB if FALSE (default) use TMB if TRUE use ADMB (note MCMC cannot be used with TMB)
#' @template dots
#' @return an \code{a4aFit} object if fit is "MP" or an \code{a4aFitSA} object if fit is "assessment"
#' @aliases a4aSCA a4aSCA-methods
#' @examples
#' data(ple4)
#' data(ple4.index)
#'
#' # fishing mortality by age and year (separable) AND catchability at age without year trend
#' fmodel <- ~factor(age) + factor(year)
#' qmodel <- list(~factor(age))
#' fit1 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#'
#' # fishing mortality as a smoother by age and year (but still separable) AND
#' # catchability at age without year trend
#' fmodel <- ~ s(age, k=4) + s(year, k=10)
#' qmodel <- list(~factor(age))
#' fit2 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#'
#' # fishing mortality as a smoother by age and year (but still separable) AND
#' # catchability as a smoother by age without year trend
#' fmodel <- ~ s(age, k=4) + s(year, k=10)
#' qmodel <- list(~s(age, k=4))
#' fit3 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#'
#' # fishing mortality as a smoother by age and year (but still separable) AND
#' # catchability as a smoother by age with year trend
#' fmodel <- ~ s(age, k=4) + s(year, k=10)
#' qmodel <- list(~s(age, k=4) + year)
#' fit4 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#'
#' # It's a statistical model
#' BIC(fit1, fit2, fit3, fit4)
#'
#' # fishing mortality as a smoother by age and year with interactions (i.e. non-separable) AND
#' # catchability as a smoother by age without year trend
#' fmodel <- ~ te(age, year, k=c(4, 10))
#' qmodel <- list(~s(age, k=4))
#' fit5 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#'
#' # fit3 + smoother in recruitment
#' fmodel <- ~ s(age, k=4) + s(year, k=20)
#' qmodel <- list(~s(age, k=4))
#' rmodel <- ~s(year, k=20)
#' fit6 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, srmodel=rmodel, ple4, FLIndices(ple4.index))
#'
#' # fit3 + bevholt
#' rmodel <- ~ bevholt(CV=0.05)
#' fit7 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, srmodel=rmodel, ple4, FLIndices(ple4.index))
setGeneric("a4aSCA", function(stock, indices, ...) standardGeneric("a4aSCA"))

#' @rdname a4aSCA
setMethod("a4aSCA", signature("FLStock", "FLIndex"),
  function(stock, indices, ...) {
    a4aSCA(stock, FLIndices(IND=indices), ...)
  }
)

#' @rdname a4aSCA
#setMethod("a4aSCA", signature("FLStock", "FLIndices"), function(stock, indices, fmodel  = ~ s(age, k = 3) + factor(year), qmodel  = lapply(seq(length(indices)), function(i) ~ 1), srmodel = ~ factor(year), n1model = ~ factor(age), vmodel  = missing, covar=missing, wkdir=missing, verbose = FALSE, fit = "assessment", center = TRUE, mcmc=missing) {

setMethod("a4aSCA", signature("FLStock", "FLIndices"),
  function(stock, indices, fmodel, qmodel, srmodel = ~ factor(year), n1model, vmodel,
           covar = missing, wkdir = missing, verbose = FALSE, fit = "assessment", center = TRUE, mcmc = missing, useADMB = FALSE)
{
  fit <- match.arg(fit, c("MP", "assessment", "MCMC"))

#  # set up default for vmodel
#  if (missing(vmodel)) {
#    vmodel  <- lapply(seq(length(indices) + 1), function(i) ~ 1)
#    vmodel[[1]] <- ~ s(age, k = 3)
#  }

  #-----------------------------------------------------------------
  # default submodels, same as sca
  stkdms <- dims(stock)
  if(stkdms$quant != "age"){
    stop("FLStock object must be age based") #Happy now, Rob?
  }

  if(missing(fmodel)){
	  # set up a robust model...  this needs some work...
	  ka <- stkdms$age
	  #ka <- if (ka < 3) ka else min(max(3, floor(.5 * ka)), 6)
	  ky <- floor(.5 * stkdms$year)
	  if (ka >= 3) {
		ka <- min(max(3, floor(.4 * ka)), 6)
		kb <- min(max(3, floor(2 * ka)), 10)
	    #fmodel <- formula(paste("~ te(age, year, k = c(", ka,",", ky,"), bs = 'tp')"))
	    fmodel <- formula(paste("~ te(age, year, k = c(", ka,",", ky,"), bs = 'tp') + s(age, k=", kb, ")"))
	  } else {
		fmodel <- formula(paste("~ age + s(year, k = ", ky,")"))
	  }
  }

  if(missing(qmodel)){
	  inddms <- lapply(indices, dims)
	  ka <- sapply(inddms, "[[", "age")
	  qmodel <- lapply(ka, function(i){
	  	# one age ~1 2 ages ~ age 3 or more s(age)
	  	if (i == 1) ~ 1 else if (i <= 3) ~ age else formula(paste("~ s(age, k = ", pmax(3, pmin(floor(.7 * ka), 7)), ")"))
	  })
  }

  if(missing(n1model)){
	  if (stkdms$age > 10) {
		n1model <- ~ s(age, k = 10)
	  } else {
		n1model <- ~ factor(age)
	  }
  }

  if(missing(vmodel)){
	  vmodel  <- lapply(seq(length(indices) + 1), function(i) ~ 1)
	  vmodel[[1]] <- ~ s(age, k = 3)
  }

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
	    outi <- a4aInternal(fmodel = fmodel, qmodel = qmodel, srmodel = srmodel, n1model = n1model, vmodel = vmodel, stock = istock, indices = iindices, covar = icovar, wkdir = wkdir, verbose = verbose, fit = ifit, center = center, mcmc=mcmc, useADMB = useADMB)
    } else if(!missing(covar) & missing(wkdir)){
      icovar <- lapply(covar, function(x) x[,, grid$unit[i], grid$area[i], , min(grid$iter[i], dims(x)$iter)])
	    outi <- a4aInternal(fmodel = fmodel, qmodel = qmodel, srmodel = srmodel, n1model = n1model, vmodel = vmodel, stock = istock, indices = iindices, covar = icovar, verbose = verbose, fit = ifit, center = center, mcmc=mcmc, useADMB = useADMB)
    } else if(missing(covar) & !missing(wkdir)){
	    outi <- a4aInternal(fmodel = fmodel, qmodel = qmodel, srmodel = srmodel, n1model = n1model, vmodel = vmodel, stock = istock, indices = iindices, wkdir=wkdir, verbose = verbose, fit = ifit, center = center, mcmc=mcmc, useADMB = useADMB)
	  } else {
	    outi <- a4aInternal(fmodel = fmodel, qmodel = qmodel, srmodel = srmodel, n1model = n1model, vmodel = vmodel, stock = istock, indices = iindices, verbose = verbose, fit = ifit, center = center, mcmc=mcmc, useADMB = useADMB)
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
          out@pars@qmodel[[j]]@centering[,i] <- outi@pars@qmodel[[j]]@centering
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

#' @title Stock assessment model advanced method
#' @name a4aInternal
#' @docType methods
#' @rdname a4aInternal
#' @description The advanced user interface to the a4a fitting routine.
#'
#' @param fmodel a formula object depicting the model for log fishing mortality at age
#' @param qmodel a list of formula objects depicting the models for log survey catchability at age
#' @param srmodel a formula object depicting the model for log recruitment
#' @param n1model a formula object depicting the model for the first year of catch data
#' @param vmodel a list of formula objects depicting the models for log survey and log fishing mortality variance
#' @param stock an FLStock object containing catch and stock information
#' @param indices an FLIndices object containing survey indices
#' @param covar a list with covariates
#' @param wkdir used to set a working directory for the admb optimiser.  If wkdir is set all admb files are saved to this folder otherwise they are deleted.
#' @param verbose if true admb fitting information is printed to the screen
#' @param fit character with type of fit: 'MP' or 'assessment'; the former doesn't require the hessian to be computed, while the latter does.
#' @param center \code{logical} specifying whether data is centered before estimating or not
#' @param mcmc \code{SCAMCMC} specifying parameters for the ADMB MCMC run, check ADMB manual for detailed description
#' @param useADMB if FALSE (default) use TMB if TRUE use ADMB (note MCMC cannot be used with TMB)
#' @return an \code{a4aFit} object if fit is "MP" or an \code{a4aFitSA} if fit is "assessment"
#' @aliases a4aInternal
a4aInternal <- function(stock, indices, fmodel  = ~ s(age, k = 3) + factor(year),
                qmodel  = lapply(seq(length(indices)), function(i) ~ 1),
                srmodel = ~ factor(year),
                n1model = ~ factor(age),
                vmodel  = lapply(seq(length(indices) + 1), function(i) ~ 1),
                covar=missing, wkdir=missing, verbose = FALSE, fit = "assessment",
                center = TRUE, mcmc=missing, useADMB = FALSE)
{

  # first check permissions of executable
  #	exeok <- check.executable()
  #	if (!exeok) stop("a4a executable has wrong permissions.")

  # if fit MCMC mcmc object must exist
  if(fit=="MCMC" & missing(mcmc)) stop("To run MCMC you need to set the mcmc argument, using the method SCAMCMC.")

  # start timer
	my.time.used <- numeric(4)
	my.time.used[1] <- Sys.time()

	# some quick friendly things
	if (!inherits(indices, "FLIndices")) indices <- FLIndex(indices)
	if (!inherits(catch.n(stock),"FLQuantDistr")) {
		varslt <- catch.n(stock)
		varslt[] <- NA
		catch.n(stock) <- FLQuantDistr(catch.n(stock), varslt)
	}

	# check survey names
	if (length(names(indices)) == 0) {
		# check that survey names don't use ':'
		if(sum(grepl(':', names(indices)))>0) stop('Indices names can\'t use the character \':\'')
		snames <- make.unique(rep("survey", length(indices)))
	} else {
		snames <- make.unique(names(indices))
	}

	for (i in seq_along(indices)) name(indices[[i]]) <- snames[i]
	names(indices) <- snames

	# what kind of run is this
	fit <- match.arg(fit, c("MP", "assessment", "MCMC", "setup"))

	#========================================================================
	# Extract observations and auxilliary info from stock and indices objects
	#========================================================================
	# first some checks
	if (any(is.infinite(log(catch.n(stock))))) stop("only non-zero catches allowed.")
	if (any(is.infinite(log(var(catch.n(stock)))))) stop("only non-zero catch variances allowed.")
	if (any(is.infinite(log( unlist(lapply(indices, function(x) c(index(x)))) ))))  stop("only non-zero survey indices allowed.")
	if (any(is.infinite(log( unlist(lapply(indices, function(x) c(index.var(x)))) ))))  stop("only non-zero survey index variances allowed.")

	# convert catches and indices to a list of named arrays
	list.obs <- c(list(catch = quant2mat(catch.n(stock)@.Data)), lapply(indices, function(x) quant2mat(index(x)@.Data)))

	# convert the variances of catches and indices to a list of named arrays
	list.var <- c(list(catch = quant2mat(catch.n(stock)@var)), lapply(indices, function(x) quant2mat(index.var(x))))

	# calculate appropriate centering for observations on log scale
	# a bit spaguetti ... if center is a numeric vector only those elements will be centered
	center.log <- sapply(list.obs, function(x) mean(log(x), na.rm = TRUE))
	if(is.numeric(center)) center.log[-center][] <- 0 else if(!isTRUE(center)) center.log[] <- 0

	# convert to dataframe. NOTE: list2df also logs the observations and centers
	df.data <- do.call(rbind, lapply(1:length(list.obs), list2df, list.obs=list.obs, list.var=list.var, center.log=center.log))

	if (any(df.data[,5] != 1)) message("Note: Provided variances will be used to weight observations.\n\tWeighting assumes variances are on the log scale or equivalently log(CV^2 + 1).")

	# extract auxilliary stock info
	fbar <-  unname(range(stock)[c("minfbar","maxfbar")])
	plusgroup <- as.integer( !is.na(range(stock)["plusgroup"]), range(stock)["plusgroup"] >= range(stock)["max"] )

	# extract auxilliary survey info - always assume oldest age is true age TODO TODO TODO !!
	surveytime <- unname(sapply(indices, function(x) mean(c(dims(x)$startf, dims(x)$endf))))
	names(surveytime) <- names(indices)
	if (any(is.na(surveytime))) stop("You need to define startf and endf for each index!!")

	#========================================================================
	# Make a full data.frame and add in covariates and observations
	#========================================================================

	# build a full data frame first (we will use this for the variance model so it is not a waste)
	full.df <- do.call(rbind, lapply(1:length(list.obs), function(i) cbind(fleet = i, make.df(i, stock=stock, indices=indices))))

	#-------------------------------------------------------------------------
	# NOTE: check covar object. Need to be consistent on how this information
	# is passed to the method. Currently sometimes through the covar object
	# others directly with vectors.
	#-------------------------------------------------------------------------
	if (!missing(covar)) {
		# add in covariates to data.frame - it is easiest to provide covariates in one list
		tmp <- lapply(seq_along(covar), function(i) {
			x <- as.data.frame(covar[[i]])[c(1,2,7)]
			if (length(unique(x$age)) == 1) x <- x[names(x) != "age"]
			if (length(unique(x$year)) == 1) x <- x[names(x) != "year"]
			names(x) <- gsub("data", names(covar)[i], names(x))
			x
		})
		covar.df <- tmp[[1]]
		for (i in seq_along(tmp[-1])) covar.df <- merge(covar.df, i, all = TRUE, sort = FALSE)

		full.df <- merge(full.df, covar.df, all.x = TRUE, all.y = FALSE)
	}
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	# add in data
	full.df <- merge(full.df, df.data, all.x = TRUE, all.y = FALSE)
	# put biomass surveys in min age position
	full.df$age <- with(full.df, replace(age, is.na(age), min(age, na.rm = TRUE)))
	full.df$fleet <- factor(names(list.obs)[full.df$fleet], levels = names(list.obs))

	# set weights for missing values to zero
	full.df$weights[is.na(full.df$obs)] <- 0
	# inform that missing values will be treated as missing at random
	if(verbose){
		if (any(is.na(full.df$obs))){
			message("Note: The following observations are treated as being missing at random:\n\t",
			paste(capture.output(print(subset(full.df,is.na(obs))[c("fleet","year","age")], row.names = FALSE)),
			collapse = "\n\t"), "\n      Predictions will be made for missing observations." )
		}
	}

	# fill out data frame for all eventualities ... except ages ....
	# or a quick one that would not mess array dims and things later is to fix ages in data frame but keep rows...
	# it would actually be easiest to present all auxilliary data and covariates as a data.frame and the observations as a matrix...
	temp.full.df <- expand.grid(lapply(full.df[3:1],function(x) sort(unique(x))),  KEEP.OUT.ATTRS = FALSE)[3:1]
	for (i in seq_along(indices)) {
		# if biomass survey skip this step
		if (is(indices[[i]], 'FLIndexBiomass')) next
			.ages <- temp.full.df$age [temp.full.df$fleet == levels(full.df$fleet)[i+1]]
			.range <- range(subset(full.df, fleet == levels(full.df$fleet)[i+1])$age)
			.ages[.ages < .range[1]] <- .range[1]
			.ages[.ages > .range[2]] <- .range[2]
			temp.full.df$age [temp.full.df$fleet == levels(full.df$fleet)[i+1]] <- .ages
	}
	full.df <- merge(temp.full.df, full.df, all.x = TRUE)


	# add in auxilliary data - maturity, natural mortality etc.
	aux.df <- cbind(as.data.frame(stock.n(stock))[c("year","age")],
		m = log(as.data.frame(m(stock))$data), mat = as.data.frame(mat(stock))$data,
		stock.wt = as.data.frame(stock.wt(stock))$data, catch.wt = as.data.frame(catch.wt(stock))$data,
		m.spwn = as.data.frame(m.spwn(stock))$data, harvest.spwn = as.data.frame(harvest.spwn(stock))$data)

	full.df <- merge(full.df, aux.df, all.x = TRUE, all.y = FALSE)
	full.df$mat.wt <- full.df$stock.wt * full.df$mat

	# set unspecified covariates and aux data to 0
	# TODO we should warn about this... or check before the merge that there are covars and m and mat for requested predictions
	full.df[is.na(full.df)] <- 0

	# order df
	full.df <- full.df[with(full.df, order(fleet, year, age)), ]
	rownames(full.df) <- NULL
	full.df <- full.df[c(3,1,2,4:ncol(full.df))]

	# TODO stop if some obs have no parameter - it should be quite rare ...

	#========================================================================
	# Process model formulas
	#========================================================================

	# F model matrix
	Xf <- getX(fmodel, subset(full.df, fleet == "catch"))
	# F model offsets
	# ...

	# Q model matrix
	fleet.names <- c("catch", names(indices))
	Xqlist <- lapply(seq_along(indices), function(i) getX(qmodel[[i]], subset(full.df, fleet == fleet.names[i+1])))
	Xq <- as.matrix(do.call(bdiag, Xqlist))
	# Q model offsets
	# ...

	# var model matrix
	Xvlist <- lapply(1:length(fleet.names), function(i) getX(vmodel[[i]], subset(full.df, fleet == fleet.names[i])))
	Xv <- as.matrix(do.call(bdiag, Xvlist))
	# var model offsets
	# ...

	# initial age structure model matrix
	Xny1 <- getX(n1model, subset(full.df, year == min(year) & age > min(age) & fleet == "catch"))

	#-------------------------------------------------------------------------
	# NOTE: how are covars being passed ?
	#-------------------------------------------------------------------------
	# now separate model elements, find SR models, stop if more than one specified
	a4as <- isPresenta4aSRmodel(srmodel)
	if (sum(a4as) > 1) stop("you can only specify one type of stock recruit relationship.")
	srrmod <- geta4aSRmodel(srmodel)
	if (sum(a4as) == 0 && max(full.df$year) > max(df.data$year)) stop("you need to specify a stock recruitment relationship to forecast with out survey information.")

	# extract a and b model formulas and add on any extra bits to amodel.
	# NB SRR models should be parametrised so that amodel is the level of recruitment!!!
	srr <- eval(parse(text = srrmod))
	if (sum(a4as) > 0 && any(!a4as)) {
		# ignore .. the following line adds these onto the end of amod
		message("Note: Trailing formula elements in the srmodel have been removed")
		#srr$amodel <- eval(parse(text = paste("~", as.character(srr$amodel)[length(srr$amodel)], "+", paste(facs[!a4as], collapse = " + ")) ))
	}
	Xsra <- getX(srr$a, subset(full.df, fleet == "catch" & age == dims(stock)$min))
	Xsrb <- getX(srr$b, subset(full.df, fleet == "catch" & age == dims(stock)$min))

	# can we do a quick check for identifiability of the srr model? ...
	if (ncol(Xsra) + ncol(Xsrb) > dims(stock)$year) stop("Stock recruitment model is over parameterised, please reduce the number parameters")

	# internal r model matrix - this is setting the X model for the recruitments,
	# it is not the same as the sr model which is the model for the relationship
	# between the recruitments and SSB
	if (sum(a4as) == 0) rmodel <- srmodel else rmodel <- ~ factor(year)
	Xr <- getX(rmodel, subset(full.df, age == min(age) & fleet == "catch"))

  #========================================================================
  # Fit the model and return a list of objects detailing the fit
  #========================================================================

  # change NA to -1 for admb
  df.data$age <- with(df.data, replace(age, is.na(age), -1))

  # build survey's max and min age vectors
  # If age is not set (NA) it will take the min and max from catch.n
  srvRange <- do.call('rbind',lapply(indices, range))
  srvMinAge <- srvRange[,'min']
  srvMinAge[is.na(srvMinAge)] <- range(full.df$age)[1]
  names(srvMinAge) <- names(indices)
  srvMaxAge <- srvRange[,'max']
  srvMaxAge[is.na(srvMaxAge)] <- range(full.df$age)[2]
  names(srvMaxAge) <- names(indices)



  if (useADMB) { # fit using the ADMB code

    #========================================================================
    # Create the directory where to store model config, data and results files
    #========================================================================
    #-------------------------------------------------------------------------
    # NOTE: move to internal funs
    #-------------------------------------------------------------------------
    if (missing(wkdir)) keep <- FALSE else keep <- TRUE # keep results if wkdir is set by user
    if (keep) {
      # create the directory locally - whereever specified by the user
      wkdir.start <- wkdir
      # if this directory already exists, try the numbered versions
      kk <- 1
      ans <- file.exists(wkdir)
      while(ans) {
        wkdir <- paste(wkdir.start,"-", kk, sep = "")
        kk <- kk + 1
        ans <- file.exists(wkdir)
      }
      # if several a4aFit()'s are run in parallel, we might have a
      # conflict. if so, create a random name
      if (file.exists(wkdir)) {
        wkdir <- paste(wkdir, "-", substring(as.character(runif(1)), 3), sep = "")
      }
      cat("Model and results are stored in working directory [", wkdir,"]\n", sep = "")
    } else {
      # no wkdir specified by user so create a temporary directory
      wkdir <- tempfile('file', tempdir())
    }

    # Create a directory where to store data and results
    # TODO check if wkdir exists
    dir.create(wkdir, showWarnings = FALSE)
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    fitList <- fitADMB(fit, wkdir, df.data, stock, indices, full.df,
                       fbar, plusgroup, surveytime, fleet.names,
                       Xf, Xq, Xv, Xny1, srr, Xsra, Xsrb, Xr, Xvlist, Xqlist,
                       my.time.used, mcmc, verbose)

    if (fit == "setup") return(fitList)
    out <- fitList$out
    my.time.used <- fitList$my.time.used
    wkdir <- fitList$wkdir
    convergence <- fitList$convergence
    pnames <- fitList$pnames
    ages <- fitList$ages
    years <- fitList$years

    # remove temporary directory - keep only true when dir is not temp dir
    if (!keep) unlink(wkdir, recursive=TRUE, force=TRUE)

  } else { # fit using the TMB code
    fitList <- fitTMB(fit, wkdir, df.data, stock, indices, full.df,
                      fbar, plusgroup, surveytime, fleet.names,
                      Xf, Xq, Xv, Xny1, srr, Xsra, Xsrb, Xr, Xvlist, Xqlist,
                      my.time.used, mcmc, verbose)

    if (fit == "setup") return(fitList)

    out <- fitList$out
    my.time.used <- fitList$my.time.used
    wkdir <- fitList$wkdir
    convergence <- fitList$convergence
    pnames <- fitList$pnames
    ages <- fitList$ages
    years <- fitList$years
  }

	#========================================================================
	# convert into appropriate class
	#========================================================================

	#------------------------------------------------------------------------
	# a4aFit
	#------------------------------------------------------------------------

	a4aout <- a4aFit()
	ind.names <- names(indices)

	stk.n <- t(out$N)
	names(dimnames(stk.n)) <- c("age","year")

	hvst <- t(out$F)
	names(dimnames(hvst)) <- c("age","year")

	logq <- lapply(1:length(indices), function(i) {
		x <- drop(out$Q[,,i])
		names(dimnames(x)) <- c("age","year")
		if (is(indices[[i]], 'FLIndexBiomass')){
			x <- with(dimnames(index(indices[[i]])), x[1, year,drop = FALSE])
			dimnames(x)[[1]] <- "all"
			x <- FLQuant(x, dimnames=dimnames(index(indices[[i]])))
		} else {
			x <- with(dimnames(index(indices[[i]])), x[age, year])
			x <- FLQuant(x, dimnames=dimnames(index(indices[[i]])))
		}
		x
	})
	names(logq) <- ind.names

	a4aout@name    <- stock@name
	a4aout@desc    <- stock@desc
	a4aout@range   <- stock@range
	a4aout@call    <- match.call()
	a4aout@stock.n <- FLQuant(stk.n) * exp(center.log[1])
	a4aout@harvest <- FLQuant(hvst, units = "f")
	Z <- a4aout@harvest + m(stock)
	a4aout@catch.n <- a4aout@harvest / Z * (1 - exp(-Z)) * a4aout@stock.n

	index <- lapply(1:length(indices), function(i) {
	 	dmns <- dimnames(logq[[i]])
	 	if (is(indices[[i]], 'FLIndexBiomass')) {
		dmns[[1]] <- ac(srvMinAge[i]:srvMaxAge[i])
	 		qq <- exp(logq[[i]] - center.log[1] + center.log[i+1])
	 		nn <- stock.n(a4aout)[dmns[[1]], dmns[[2]]]
	 		zz <- exp(-Z[dmns[[1]], dmns[[2]]]*surveytime[i])
	 		# there was a bug here, nn was not being reduced by Z
	 		bb <- apply(nn*zz*stock.wt(stock)[dmns[[1]], dmns[[2]]], 2:6, sum, na.rm=T)
	 		ii <- qq*bb
	 	} else {
	 		nn <- stock.n(a4aout)[dmns[[1]], dmns[[2]]]
	 		qq <- exp(logq[[i]] - center.log[1] + center.log[i+1])
	 		zz <- exp(-Z[dmns[[1]], dmns[[2]]]*surveytime[i])
	 		ii <- qq*nn*zz
	 	}
	 	ii
	})
	names(index) <- ind.names
	a4aout@index <- FLQuants(index)

	# GCV (Wood, 2006, pag. 132)
	flev <- diag(Xf %*% solve(t(Xf) %*% Xf) %*% t(Xf))
	cgcv <- prod(dim(a4aout@catch.n)) * sum(c(log(catch.n(stock)/a4aout@catch.n))^2)/sum(1-flev)^2
	tmpSumm <- with(out, c(nopar, nlogl, maxgrad, nrow(df.data), cgcv, convergence, NA))
	a4aout@fitSumm <- array(tmpSumm, dimnames = list(c("nopar","nlogl","maxgrad","nobs","gcv", "convergence", "accrate")))

	#------------------------------------------------------------------------
	# a4aFitSA
	#------------------------------------------------------------------------

	if (fit %in% c("assessment", "MCMC")) {

		# coerce into a4aFitSA
		a4aout <- a4aFitSA(a4aout)

		# fill up stkmodel
		a4aout@pars@stkmodel@name      <- a4aout@name
		a4aout@pars@stkmodel@desc      <- a4aout@desc
		a4aout@pars@stkmodel@range     <- a4aout@range
		a4aout@pars@stkmodel@centering <- FLPar(centering = center.log[1])
		a4aout@pars@stkmodel@fMod      <- fmodel
		a4aout@pars@stkmodel@n1Mod     <- n1model
		a4aout@pars@stkmodel@srMod     <- srmodel
		a4aout@pars@stkmodel@m         <- m(stock)
		a4aout@pars@stkmodel@wt        <- stock.wt(stock)
    a4aout@pars@stkmodel@mat       <- mat(stock)
		a4aout@pars@stkmodel@link      <- log
		a4aout@pars@stkmodel@linkinv   <- exp
		pars <- out$par.est
		active <- sapply(out$par.std, length) > 0
		# add in iter dimension
		dim(out$cov) <- c(dim(out$cov), 1)
		dimnames(out$cov) <- list(unlist(pnames), unlist(pnames), 1)
		stkactive <- active
		stkactive[2:3] <- FALSE
		a4aout@pars@stkmodel@coefficients <- FLPar(structure(unlist(pars[stkactive]), names = unlist(pnames[stkactive])))
		units(a4aout@pars@stkmodel@coefficients) <- "NA"
		a4aout@pars@stkmodel@distr <- "norm"
		whichcol <-  grep(paste("(^",c("f","n1","r","sra","srb"),"Mod:)",collapse="|",sep=""), unlist(pnames))
		# we can use the inverse of a subset of the precision matrix
		# if we want the vcov matrix conditional on the
		# other (qmodel etc.) parameter estimates.
		##a4aout@pars@stkmodel@vcov <- solve(out$prec[whichcol, whichcol])
		# or just the full vcov matrix, unconditional on the other things...
		a4aout@pars@stkmodel@vcov <- out$cov[whichcol, whichcol, 1, drop=FALSE]

		# fill up qmodel
		qmodels <- lapply(seq_along(indices), function(i){
			which <- sapply(strsplit(pnames[[2]], split=":"), "[[", 2) %in% fleet.names[i+1]
			submodel(formula = qmodel[[i]],
				       coefficients = FLPar(structure(pars[[2]][which], names = pnames[[2]][which])),
				       vcov = out$cov[pnames[[2]][which], pnames[[2]][which], 1, drop = FALSE],
				       distr = "norm",
				       centering = FLPar(centering = center.log[i+1] - center.log[1]), # must also subtract catch scaling
				       name = fleet.names[i+1],
				       desc = indices[[i]]@desc,
				       range = indices[[i]]@range,
				       link = log,
				       linkinv = exp
			)
		})

		names(qmodels) <- fleet.names[-1]
		a4aout@pars@qmodel <- submodels(qmodels)

		# fill up vmodel
		vmodels <- lapply(seq_along(fleet.names), function(i){
			which <- sapply(strsplit(pnames[[3]], split=":"), "[[", 2) %in% fleet.names[i]
			submodel(formula = vmodel[[i]],
				       coefficients = FLPar(structure(pars[[3]][which], names = pnames[[3]][which])),
				       vcov = out$cov[pnames[[3]][which], pnames[[3]][which], 1, drop = FALSE],
				       distr = "norm",
				       centering = FLPar(centering = 0),
				       name = fleet.names[i],
				       desc = "",
				       range = if (i==1) stock@range else indices[[i-1]]@range,
				       link = log,
				       linkinv = exp
			)
		})

		names(vmodels) <- fleet.names
		a4aout@pars@vmodel <- submodels(vmodels)
	}

	#------------------------------------------------------------------------
	# a4aFitMCMC
	# To build this class we're making use of the SA class, whih is not the
	# most efficient code but it's easier to implement.
	#------------------------------------------------------------------------

	if (fit == "MCMC") {

    if (!useADMB) stop("not implemented in TMB version")

		# coerce into a4aFitMCMC
		a4aout <- a4aFitMCMC(a4aout, mcmc=mcmc)

		# fill parameters
		a4aout@pars@stkmodel@params <- propagate(a4aout@pars@stkmodel@params, out$mcmc$nit)
		a4aout@pars@stkmodel@params[] <- t(out$mcmc$mcmcout[,dimnames(a4aout@pars@stkmodel@params)[[1]]])
		dimnames(a4aout@pars@stkmodel@params)[[2]] <- 1:out$mcmc$nit
		a4aout@pars@stkmodel@vcov[] <- NA

		# fill derived quantities
		a4aout@stock.n <- propagate(a4aout@stock.n, out$mcmc$nit)
		a4aout@stock.n[] <- c(t(out$mcmc$N * exp(center.log[1])))

		a4aout@harvest <- propagate(a4aout@harvest, out$mcmc$nit)
		a4aout@harvest[] <- c(t(out$mcmc$F))

		Z <- a4aout@harvest + m(stock)
		a4aout@catch.n <- a4aout@harvest / Z * (1 - exp(-Z)) * a4aout@stock.n

		# fill catchability and variance model
		for(i in fleet.names){
			a4aout@pars@vmodel[[i]]@params <- propagate(a4aout@pars@vmodel[[i]]@params, out$mcmc$nit)
			a4aout@pars@vmodel[[i]]@params[] <- t(out$mcmc$mcmcout[,dimnames(a4aout@pars@vmodel[[i]]@params)[[1]]])
			a4aout@pars@vmodel[[i]]@vcov[] <- NA

			if(i!="catch"){
				a4aout@pars@qmodel[[i]]@params <- propagate(a4aout@pars@qmodel[[i]]@params, out$mcmc$nit)
				a4aout@pars@qmodel[[i]]@params[] <- t(out$mcmc$mcmcout[,dimnames(a4aout@pars@qmodel[[i]]@params)[[1]]])
				a4aout@pars@qmodel[[i]]@vcov[] <- NA
				dmns <- dimnames(a4aout@index[[i]])
				a4aout@index[[i]] <- propagate(a4aout@index[[i]], out$mcmc$nit)
			 	if (is(indices[[i]], 'FLIndexBiomass')) {
					a4aout@index[[i]][] <- t(out$mcmc$QQ[out$mcmc$idq$s==i & out$mcmc$idq$y %in% as.numeric(dmns[[2]]),1, drop=FALSE])
					dmns[[1]] <- ac(srvMinAge[i]:srvMaxAge[i])
					a4aout@index[[i]] <- exp(a4aout@index[[i]] - center.log["catch"] + center.log[i])
					nn <- stock.n(a4aout)[dmns[[1]], dmns[[2]]]
					zz <- exp(-Z[dmns[[1]], dmns[[2]]]*surveytime[i])
			 		bb <- apply(nn*zz*stock.wt(stock)[dmns[[1]], dmns[[2]]], 2:6, sum, na.rm=TRUE)
					a4aout@index[[i]] <- a4aout@index[[i]]*bb
				} else {
					a4aout@index[[i]][] <- t(out$mcmc$QQ[out$mcmc$idq$s==i & out$mcmc$idq$y %in% as.numeric(dmns[[2]]),ages %in% dmns[[1]]])
					a4aout@index[[i]] <- exp(a4aout@index[[i]] - center.log["catch"] + center.log[i])
					nn <- stock.n(a4aout)[dmns[[1]], dmns[[2]]]
					zz <- exp(-Z[dmns[[1]], dmns[[2]]]*surveytime[i])
					a4aout@index[[i]] <- a4aout@index[[i]]*nn*zz
				}
			}
		}


		tmpSumm <- with(out, c(nopar, NA, NA, nrow(df.data), NA, NA, as.numeric(mcmc$logfile[length(mcmc$logfile)-2])))
		a4aout@fitSumm <- array(tmpSumm, dimnames = list(c("nopar","nlogl","maxgrad","nobs","gcv", "convergence", "accrate")))
	}

	#========================================================================
	# return
	#========================================================================

	# end time
	my.time.used[4] <- Sys.time()
	a4aout@clock <- c("Pre-processing" = diff(my.time.used)[1], "Running a4a" = diff(my.time.used)[2], "Post-processing" = diff(my.time.used)[3], Total = my.time.used[4] - my.time.used[1])

	units(a4aout@harvest) <- "f"

	return(a4aout)
}

#' @title Breakpoints
#' @name breakpts
#' @rdname breakpts
#' @description Method to set breakpoints in submodels
#' @param var a \code{numeric} object that defines the variable to be "broken"
#' @param breaks a \code{numeric} object that defines the breakpoints
#' @template dots
#' @return a \code{factor} with levels according to the defined breaks
#' @aliases breakpts breakpts-methods
setGeneric("breakpts", function(var, ...) standardGeneric("breakpts"))
#' @rdname breakpts
setMethod("breakpts", "numeric", function(var, breaks, ...) {
  if (min(var, na.rm = TRUE) < min(breaks)) breaks <- c(min(var, na.rm = TRUE) - 1, breaks)
  if (max(var, na.rm = TRUE) > max(breaks)) breaks <- c(breaks, max(var, na.rm = TRUE))
  label <- paste0("(",breaks[-length(breaks)], ",", breaks[-1], "]")
  cut(var, breaks = breaks, label = label)
})


fitTMB <- function(fit, wkdir, df.data, stock, indices, full.df,
                    fbar, plusgroup, surveytime, fleet.names,
                    Xf, Xq, Xv, Xny1, srr, Xsra, Xsrb, Xr, Xvlist, Xqlist,
                    my.time.used, mcmc, verbose)
{
  # not used:
  # * wkdir

  #========================================================================
  # Prepare data
  #========================================================================
  # change NA to -1 for admb
  #df.data$age <- with(df.data, replace(age, is.na(age), -1))

  # set up variable names
  pnames <- list(paste0("fMod:",colnames(Xf)),
    paste0("qMod:", unlist(sapply(1:length(indices), function(i) paste0(fleet.names[i+1],":",colnames(Xqlist[[i]]))))),
    paste0("vMod:", unlist(sapply(1:length(fleet.names), function(i) paste0(fleet.names[i],":",colnames(Xvlist[[i]]))))),
    paste0("n1Mod:",colnames(Xny1)),
    paste0("rMod:",colnames(Xr)),
    paste0("sraMod:",colnames(Xsra)),
    paste0("srbMod:",colnames(Xsrb)))

  if (srr$srrCV < 0) pnames <- pnames[-c(6,7)]
  if (srr$ID==4) pnames <- pnames[-7]

  # build survey's max and min age vectors
  # If age is not set (NA) it will take the min and max from catch.n
  srvRange <- do.call('rbind',lapply(indices, range))
  srvMinAge <- srvRange[,'min']
  srvMinAge[is.na(srvMinAge)] <- range(full.df$age)[1]
  names(srvMinAge) <- names(indices)
  srvMaxAge <- srvRange[,'max']
  srvMaxAge[is.na(srvMaxAge)] <- range(full.df$age)[2]
  names(srvMaxAge) <- names(indices)

  df.aux <- unique(full.df[c("year","age","m","m.spwn","harvest.spwn","mat.wt","stock.wt")])

  Ldat <- list(
      ageRange = range(full.df$age),
      yearRange = range(full.df$year),
      surveyMinAge = srvMinAge,
      surveyMaxAge = srvMaxAge,
      surveyTimes = surveytime,
      fbarRange = fbar,
      obs = as.matrix(df.data),
      aux = as.matrix(df.aux),
      designF = Xf,
      designQ = Xq,
      designV = Xv,
      designNy1 = Xny1,
      designR = Xr,
      designRa = Xsra,
      designRb = Xsrb,
      srCV = srr$srrCV, # if (srr$srrCV < 0) 100 else srr$srrCV,
      spr0 = ifelse(!is.null(srr$SPR0), srr$SPR0, 0),
      Rmodel = srr$ID,
      isPlusGrp = plusgroup
    )

  Ldat$locFleetVec=Ldat$obs[,1]
  Ldat$locYearVec=Ldat$obs[,2]
  Ldat$locAgeVec=Ldat$obs[,3]

  Lpin <- list(
    fpar = rep(0, ncol(Ldat$designF)),
    qpar = rep(0, ncol(Ldat$designQ)),
    vpar = rep(0, ncol(Ldat$designV)),
    ny1par = rep(0, ncol(Ldat$designNy1)),
    rpar = rep(0, ncol(Ldat$designR)),
    rapar = rep(0, ncol(Ldat$designRa)),
    rbpar = rep(0, ncol(Ldat$designRb))
  )


  res <- list(Ldat=Ldat, Lpin=Lpin)

  #========================================================================
  # run model
  #========================================================================
  # run a4a split
  my.time.used[2] <- Sys.time()

  if (srr$srrCV < 0) {
    # we have no SR model
    map <- list()
    map$rapar <- factor(NA)
    map$rbpar <- factor(NA)
    obj <- MakeADFun(res$Ldat, res$Lpin, DLL="FLa4a", silent = !verbose, map = map)
  } else
  if (srr$ID == 4) {
    # its the geomean model
    map <- list()
    map$rbpar <- factor(NA)
    obj <- MakeADFun(res$Ldat, res$Lpin, DLL="FLa4a", silent = !verbose, map = map)
  } else {
    obj <- MakeADFun(res$Ldat, res$Lpin, DLL="FLa4a", silent = !verbose)
  }

  if (fit == "setup") return(obj)

  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 10000, eval.max = 10000))
  if (fit != "MP") {
    # The accuracy is good without these steps
    # but more consistent with these extra iterations
    newtonsteps <- 3 # fix to 3 for now...
    # this is required, otherwise we do not get the same
    # accuracy as the ADMB optimiser
    for (i in seq_len(newtonsteps)) { # Take a few extra newton steps borrowed from Anders/Casper!
      g <- as.numeric( obj$gr(opt$par) )
      h <- optimHess(opt$par, obj$fn, obj$gr)
      opt$par <- opt$par - solve(h, g)
      opt$objective <- obj$fn(opt$par)
    }
  }

  #rep <- obj$report()
  sdrep <- sdreport(obj,opt$par)
  sdrep$cov <- NULL # this is the cvov of the model predictions


  #========================================================================
  # Post process
  #========================================================================
  # post processing split
  my.time.used[3] <- Sys.time()

  if (fit %in% c("MP","assessment", "MCMC")) {

    # read admb output from file
    out <- list()

    # get .par file equivalents
    out$nopar <- length(opt$par)
    out$nlogl <- opt$objective
    out$maxgrad <- max(abs(sdrep$gradient.fixed))

    if (srr$srrCV < 0) {
      # we have no SR model
      out$par.est <- as.list(sdrep,"Est")[-c(6,7)]
    } else
    if (srr$ID == 4) {
      # geomean
      out$par.est <- as.list(sdrep,"Est")[-7]
    } else {
      out$par.est <- as.list(sdrep,"Est")
    }

    ages <- sort(unique(full.df$age))
    years <- sort(unique(full.df$year))

    # set up basic return values
    convergence <- opt$convergence
    out$N <- matrix(sdrep$value[names(sdrep$value) == "expn"], length(years), length(ages))
    out$F <- matrix(sdrep$value[names(sdrep$value) == "expf"], length(years), length(ages))
    out$Q <- matrix(sdrep$value[names(sdrep$value) == "q"], length(years)*length(indices), length(ages))
    dim(out$Q) <- c(length(years), length(indices), length(ages))
    out$Q <- aperm(out$Q, c(3,1,2))
    colnames(out$N) <- colnames(out$F) <- ages
    rownames(out$N) <- rownames(out$F) <- years
    dimnames(out$Q) <- list(ages, years, names(indices))

    if (fit != "MP") {
      # return hessian and standard errors
      if (!sdrep$pdHess) {
        # if model was not identifiable make sure return values are NAs
        warning("Hessian was not positive definite.")
        out$N[] <- NA
        out$F[] <- NA
        out$Q[] <- NA
        out$cov <- sdrep$cov.fixed
        out$cov[] <- NA
      } else {
        out$logDetHess <- NA
        if (srr$srrCV < 0) {
          # we have no SR model
          out$par.std <- as.list(sdrep,"Std")[-c(6,7)]
        } else
        if (srr$ID == 4) {
          # geomean
          out$par.std <- as.list(sdrep,"Std")[-7]
        } else {
          out$par.std <- as.list(sdrep,"Std")
        }
        out$cov <- sdrep$cov.fixed
      }
    }
  }

  list(out = out, my.time.used = my.time.used, wkdir = NULL,
       convergence = convergence, pnames = pnames,
       ages = ages, years = years)
}


fitADMB <- function(fit, wkdir, df.data, stock, indices, full.df,
                    fbar, plusgroup, surveytime, fleet.names,
                    Xf, Xq, Xv, Xny1, srr, Xsra, Xsrb, Xr, Xvlist, Xqlist,
                    my.time.used, mcmc, verbose)
{

  #========================================================================
  # Write model matrices and model info to files in wkdir
  #========================================================================

  #-------------------------------------------------------------------------
  # write to data file
  #-------------------------------------------------------------------------
  filename <- paste0(wkdir,'/a4a.dat')

  # change NA to -1 for admb
  df.data$age <- with(df.data, replace(age, is.na(age), -1))

  # build survey's max and min age vectors
  # If age is not set (NA) it will take the min and max from catch.n
  srvRange <- do.call('rbind',lapply(indices, range))
  srvMinAge <- srvRange[,'min']
  srvMinAge[is.na(srvMinAge)] <- range(full.df$age)[1]
  names(srvMinAge) <- names(indices)
  srvMaxAge <- srvRange[,'max']
  srvMaxAge[is.na(srvMaxAge)] <- range(full.df$age)[2]
  names(srvMaxAge) <- names(indices)

  cat("# Data for the a4a model",
    "\n# Full age range\n", range(full.df$age),
    "\n# Full year range\n", range(full.df$year),
    "\n# Number of surveys\n", length(unique(full.df$fleet)) - 1,
    "\n# Survey min ages\n", paste(srvMinAge, collapse = " "),
    "\n# Survey max ages\n", paste(srvMaxAge, collapse = " "),
    "\n# Survey time as a fraction into the year (one for each survey)\n", paste(surveytime, collapse = " "),
    "\n# fbar range\n", paste(fbar, collapse = " "),
    "\n# Last age group considered plus group 0=no 1=yes\n", plusgroup,
    "\n# Number of observations\n", nrow(df.data),
    "\n# Observation data frame",
    "\n# fleet\tyear\tage\tobservation\tweights\n", file=filename); write.t(df.data, file=filename)
  df.aux <- unique(full.df[c("year","age","m","m.spwn","harvest.spwn","mat.wt","stock.wt")])
  cat("# Auxilliary data frame", # should include offsets here?!
    "\n# Number of auxilliary data\n", nrow(df.aux),
    "\n# year\tage\tm\tm.spwn\tharvest.spwn\tmat.wt\n", file=filename, append = TRUE); write.t(df.aux, file=filename)

  #-------------------------------------------------------------------------
  # write config files
  #-------------------------------------------------------------------------

  # fmodel
  filename <- paste0(wkdir,'/fmodel.cfg')

  cat("# F model config for the a4a model",
    "\n# model params\n", ncol(Xf),
    "\n# number of rows\n", nrow(Xf),
    "\n# design matrix\n", file = filename); write.t(Xf, file=filename)

  Covf <- getCov(nrow(Xf), model = "iid", tau = 1)
  cat("# prior model (in sparse format)",
    "\n# flag to turn F-deviations on and off 0=off 1=on\n", 0, # off for now - until we work on the interface for randomness
    "\n# var-cov matrix\n", file = filename, append = TRUE); write.t.sparse(Covf, file=filename)

  # qmodel
  filename <- paste0(wkdir,'/qmodel.cfg')

  cat("# Q model config for the a4a model",
    "\n# model params\n", ncol(Xq),
    "\n# number of rows\n", nrow(Xq),
    "\n# design matrix\n", file = filename); write.t(Xq, file=filename)

  Covq <- getCov(nrow(Xq), model = "iid", tau = 1)
  cat("# prior model (in sparse format)",
    "\n# flag to turn Q-deviations on and off 0=off 1=on\n", 0, # off for now - until we work on the interface for randomness
    "\n# var-cov matrix\n", file = filename, append = TRUE); write.t.sparse(Covq, file=filename)

  # vmodel no random effects for variances
  filename <- paste0(wkdir,'/vmodel.cfg')

  cat("# variance model config for the a4a model",
    "\n# model params\n", ncol(Xv),
    "\n# number of rows\n", nrow(Xv),
    "\n# design matrix\n", file = filename); write.t(Xv, file=filename)

  # n1model no random effects for initial ages
  filename <- paste0(wkdir,'/ny1model.cfg')

  cat("# initial age structure model config for the a4a model",
    "\n# model params\n", ncol(Xny1),
    "\n# number of rows\n", nrow(Xny1),
    "\n# design matrix\n", file = filename); write.t(Xny1, file=filename)

  # rmodel
  filename <- paste0(wkdir,'/srrmodel.cfg')

  cat("# R model config for the a4a model",
    "\n# SR model ID:",srr$srr,"\n", srr$ID,
    "\n# SR CV:\n", srr$srrCV,
    "\n# SPR0 :\n", srr$SPR0,
    "\n# a model params\n", ncol(Xsra),
    "\n# a model number of rows\n", nrow(Xsra),
    "\n# a model design matrix\n", file = filename); write.t(Xsra, file=filename)
  cat("# b model params\n", ncol(Xsrb),
    "\n# b model number of rows\n", nrow(Xsrb),
    "\n# b model design matrix\n", file = filename, append = TRUE); write.t(Xsrb, file=filename)

  Covr <- getCov(nrow(Xsra), model = "iid", tau = 1)
  cat("# prior model for the SRR a param",
    "\n# flag to turn SRR a param deviations on and off 0=off 1=on\n", 0, # off for now - until we work on the interface for randomness
    "\n# var-cov matrix\n", file = filename, append = TRUE); write.t.sparse(Covr, file=filename)

  # r internal model
  filename <- paste0(wkdir,'/rmodel.cfg')

  cat("# internal model for recruits: orthoganal design",
    "\n# model params\n", ncol(Xr),
    "\n# number of rows\n", nrow(Xr),
    "\n# design matrix\n", file = filename); write.t(Xr, file=filename)

   # set up variable names
  pnames <- list(paste0("fMod:",colnames(Xf)),
    paste0("qMod:", unlist(sapply(1:length(indices), function(i) paste0(fleet.names[i+1],":",colnames(Xqlist[[i]]))))),
    paste0("vMod:", unlist(sapply(1:length(fleet.names), function(i) paste0(fleet.names[i],":",colnames(Xvlist[[i]]))))),
    paste0("n1Mod:",colnames(Xny1)),
    paste0("rMod:",colnames(Xr)),
    paste0("sraMod:",colnames(Xsra)),
    paste0("srbMod:",colnames(Xsrb)))

  if (srr$srrCV < 0) pnames <- pnames[-c(6,7)]
  if (srr$ID==4) pnames <- pnames[-7]

  # end here if we just want to write the data and model files
  if (fit == "setup") return(wkdir)


  #========================================================================
  # Run the ADMB executable (build up argument list first)
  #========================================================================

  # run a4a split
  my.time.used[2] <- Sys.time()

  # arguments
  args <- character(0)
  # MCMC
  if (fit == "MCMC") args <- getADMBCallArgs(mcmc)
  # if running an MSE no need to work out hessian
  if (fit == "MP") args <- c(args, "-est")
  args <- paste(args, collapse = " ")

  # run executable in wkdir directory
  if (os.type("linux")) {
    if (verbose) {
      echoc <- system(paste0("cd ", shQuote(wkdir), ";a4a ", args))
    } else {
      echoc <- system(paste0("cd ", shQuote(wkdir), ";a4a ", args, " > logfile.txt"))
    }
    if(fit=="MCMC") system(paste0("cd ", shQuote(wkdir), ";a4a -mceval"))
  } else if (os.type("osx")) {
    if (verbose) {
      echoc <- system(paste0("cd ", shQuote(wkdir), ";a4a ", args))
    } else {
      echoc <- system(paste0("cd ", shQuote(wkdir), ";a4a ", args, " > logfile.txt"))
    }
    if(fit=="MCMC") system(paste0("cd ", shQuote(wkdir), ";a4a -mceval"))
  } else if (os.type("windows")) {
    if (verbose) {
      echoc <- shell(paste0("cd /D", shQuote(wkdir), " & a4a", args))
    } else {
      echoc <- shell(paste0("cd /D", shQuote(wkdir), " & a4a ", args, " > logfile.txt"))
    }
    if(fit=="MCMC") shell(paste0("cd /D", shQuote(wkdir), " & a4a -mceval"))
  }

  #========================================================================
  # Read in ADMB output: how this is done will depend on the 'fit' argument
  #========================================================================

  # post processing split
  my.time.used[3] <- Sys.time()

  if (fit %in% c("MP","assessment", "MCMC")) {

    # read admb output from file
    out <- list()

    # read .par file
    out[c("nopar","nlogl","maxgrad")] <- as.numeric(scan(paste0(wkdir, '/a4a.par'), what = '', nlines = 1, quiet = TRUE)[c(6, 11, 16)])
    lin <- matrix(readLines(paste0(wkdir, '/a4a.par'))[-1], ncol = 2, byrow = TRUE)
    out$par.est <- lapply(strsplit(sub(" ", "",lin[,2]), " "), as.numeric)
    names(out$par.est) <- gsub("[# |:]", "", lin[,1])

    # read derived model quantities
    if(fit=="MP"){
      convergence <- NA
      ages <- sort(unique(full.df$age))
      years <- sort(unique(full.df$year))
      out$N <- as.matrix(read.table(paste0(wkdir, '/n.out'), header = FALSE))
      out$F <- as.matrix(read.table(paste0(wkdir, '/f.out'), header = FALSE))
      out$Q <- as.matrix(read.table(paste0(wkdir, '/q.out'), header = FALSE))
      dim(out$Q) <- c(length(years), length(indices), length(ages))
      out$Q <- aperm(out$Q, c(3,1,2))
        colnames(out$N) <- colnames(out$F) <- ages
        rownames(out$N) <- rownames(out$F) <- years
        dimnames(out$Q) <- list(ages, years, names(indices))
    } else {
    # read derived model quantities and hessian
      if (!file.exists(paste0(wkdir, '/a4a.cor'))) {
        # this is all a bit spaghetti coding to reverse engineering Colin's code
        # the idea is to return a valid object with the correct dimensions
        # with NA for N and F but with the parameters estimated by ADMB
        # although if cor is not present means that it didn't converge ...
        warning("Hessian was not positive definite.")
        convergence <- 1
        ages <- sort(unique(full.df$age))
        years <- sort(unique(full.df$year))
        rr <- sum(!(out$par.est$rapar>0))==length(out$par.est$rapar)
        out$par.est <- lapply(out$par.est, function(x) x[] <- rep(NA, length(x)))
        out$par.std <- out$par.est
        if(isTRUE(rr)) out$par.std$rapar <- out$par.std$rbpar <- vector(mode="numeric")
        out$cov <- array(NA, dim=rep(length(unlist(pnames)), 2))
        out$prec <- array(NA, dim=rep(length(unlist(pnames)), 2))
        out$nopar <- NA
        out$nlogl <- NA
        out$maxgrad <- NA
        out$cgcv <- NA
        flqNA <- stock.n(stock)
        flqNA[] <- NA
        out$N <- t(flqNA[drop=T])
        out$F <- t(flqNA[drop=T])
        out$Q <- array(NA, dim=c(length(ages), length(years), length(indices)), dimnames=list(ages, years, names(indices)))
      } else {
        convergence <- 0
        # read .cor file
        lin <- readLines(paste0(wkdir, '/a4a.cor'))
        npar <- length(lin) - 2
        out$logDetHess <- as.numeric(strsplit(lin[1], '=')[[1]][2])

        sublin <- lapply(strsplit(lin[1:npar + 2], ' '), function(x) x[x != ''])
        par.names <- unlist(lapply(sublin, function(x) x[2]))
        par.std <- as.numeric(unlist(lapply(sublin, function(x) x[4])))

        out$par.std <- lapply(names(out$par.est), function(x) par.std[which(par.names==x)])
        names(out$par.std) <- names(out$par.est)

        # use this as it seems to be more robust.. strangely...
        # I think it is because the solve method that ADMB uses is not
        # as good as the R one.... small numerical errors are
        # resulting in non-positive definite vcov mats for subsets of parameters.
        #hess <- getADMBHessian(wkdir)$hes
        #out$cov <- solve(hess, tol=1e-50)
        #out$prec <- hess
        #out$nopar <- ncol(hess)

        out$cov <- getADMBCovariance(wkdir)$cov
        out$nopar <- ncol(out$cov)

        # read derived model quantities
        ages <- sort(unique(full.df$age))
        years <- sort(unique(full.df$year))

        out$N <- as.matrix(read.table(paste0(wkdir, '/n.out'), header = FALSE))
        out$F <- as.matrix(read.table(paste0(wkdir, '/f.out'), header = FALSE))
        out$Q <- as.matrix(read.table(paste0(wkdir, '/q.out'), header = FALSE))
        dim(out$Q) <- c(length(years), length(indices), length(ages))
        out$Q <- aperm(out$Q, c(3,1,2))
          colnames(out$N) <- colnames(out$F) <- ages
          rownames(out$N) <- rownames(out$F) <- years
          dimnames(out$Q) <- list(ages, years, names(indices))
      }
    }
  }

  if (fit == "MCMC") {

    # get number of iters (CHECK WITH NITERS  which are not being used )
    nit <- getN(mcmc)

    # read files
    filen <- file(paste0(wkdir, '/a4a.psv'), "rb")
    nopar <- readBin(filen, what = integer(), n = 1)
    mcmcout <- readBin(filen, what = numeric(), n = nopar * nit)# TODO check this line
    close(filen)
    out$mcmc$mcmcout <- matrix(mcmcout, byrow = TRUE, ncol = nopar)
    colnames(out$mcmc$mcmcout) <- unlist(pnames)

    out$mcmc <- list()

    out$mcmc$nit <- nit
    out$mcmc$N <- read.table(paste0(wkdir, "/NMCMCreport.csv"), sep=" ", header=FALSE)[-1]
    out$mcmc$F <- read.table(paste0(wkdir, "/FMCMCreport.csv"), sep=" ", header=FALSE)[-1]
    out$mcmc$QQ <- read.table(paste0(wkdir, "/QMCMCreport.csv"), sep=" ", header=FALSE)[-1]
    out$mcmc$idq <- expand.grid(y=years, s=fleet.names[-1], i=1:nit)
    out$mcmc$V <- read.table(paste0(wkdir, "/VMCMCreport.csv"), sep=" ", header=FALSE)[-1]
    out$mcmc$idv <- expand.grid(y=years, s=fleet.names, i=1:nit)
    out$mcmc$logfile <- readLines(paste0(wkdir, "/logfile.txt"))
  }

  # return
  list(out = out, my.time.used = my.time.used, wkdir = wkdir,
       convergence = convergence, pnames = pnames,
       ages = ages, years = years)
}

