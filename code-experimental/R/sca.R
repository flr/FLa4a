a4a <- function(…) {

  message(“The a4a function has been renamed: please use the function a4aSCA in place of a4a. “) # or something

}


#' @title statistical catch-at-age method
#' @name sca
#' @docType methods
#' @rdname sca
#' @description User interface to the statistical catch-at-age method of the a4a stock assessment framework.   
#' @param stock a \code{FLStock} object
#' @param indices a \code{FLIndices} object
#' @param fmodel a formula object depicting the model for log fishing mortality at age
#' @param qmodel a list of formula objects depicting the models for log survey catchability at age
#' @param srmodel a formula object depicting the model for log recruitment
#' @param fit Character with type of fit, 'MP' or 'assessment', the first doesn't require the hessian to be computed, while the former does.
#' @return a \code{a4aFit} or \code{a4aFitSA} object with the fit results. 
#' @aliases sca sca-methods sca,FLStock,FLIndices-method sca,FLStock,FLIndex-method
setGeneric("sca", function(stock, indices, ...) standardGeneric("sca"))
setMethod("sca", signature("FLStock", "FLIndices"), function(stock, indices, fmodel  = ~ factor(age) + factor(year), qmodel  = lapply(seq(length(indices)), function(i) ~ s(age, k=3)), srmodel = ~ factor(year), fit = "MP")
{
	a4aSCA(stock=stock, indices=indices, fmodel=fmodel, qmodel=qmodel, srmodel=srmodel, fit=fit)

})

setMethod("sca", signature("FLStock", "FLIndex"), function(stock, indices, fmodel  = ~ factor(age) + factor(year), qmodel  = lapply(seq(length(indices)), function(i) ~ s(age, k=3)), srmodel = ~ factor(year), fit = "MP")
{
	a4aSCA(stock=stock, indices=FLIndices(indices), fmodel=fmodel, qmodel=qmodel, srmodel=srmodel, fit=fit)

})

#' @title statistical catch-at-age advanced method
#' @name a4aSCA
#' @docType methods
#' @rdname a4aSCA
#' @description Advanced user interface to the statistical catch-at-age method of the a4a stock assessment framework.   
#'
#' @param fmodel a formula object depicting the model for log fishing mortality at age
#' @param qmodel a list of formula objects depicting the models for log survey catchability at age
#' @param srmodel a formula object depicting the model for log recruitment
#' @param stock an FLStock object containing catch and stock information
#' @param indices an FLIndices object containing survey indices 
#' @param covar a list with covariates 
#' @param wkdir used to set a working directory for the admb optimiser.  If wkdir is set all admb files are saved to this folder otherwise they are deleted.
#' @param verbose if true admb fitting information is printed to the screen
#' @param fit Character with type of fit, 'MP' or 'assessment', the first doesn't require the hessian to be computed, while the former does.
#' @return an \code{a4aFit} object if fit is "MP" or an \code{a4aFitSA} if fit is "assessment"
#' @aliases a4aSCA a4aSCA-methods a4aSCA,FLStock,FLIndices-method
#' @template Example-a4a
setGeneric("a4aSCA", function(stock, indices, ...) standardGeneric("a4aSCA"))
setMethod("a4aSCA", signature("FLStock", "FLIndices"), function(stock, indices, fmodel  = ~ s(age, k = 3) + factor(year), qmodel  = lapply(seq(length(indices)), function(i) ~ 1), srmodel = ~ factor(year), n1model = ~factor(age), vmodel  = lapply(seq(length(indices) + 1), function(i) ~ 1), covar = NULL, wkdir = NULL, verbose = FALSE, fit = "assessment", center = TRUE) {

	# NOTE: what is niters doing in a4aInternal ??

  fit <- match.arg(fit, c("MP", "assessment", "sim", "Ext"))

  # now to deal with iterations ...
  
  # create a df for dimension information:
  dms <- do.call(rbind.data.frame, c(list(catch = c(dims(stock), startf = NA, endf = NA)), lapply(indices, dims)))

  # average stock over seasons
  # NOTE: Do we have a warning msg about this ?
  stock <- collapseSeasons(stock)
  
  # only allow 1 season for surveys
  if (any(dms $ season[-1] > 1)) stop("only one season per survey - please split into seperate surveys.")
  
  # now do a fit for each combination of unit, area and iter...
  # if fit = MP then we return an a4aFit with the same dimensions as stock
  # if fit = assessment then we return a4aFitSA with same dimensions as stock....  \TODO only true with iters so far
  
  grid <- do.call(expand.grid, c(dimnames(catch.n(stock))[c(3,5)], list(iter = 1:max(dms $ iter))))
  #if (!identical(sort(unique(dms $ iter)), sort(unique(c(1L, max(dms $ iter))))))
  if(length(unique(dms$iter[dms$iter>1]))>1) 
  	stop("incosistent number of iterations in stock and indices")
  it <- max(dms$iter)

  # set up objects
  # stk
  dms <- dimnames(stock.n(stock))
  dms$iter <- 1:it
  ini <- FLQuant(NA, dimnames=dms)
  out <- if (fit %in% c("MP", "sim")) a4aFit() else a4aFitSA()
  out @ desc <- desc(stock)
  out @ name <- name(stock)
  out @ range <- range(stock)
  out @ call <- match.call()
  out @ harvest <- ini
  out @ stock.n <- ini
  out @ catch.n <- ini
  # idx
  ini <- lapply(indices, function(x){
  	dms <- dimnames(index(x))
  	dms$iter <- 1:it
	FLQuant(NA, dimnames=dms)
  })
  #out @ index <- FLQuants(lapply(indices, index))
  out @ index <- FLQuants(ini)

  if (fit == "assessment") {
    out @ pars @ stkmodel @ fMod <- fmodel
    out @ pars @ stkmodel @ n1Mod <- n1model
    out @ pars @ stkmodel @ srMod <- srmodel
    # and the same for indices    
  }

  time.used <- matrix(NA, nrow = 4, ncol = nrow(grid))
  ifit <- if (fit == "sim") "assessment" else fit

  niters <- nrow(grid)
  for (i in seq(niters)) {
    # subset the stock
    istock <- stock[,, grid $ unit[i], grid $ area[i], , min(grid $ iter[i], dims(stock)$iter)]

	# check: do we need indices to have matching units, areas?
    iindices <- lapply(indices, function(x) x[,, grid $ unit[i], grid $ area[i], , min(grid $ iter[i], dims(x)$iter)])
    iindices <- FLIndices(iindices)
    # check: do we need indices to have matching units, areas?
    if (!is.null(covar)) {
      icovar <- lapply(covar, function(x) x[,, grid $ unit[i], grid $ area[i], , min(grid $ iter[i], dims(x)$iter)])
    } else {
      icovar <- NULL
    }

    # run model      
    outi <- a4aInternal(fmodel = fmodel, qmodel = qmodel, srmodel = srmodel,
                        n1model = n1model, vmodel = vmodel,
                        stock = istock, indices = iindices, covar = icovar, 
                        wkdir = wkdir, verbose = verbose, 
                        fit = ifit, center = center)

    if (i == 1) {
      tmpSumm <- outi @ fitSumm
      out @ fitSumm <- array(0, c(dim(tmpSumm), niters), c(dimnames(tmpSumm), list(iters = 1:niters)))
    }
    out @ fitSumm[,i] <- outi @ fitSumm
      
    if (fit == "MP") {
      # copy results
      out @ harvest[,, grid $ unit[i], grid $ area[i], , grid $ iter[i]] <- harvest(outi)
      out @ stock.n[,, grid $ unit[i], grid $ area[i], , grid $ iter[i]] <- stock.n(outi)
      out @ catch.n[,, grid $ unit[i], grid $ area[i], , grid $ iter[i]] <- catch.n(outi)
      # add indices
      for (j in 1:length(iindices)) {
        out @ index[[j]][,, grid $ unit[i], grid $ area[i], , grid $ iter[i]] <- index(outi)[[j]]
      }
    }
    
    if (fit == "sim") {

      # copy results with noise
      istock <- istock + outi # this automatically adds noise
      out @ harvest[,, grid $ unit[i], grid $ area[i], , grid $ iter[i]] <- harvest(istock)
      out @ stock.n[,, grid $ unit[i], grid $ area[i], , grid $ iter[i]] <- stock.n(istock)
      out @ catch.n[,, grid $ unit[i], grid $ area[i], , grid $ iter[i]] <- catch.n(istock)
      # add indices
      for (j in 1:length(iindices)) {
        out @ index[[j]][,, grid $ unit[i], grid $ area[i], , grid $ iter[i]] <- index(outi)[[j]]
      }

    }

    if (fit %in% c("assessment", "Ext")) {

      # store everything in a a4aFit SA object
      out @ harvest[,, grid $ unit[i], grid $ area[i], , grid $ iter[i]] <- harvest(outi)
      out @ stock.n[,, grid $ unit[i], grid $ area[i], , grid $ iter[i]] <- stock.n(outi)
      out @ catch.n[,, grid $ unit[i], grid $ area[i], , grid $ iter[i]] <- catch.n(outi)
      # add indices
      for (j in 1:length(iindices)) {
        out @ index[[j]][,, grid $ unit[i], grid $ area[i], , grid $ iter[i]] <- index(outi)[[j]]
      }

      # fill up models
      if (i == 1) {
        # stkmodel
        out @ pars @ stkmodel @ name      <- outi @ name
        out @ pars @ stkmodel @ desc      <- outi @ desc
        out @ pars @ stkmodel @ range     <- outi @ range
        out @ pars @ stkmodel @ centering <- rep(0, niters)  
        out @ pars @ stkmodel @ params    <- propagate(outi @ pars @ stkmodel @ params, niters)
        tmpvcov <- outi @ pars @ stkmodel @ vcov
        out @ pars @ stkmodel @ vcov      <- array(0, c(dim(tmpvcov), niters), c(dimnames(tmpvcov), list(iters = 1:niters)))
        out @ pars @ stkmodel @ m         <- propagate(outi @ pars @ stkmodel @ m, niters)
        out @ pars @ stkmodel @ units     <- units(catch.n(stock))
        # qmodel
        out @ pars @ qmodel               <- outi @ pars @ qmodel
        for (j in seq(length(indices))) {
          out @ pars @ qmodel[[j]] @ params <- propagate(outi @ pars @ qmodel[[j]] @ params, niters)
          tmpvcov <- outi @ pars @ qmodel[[j]] @ vcov
          out @ pars @ qmodel[[j]] @ vcov        <- array(0, c(dim(tmpvcov), niters), c(dimnames(tmpvcov), list(iters = 1:niters)))
        }
        # vmodel
        out @ pars @ vmodel               <- outi @ pars @ vmodel
        for (j in seq(length(indices)+1)) {
          out @ pars @ vmodel[[j]] @ params <- propagate(outi @ pars @ vmodel[[j]] @ params, niters)
          tmpvcov <- outi @ pars @ vmodel[[j]] @ vcov
          out @ pars @ vmodel[[j]] @ vcov        <- array(0, c(dim(tmpvcov), niters), c(dimnames(tmpvcov), list(iters = 1:niters)))
        }
      }
    
      # now the a4aFitSA bits                                 
      out @ pars @ stkmodel @ centering[i] <- outi @ pars @ stkmodel @ centering     
      out @ pars @ stkmodel @ params[,i]   <- outi @ pars @ stkmodel @ params
      out @ pars @ stkmodel @ vcov[,,i]    <- outi @ pars @ stkmodel @ vcov
      out @ pars @ stkmodel @ m[,,,,,i]    <- outi @ pars @ stkmodel @ m
      # qmodel
      for (j in seq(length(indices))) {
        out @ pars @ qmodel[[j]] @ params[,i] <- outi @ pars @ qmodel[[j]] @ params
        out @ pars @ qmodel[[j]] @ vcov[,,i]  <- outi @ pars @ qmodel[[j]] @ vcov
      }
      # vmodel
      for (j in seq(length(indices)+1)) {
        out @ pars @ vmodel[[j]] @ params[,i] <- outi @ pars @ vmodel[[j]] @ params
        out @ pars @ vmodel[[j]] @ vcov[,,i]   <- outi @ pars @ vmodel[[j]] @ vcov
      }

    }
    
    if (fit == "Ext") {
      ### ??
    }    
    
    # keep timing info
    time.used[,i] <- outi @ clock
  }
   
  units(out @ harvest) <- "f"  
   
  # add in combined timings
  out @ clock <- outi @ clock # to get names
  out @ clock[] <- rowSums(time.used)
  
  # return out
  out
})

