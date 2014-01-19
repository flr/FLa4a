
#' The user interface to the a4a fitting routine.
#'
#'
#' @param fmodel a formula object depicting the model for log fishing mortality at age
#' @param qmodel a list of formula objects depicting the models for log survey catchability at age
#' @param rmodel a formula object depicting the model for log recruitment
#' @param stock an FLStock object containing catch and stock information
#' @param indices an FLIndices object containing survey indices 
#' @param fmodel.extra a named list of expressions used to add covariates to the model for log fishing mortality at age
#' @param qmodel.extra a named list of expressions used to add covariates the models for log survey catchability at age
#' @param wkdir used to set a working directory for the admb optimiser.  If wkdir is set all admb files are saved to this folder otherwise they are deleted.
#' @param verbose if true admb fitting information is printed to the screen
#' @return an FLa4aFit object
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @author Ernesto Jardim \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
collapseSeasons <- function (stock) {
  
  if (dims(stock) $ season == 1) return (stock) # do nothing
  
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


#' The user interface to the a4a fitting routine.
#'
#'
#' @param fmodel a formula object depicting the model for log fishing mortality at age
#' @param qmodel a list of formula objects depicting the models for log survey catchability at age
#' @param rmodel a formula object depicting the model for log recruitment
#' @param stock an FLStock object containing catch and stock information
#' @param indices an FLIndices object containing survey indices 
#' @param fmodel.extra a named list of expressions used to add covariates to the model for log fishing mortality at age
#' @param qmodel.extra a named list of expressions used to add covariates the models for log survey catchability at age
#' @param wkdir used to set a working directory for the admb optimiser.  If wkdir is set all admb files are saved to this folder otherwise they are deleted.
#' @param verbose if true admb fitting information is printed to the screen
#' @return an FLa4aFit object
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @author Ernesto Jardim \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
a4a <- function(fmodel  = ~ s(age, k = 3) + factor(year), 
                qmodel  = lapply(seq(length(indices)), function(i) ~ 1), 
                srmodel = ~ factor(year),
                stock, indices, covar = NULL, 
                wkdir = NULL, verbose = FALSE, fit = "MP")
{

  fit <- match.arg(fit, c("MP", "assessment", "sim", "Ext"))

  # hard coded settings
  n1model <- ~ factor(age)
  vmodel  <- lapply(seq(length(indices) + 1), function(i) ~ 1)
  vmodel[[1]] <- ~ s(age, k = 3)
  
  #
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
#  out @ harvest <- harvest(stock)[,,,,,rep(1,it)]
#  out @ stock.n <- stock.n(stock)
#  out @ catch.n <- catch.n(stock)
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
    #istock <- stock[,, grid $ unit[i], grid $ area[i], , grid $ iter[i]]
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
                        fit = ifit)

    if (i == 1) {
      tmpSumm <- outi @ fitSumm
      out @ fitSumm                     <- array(0, c(dim(tmpSumm), niters), c(dimnames(tmpSumm), list(iters = 1:niters)))
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
}


#' The user interface to the a4a fitting routine.
#'
#'
#' @param fmodel a formula object depicting the model for log fishing mortality at age
#' @param qmodel a list of formula objects depicting the models for log survey catchability at age
#' @param rmodel a formula object depicting the model for log recruitment
#' @param stock an FLStock object containing catch and stock information
#' @param indices an FLIndices object containing survey indices 
#' @param fmodel.extra a named list of expressions used to add covariates to the model for log fishing mortality at age
#' @param qmodel.extra a named list of expressions used to add covariates the models for log survey catchability at age
#' @param wkdir used to set a working directory for the admb optimiser.  If wkdir is set all admb files are saved to this folder otherwise they are deleted.
#' @param verbose if true admb fitting information is printed to the screen
#' @return an FLa4aFit object
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @author Ernesto Jardim \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
#' @examples
#' data(ple4)
#' data(ple4.indices)
#'
#' # define sub models
#' fmodel <- ~ s(age, k=4) + factor(year)
#' # three surveys - try 'length(ple4.indices)' - so three formulas
#' qmodel <- list(~ s(age, k=4), ~ s(age, k=4), ~ age)
#' fit <- a4a(fmodel = fmodel, qmodel = qmodel, stock = ple4, indices = ple4.indices, fit = "assessment")
#'
#' # add the fit to the stock, but with fitting error
#' fitstk <- propagate(ple4, 1000) + fit
#' 
#' doOne <- function(p, object) cbind(as.data.frame(apply(object, 1:5, quantile, p)), p = p[1])
#' dat <- do.call(rbind, lapply(c(0.025, 0.5, 0.975), doOne, object = ssb(fitstk)))
#' xyplot(I(data/1000) ~ year, group = p, data = dat, 
#'        type = c("l","g"), lwd = c(1,2,1), lty = c(2,1,2), col = 1,
#'        ylab = "SSB (thousand tonnes)", xlab = "Year")
#' # end of example
a4aInternal <- function(fmodel  = ~ s(age, k = 3) + factor(year), 
                qmodel  = lapply(seq(length(indices)), function(i) ~ 1), 
                srmodel = ~ factor(year),
                n1model = ~ factor(age), 
                vmodel  = lapply(seq(length(indices) + 1), function(i) ~ 1),
                stock, indices, covar = NULL, 
                wkdir = NULL, verbose = FALSE, fit = "assessment", 
                niters = 1000, center = TRUE)
{


  # first check permissions of executable
  exeok <- check.executable()
  if (!exeok) stop("a4a executable has wrong permissions.")

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
    snames <- make.unique(rep("survey", length(indices)))
  } else {
    snames <- make.unique(names(indices))
  }
  for (i in seq_along(indices)) name(indices[[i]]) <- snames[i]
  names(indices) <- snames
 

  # what kind of run is this: 'setup' just writes data files - usefull when developing code
  fit <- match.arg(fit, c("MP", "assessment", "debug", "setup", "MCMC", "Ext")) # MCMC is experimental


  # ------------------------------------------------------------------------
  #
  # Extract observations and auxilliary info from stock and indices objects
  #
  # ------------------------------------------------------------------------
  
  # first some checks
  if (any(is.infinite(log(catch.n(stock))))) stop("only non-zero catches allowed.")
  if (any(is.infinite(log( unlist(lapply(indices, function(x) c(index(x)))) ))))  stop("only non-zero survey indices allowed.")
  if (any(is.infinite(log( unlist(lapply(indices, function(x) c(index.var(x)))) ))))  stop("only non-zero survey index variances allowed.")
   
  # utility to convert to a 2d array
  quant2mat <- 
  function(x) {
    out <- x[drop=TRUE]
    dim(out) <- dim(x)[1:2]
    dimnames(out) <- dimnames(x)[1:2]
    if (nrow(out) == 1 && dimnames(out)[[1]] == "all") dimnames(out)[[1]] <- NA_character_  # "all" denotes a biomass survey
    out 
  }

  # convert catches and indices to a list of named arrays
  list.obs <- c(list(catch = quant2mat(catch.n(stock))),
                     lapply(indices, function(x) quant2mat(index(x)) ))
                     
  # convert the variances of catches and indices to a list of named arrays
  list.var <- c(list(catch = quant2mat(catch.n(stock) @ var)),
                     lapply(indices, function(x) quant2mat(index.var(x)) ))                
                     
  # calculate appropriate centering for observations on log scale
  center.log <- sapply(list.obs, function(x) mean(log(x), na.rm = TRUE))
  if (!center) center.log[] <- 0

  # convert to dataframe
  list2df <- function(fleet)
  {
    x <- list.obs[[fleet]]
    v <- as.vector(list.var[[fleet]])
    year <- as.numeric(colnames(x)[col(x)])
    age <- as.numeric(rownames(x)[row(x)])
    obs <- log(as.vector(x)) - center.log[fleet]
    if (all(is.na(v))) {
      wts <- 1
    } else # use inverse variance weighting
    {
      wts <-  1 / v # inverse variance weigting
    }
    ret <- data.frame(fleet = fleet, year = year, age = age, obs = obs, weights = wts)
    ret <- ret[!is.na(ret $ obs), ]
    if (any(is.na(ret[,5])) || any(ret[,5] <= 0)) {
      ret[,5] <- 1
      warning("*** NA and/or non-positive variances found in: ", names(list.obs)[fleet], " - all variances set to 1", call. = FALSE)
    }
    ret
  }

  df.data <- do.call(rbind, lapply(1:length(list.obs), list2df))
  
  if (any(df.data[,5] != 1)) message("Note: Provided variances will be used to weight observations.\n\tWeighting assumes variances are on the log scale or equivalently log(CV^2 + 1).")

  # extract auxilliary stock info
  fbar <-  unname(range(stock)[c("minfbar","maxfbar")])
  plusgroup <- as.integer( !is.na(range(stock)["plusgroup"]), range(stock)["plusgroup"] >= range(stock)["max"] )

  # extract auxilliary survey info - always assume oldest age is true age TODO TODO TODO !! 
  surveytime <- unname(sapply(indices, function(x) mean(c(dims(x) $ startf, dims(x) $ endf))))
  if (any(is.na(surveytime))) stop("You need to define startf and endf for each index!!")

  # ------------------------------------------------------------------------
  #
  # Make a full data.frame and add in covariates and observations
  #
  # ------------------------------------------------------------------------
 
  # build a full data frame first (we will use this for the variance model so it is not a waste)
  make.df <- function(fleet) {
    thing <- if (fleet == 1) stock else indices[[fleet - 1]]
    expand.grid(age = if (is.na(range(thing)["min"])) NA else range(thing)["min"]:range(thing)["max"], 
                year = range(thing)["minyear"]:range(thing)["maxyear"])[2:1]
  }
    
  full.df <- do.call(rbind, lapply(1:length(list.obs), function(i) cbind(fleet = i, make.df(i))))

  if (!is.null(covar)) {
  # add in covariates to data.frame - it is easiest to provide covariates in one list
  tmp <- 
    lapply(seq_along(covar), 
      function(i) {
        x <- as.data.frame(covar[[i]])[c(1,2,7)]
        if (length(unique(x $ age)) == 1) x <- x[names(x) != "age"]
        if (length(unique(x $ year)) == 1) x <- x[names(x) != "year"]
        names(x) <- gsub("data", names(covar)[i], names(x))
        x
      })
  covar.df <- tmp[[1]]
  for (i in seq_along(tmp[-1])) covar.df <- merge(covar.df, i, all = TRUE, sort = FALSE)

  full.df <- merge(full.df, covar.df, all.x = TRUE, all.y = FALSE)
  } 
   
  # add in data
  full.df <- merge(full.df, df.data, all.x = TRUE, all.y = FALSE)
  # put biomass surveys in min age position
  full.df $ age <- with(full.df, replace(age, is.na(age), min(age, na.rm = TRUE)))

  full.df $ fleet <- factor(names(list.obs)[full.df $ fleet], levels = names(list.obs))

  # set weights for missing values to zero
  full.df $ weights[is.na(full.df $ obs)] <- 0
  # inform that missing values will be treated as missing at random
  if (any(is.na(full.df $ obs))) 
    message("Note: The following observations are treated as being missing at random:\n\t", 
             paste(capture.output(print(subset(full.df,is.na(obs))[c("fleet","year","age")], row.names = FALSE)), collapse = "\n\t"),
          "\n      Predictions will be made for missing observations." )
  
  # fill out data frame for all eventualities ... except ages .... 
  # or a quick one that would not mess array dims and things later is to fix ages in data frame but keep rows...
  # it would actually be easiest to present all auxilliary data and covariates as a data.frame and the observations as a matrix...
  temp.full.df <- expand.grid(lapply(full.df[3:1],function(x) sort(unique(x))),  KEEP.OUT.ATTRS = FALSE)[3:1]
  for (i in seq_along(indices)) {
    # if biomass survey skip this step
    if (is.na(range(indices[[i]])["min"])) next
    .ages <- temp.full.df $ age [temp.full.df $ fleet == levels(full.df $ fleet)[i+1]]
    .range <- range(subset(full.df, fleet == levels(full.df $ fleet)[i+1]) $ age)
    .ages[.ages < .range[1]] <- .range[1]
    .ages[.ages > .range[2]] <- .range[2] 
    temp.full.df $ age [temp.full.df $ fleet == levels(full.df $ fleet)[i+1]] <- .ages 
  }
  full.df <- merge(temp.full.df, full.df, all.x = TRUE)
  

  # add in auxilliary data - maturity, natural mortality etc.
  aux.df <-
    cbind(as.data.frame(stock.n(stock))[c("year","age")], 
          m = log(as.data.frame(m(stock)) $ data), mat = as.data.frame(mat(stock)) $ data,
          stock.wt = as.data.frame(stock.wt(stock)) $ data, catch.wt = as.data.frame(catch.wt(stock)) $ data,
          m.spwn = as.data.frame(m.spwn(stock)) $ data, harvest.spwn = as.data.frame(harvest.spwn(stock)) $ data)
  
  full.df <- merge(full.df, aux.df, all.x = TRUE, all.y = FALSE)
  full.df $ mat.wt <- full.df $ stock.wt * full.df $ mat 

  # set unspecified covariates and aux data to 0
  # TODO we should warn about this... or check before the merge that there are covars and m and mat for requested predictions
  full.df[is.na(full.df)] <- 0 

  # order df
  full.df <- full.df[with(full.df, order(fleet, year, age)), ]
  rownames(full.df) <- NULL
  full.df <- full.df[c(3,1,2,4:ncol(full.df))]
   
  # TODO stop if some obs have no parameter - it should be quite rare ...  

  # ------------------------------------------------------------------------
  #
  # Process model formulas
  #
  # ------------------------------------------------------------------------


  # make sure contrasts are set to sumto zero for better performance
  opts <- options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly")) 

  # F model matrix
  Xf <- getX(fmodel, subset(full.df, fleet == "catch"))
  # F model offsets
  # ...
  
  # Q model matrix  
  fleet.names <- c("catch", names(indices))
  Xqlist <- lapply(seq_along(indices), function(i) getX(qmodel[[i]], subset(full.df, fleet == fleet.names[i+1])))
  Xq <- as.matrix(do.call(bdiag, Xqlist))  
  # if model is one formula:
  #Xq <- getX(qmodel, subset(full.df, fleet != "catch"))
  # Q model offsets
  # ...
  
  # var model matrix
  Xvlist <- lapply(1:length(fleet.names), function(i) getX(vmodel[[i]], subset(full.df, fleet == fleet.names[i])))
  Xv <- as.matrix(do.call(bdiag, Xvlist))  
  # if model is one formula:  
  #Xv <- getX(vmodel, full.df)
  # var model offsets
  # ...
    
  # initial age structure model matrix
  Xny1 <- getX(n1model, subset(full.df, year == min(year) & age > min(age) & fleet == "catch"))
  # var model offsets
  # ...  
    
  # process recruitment formula: 
  # builder functions - could be more hadley...
  bevholt <- function(a = ~ 1, b = ~ 1, CV = 0.5) {
    if (CV <= 0) stop ("CV in stock recruit relationship cannot be les than zero")
    list(srr = "bevholt", a = a, b = b, SPR0 = 1, srrCV = CV, ID = 1)
  }
  bevholtSV <- function(h = ~ 1, v = ~ 1, SPR0 = 1, CV = 0.5) {
    if (CV <= 0) stop ("CV in stock recruit relationship cannot be les than zero")
    list(srr = "bevholtSV", a = h, b = v, SPR0 = SPR0, srrCV = CV, ID = 5)
  }
  ricker <- function(a = ~ 1, b = ~ 1, CV = 0.5) {
    if (CV <= 0) stop ("CV in stock recruit relationship cannot be les than zero")
    list(srr = "ricker", a = a, b = b, SPR0 = 1, srrCV = CV, ID = 2)
  }
  hockey <- function(a = ~ 1, b = ~ 1, CV = 0.5) {
    if (CV <= 0) stop ("CV in stock recruit relationship cannot be les than zero")
    list(srr = "hockey", a = a, b = b, SPR0 = 1, srrCV = CV, ID = 3)
  }
  geomean <- function(a = ~ 1, CV = 0.5) {
    if (CV <= 0) stop ("CV in stock recruit relationship cannot be les than zero")
    list(srr = "geomean", a = a, b = ~ 1, SPR0 = 1, srrCV = CV, ID = 4)
  }
  none <- function() list(srr = "geomean", a = ~ 1, b = ~ 1, srrCV = -1, ID = 4)

  # now separate model elements, find SR models, stop if more than one specified
  facs <- strsplit(as.character(srmodel)[length(srmodel)], "[+]")[[1]]
  facs <- gsub("(^ )|( $)", "", facs) # remove leading and trailing spaces
  a4as <- grepl(paste("(^",c("bevholt", "bevholtSV", "ricker","hockey","geomean"),"[(])", collapse = "|", sep = ""), facs)
  if (sum(a4as) > 1) stop("you can only specify one type of stock recruit relationship.")
  srrmod <- if (sum(a4as) == 0) "none()" else facs[a4as]
  if (sum(a4as) == 0 && max(full.df $ year) > max(df.data $ year)) stop("you need to specify a stock recruitment relationship to forecast with out survey information.")
    
  # extract a and b model formulas and add on any extra bits to amodel. 
  # NB SRR models should be parametrised so that amodel is the level of recruitment!!!
  srr <- eval(parse(text = srrmod))
  if (sum(a4as) > 0 && any(!a4as)) {
    # ignore .. the following line adds these onto the end of amod
    message("Note: Trailing formula elements in the srmodel have been removed")
    #srr $ amodel <- eval(parse(text = paste("~", as.character(srr $ amodel)[length(srr$amodel)], "+", paste(facs[!a4as], collapse = " + ")) ))
  }
  Xsra <- getX(srr $ a, subset(full.df, fleet == "catch" & age == dims(stock)$min))  
  Xsrb <- getX(srr $ b, subset(full.df, fleet == "catch" & age == dims(stock)$min))  
 
  # can we do a quick check for identifiability of the srr model? ...
  if (ncol(Xsra) + ncol(Xsrb) > dims(stock)$year) stop("Stock recruitment model is over parameterised, please reduce the number parameters")

  # internal r model matrix
  if (sum(a4as) == 0) rmodel <- srmodel else rmodel <- ~ factor(year) 
  Xr <- getX(rmodel, subset(full.df, age == min(age) & fleet == "catch"))

 
  # reset contrast options
  options(opts)

  # ------------------------------------------------------------------------
  #
  # Create the directory where to store model config, data and results files
  #
  # ------------------------------------------------------------------------
  
  if (is.null(wkdir)) keep <- FALSE else keep <- TRUE # keep results if wkdir is set by user
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

  # ------------------------------------------------------------------------
  #
  # Write model matrices and model info to files in wkdir
  #
  # ------------------------------------------------------------------------

  # local utility
  write.t <- function(x, ...) write.table(x, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t', file = filename, append = TRUE)
  write.t.sparse <- function(x, ...) {
     x <- as(x, "dsCMatrix")
     cat("\n# i\n", x @ i, 
         "\n# p\n", x @ p,
         "\n# x\n", x @ x, file = filename, append = TRUE)  
  }  

  # write to data file 
  # ------------------
   
  filename <- paste0(wkdir,'/a4a.dat')

  # change NA to -1 for admb
  df.data $ age <- with(df.data, replace(age, is.na(age), -1))

  cat("# Data for the a4a model",
    "\n# Full age range\n", range(full.df $ age),
    "\n# Full year range\n", range(full.df $ year),
    "\n# Number of surveys\n", length(unique(full.df $ fleet)) - 1,
    "\n# Survey time as a fraction into the year (one for each survey)\n", paste(surveytime, collapse = " "), 
    "\n# fbar range\n", paste(fbar, collapse = " "),
    "\n# Last age group considered plus group 0=no 1=yes\n", plusgroup,
    "\n# Number of observations\n", nrow(df.data),
    "\n# Observation data frame",
    "\n# fleet\tyear\tage\tobservation\tweights\n", file=filename); write.t(df.data)
  df.aux <- unique(full.df[c("year","age","m","m.spwn","harvest.spwn","mat.wt","stock.wt")])
  cat("# Auxilliary data frame", # should include offsets here?!
    "\n# Number of auxilliary data\n", nrow(df.aux),
    "\n# year\tage\tm\tm.spwn\tharvest.spwn\tmat.wt\n", file=filename, append = TRUE); write.t(df.aux)
  
  # write config files
  # ------------------
  
  # fmodel
  filename <- paste0(wkdir,'/fmodel.cfg')

  cat("# F model config for the a4a model",
    "\n# model params\n", ncol(Xf), 
    "\n# number of rows\n", nrow(Xf),
    "\n# design matrix\n", file = filename); write.t(Xf)

  Covf <- getCov(nrow(Xf), model = "iid", tau = 1)
  cat("# prior model (in sparse format)",
    "\n# flag to turn F-deviations on and off 0=off 1=on\n", 0, # off for now - until we work on the interface for randomness
    "\n# var-cov matrix\n", file = filename, append = TRUE); write.t.sparse(Covf)
 
  # qmodel
  filename <- paste0(wkdir,'/qmodel.cfg')

  cat("# Q model config for the a4a model",
    "\n# model params\n", ncol(Xq), 
    "\n# number of rows\n", nrow(Xq),
    "\n# design matrix\n", file = filename); write.t(Xq)

  Covq <- getCov(nrow(Xq), model = "iid", tau = 1)
  cat("# prior model (in sparse format)",
    "\n# flag to turn Q-deviations on and off 0=off 1=on\n", 0, # off for now - until we work on the interface for randomness
    "\n# var-cov matrix\n", file = filename, append = TRUE); write.t.sparse(Covq)

  
  # vmodel no random effects for variances
  filename <- paste0(wkdir,'/vmodel.cfg')

  cat("# variance model config for the a4a model",
    "\n# model params\n", ncol(Xv), 
    "\n# number of rows\n", nrow(Xv),
    "\n# design matrix\n", file = filename); write.t(Xv)


  # n1model no random effects for initial ages
  filename <- paste0(wkdir,'/ny1model.cfg')

  cat("# initial age structure model config for the a4a model",
    "\n# model params\n", ncol(Xny1), 
    "\n# number of rows\n", nrow(Xny1),
    "\n# design matrix\n", file = filename); write.t(Xny1)

   
  # rmodel
  filename <- paste0(wkdir,'/srrmodel.cfg')

  cat("# R model config for the a4a model",
    "\n# SR model ID:",srr $ srr,"\n", srr $ ID,
    "\n# SR CV:\n", srr $ srrCV,
    "\n# SPR0 :\n", srr $ SPR0,
    "\n# a model params\n", ncol(Xsra), 
    "\n# a model number of rows\n", nrow(Xsra),
    "\n# a model design matrix\n", file = filename); write.t(Xsra)
  cat("# b model params\n", ncol(Xsrb), 
    "\n# b model number of rows\n", nrow(Xsrb),
    "\n# b model design matrix\n", file = filename, append = TRUE); write.t(Xsrb)

  Covr <- getCov(nrow(Xsra), model = "iid", tau = 1)
  cat("# prior model for the SRR a param",
    "\n# flag to turn SRR a param deviations on and off 0=off 1=on\n", 0, # off for now - until we work on the interface for randomness
    "\n# var-cov matrix\n", file = filename, append = TRUE); write.t.sparse(Covr)
 
  # r internal model 
  filename <- paste0(wkdir,'/rmodel.cfg')

  cat("# internal model for recruits: orthoganal design",
    "\n# model params\n", ncol(Xr), 
    "\n# number of rows\n", nrow(Xr),
    "\n# design matrix\n", file = filename); write.t(Xr)
     
   # set up variable names
  pnames <- list(paste0("fMod:",colnames(Xf)), 
                 paste0("qMod:", unlist(sapply(1:length(indices), function(i) paste0(fleet.names[i+1],":",colnames(Xqlist[[i]]))))),
                 paste0("vMod:", unlist(sapply(1:length(fleet.names), function(i) paste0(fleet.names[i],":",colnames(Xvlist[[i]]))))),
                 paste0("n1Mod:",colnames(Xny1)), 
                 paste0("rMod:",colnames(Xr)), 
                 paste0("sraMod:",colnames(Xsra)), 
                 paste0("srbMod:",colnames(Xsrb)))    
  
  if (srr $ srrCV < 0) pnames <- pnames[-c(6,7)]   
     
  # end here if we just want to write the data and model files
  if (fit == "setup") return(wkdir)
 
 
  # ------------------------------------------------------------------------
  #
  # Run the ADMB executable (build up argument list first)
  #
  # ------------------------------------------------------------------------

  # run a4a split
  my.time.used[2] <- Sys.time() 
 
  # arguments
  args <- character(0)
  #TODO do we want to allow the use of this?
  if (fit == "MCMC") args <- c(args, paste0("-mcmc ", niters * 10," -mcsave 10"))
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
  } else if (os.type("windows")) {
    if (verbose) {
      echoc <- shell(paste0("cd /D", shQuote(wkdir), " & a4a", args))
    } else {
      echoc <- shell(paste0("cd /D", shQuote(wkdir), " & a4a ", args, " > logfile.txt"))
    }
  } else if (os.type("mac")) {
    stop("The FLa4a package is not developed for Macs yet...  Sorry!")
  }
  # hold off on the error for now - we can still get the hessian and perhaps do somehting with it...
  #if (echoc != 0) stop("Bad return from a4a executable: probably solution is not unique, try reducing parameters.\n")



  # ------------------------------------------------------------------------
  #
  # Read in ADMB output: how this is done will depend on the 'fit' argument
  #
  # ------------------------------------------------------------------------

  # post processing split
  my.time.used[3] <- Sys.time()


  if (fit == "MCMC") {
  
    filen <- file(paste0(wkdir, '/a4a.psv'), "rb")
    nopar <- readBin(filen, what = integer(), n = 1)
    out <- readBin(filen, what = numeric(), n = nopar * niters * 10)  ## TODO check this line
    close(filen)
    out <- matrix(out, byrow = TRUE, ncol = nopar)
    colnames(out) <- unlist(pnames)    

  } else if (fit %in% c("MP","assessment","Ext")) {

    # read admb output from file
    out <- list()

    # read .par file
    out[c("nopar","nlogl","maxgrad")] <- 
        as.numeric(scan(paste0(wkdir, '/a4a.par'), what = '', nlines = 1, quiet = TRUE)[c(6, 11, 16)])
    lin <- matrix(readLines(paste0(wkdir, '/a4a.par'))[-1], ncol = 2, byrow = TRUE)
    out $ par.est <- lapply(strsplit(sub(" ", "",lin[,2]), " "), as.numeric)
    names(out $ par.est) <- gsub("[# |:]", "", lin[,1])


    if (fit %in% c("assessment", "Ext")) {
    
      # check .cor file exisits
      if (!file.exists(paste0(wkdir, '/a4a.cor'))) stop("Hessian was not positive definate.")
      
      # read .cor file
      lin <- readLines(paste0(wkdir, '/a4a.cor'))
      npar <- length(lin) - 2
      out $ logDetHess <- as.numeric(strsplit(lin[1], '=')[[1]][2])
      
      sublin <- lapply(strsplit(lin[1:npar + 2], ' '), function(x) x[x != ''])
      par.names <- unlist(lapply(sublin, function(x) x[2]))
      par.std <- as.numeric(unlist(lapply(sublin, function(x) x[4])))

      out $ par.std <- lapply(names(out $ par.est), function(x) par.std[which(par.names==x)])
      names(out $ par.std) <- names(out $ par.est)
      
      #out $ cor <- matrix(NA, npar, npar, dimnames = list(par.names, par.names))
      #corvec <- unlist(sapply(1:length(sublin), function(i) sublin[[i]][5:(4+i)]))
      #out $ cor[upper.tri(out $ cor, diag = TRUE)] <- as.numeric(corvec)
      #out $ cor[lower.tri(out $ cor)] <- t(out $ cor)[lower.tri(out $ cor)]
      #out $ cov <- out $ cor * (par.std %o% par.std)

      # remove ssb and fbar final year from vcov and correlation matrices
      #out $ cov <- out $ cov[1:(npar-2), 1:(npar-2)]
      #out $ cor <- out $ cor[1:(npar-2), 1:(npar-2)]

      # use this as it seems to be more robust.. strangely...
      # I think it is because the solve method that ADMB uses is not
      # as good as the R one.... small numerical errors are
      # resulting in non-poitive definate vcov mats for subsets of parameters.
      hess <- getADMBHessian(wkdir) $ hes
      out $ cov <- solve(hess)
      out $ prec <- hess
      out $ nopar <- ncol(hess)
    }

    # read residuals
    #out $ res <- read.table(paste0(wkdir, '/a4a.res'), header = TRUE)

    # read derived model quantities
    ages <- sort(unique(full.df $ age))
    years <- sort(unique(full.df $ year))

    out $ N <- as.matrix(read.table(paste0(wkdir, '/n.out'), header = FALSE))
    out $ F <- as.matrix(read.table(paste0(wkdir, '/f.out'), header = FALSE))
    out $ Q <- as.matrix(read.table(paste0(wkdir, '/q.out'), header = FALSE))
    dim(out $ Q) <- c(length(years), length(indices), length(ages))
    out $ Q <- aperm(out $ Q, c(3,1,2))
  
    colnames(out $ N) <- colnames(out $ F) <- ages
    rownames(out $ N) <- rownames(out $ F) <- years
    dimnames(out $ Q) <- list(ages, years, names(indices))

  }

  # ------------------------------------------------------------------------
  #
  # convert into appropriate class
  #
  # ------------------------------------------------------------------------

  if (fit %in% c("MP", "assessment", "debug", "Ext")) {
  
    a4aout <- a4aFit()
 
    ind.names <- names(indices)

    stk.n <- t(out $ N)
    names(dimnames(stk.n)) <- c("age","year")

    hvst <- t(out $ F)
    names(dimnames(hvst)) <- c("age","year")
    
    logq <- lapply(1:length(indices), function(i) { 
                                        x <- drop(out $ Q[,,i])
                                        names(dimnames(x)) <- c("age","year")
                                        if (is.na(range(indices[[i]])["min"])) {
                                          x <- with(dimnames(index(indices[[i]])), x[1, year,drop = FALSE])
                                          dimnames(x)[[1]] <- "all"
                                        } else {
                                          x <- with(dimnames(index(indices[[i]])), x[age, year])
                                        }
                                        FLQuant(x)}) 
    names(logq) <- ind.names
 
    # First the a4aFit bits
    a4aout @ name    <- stock @ name
    a4aout @ desc    <- stock @ desc
    a4aout @ range   <- stock @ range
    a4aout @ call    <- match.call()
    a4aout @ stock.n <- FLQuant(stk.n) * exp(center.log[1])
    a4aout @ harvest <- FLQuant(hvst, units = "f")
    
    Z <- a4aout @ harvest + m(stock)  # TODO what if surveys are bigger... require that stock is bigger but with NA catches!
    
    a4aout @ catch.n <- a4aout @ harvest / Z * (1 - exp(-Z)) * a4aout @ stock.n

    index <- lapply(1:length(indices), function(i) {
                      dmns <- dimnames(logq[[i]]) 
                      if (dmns $ age[1] == "all") {
                        exp(logq[[i]] - center.log[1] + center.log[i+1]) * yearTotals(stock.n(a4aout)[, dmns[[2]]] * stock.wt(stock)[, dmns[[2]]])[1,]  
                      } else {
                        stock.n(a4aout)[dmns[[1]], dmns[[2]]] * exp(logq[[i]]  - center.log[1] + center.log[i+1])
                      }})
    names(index) <- ind.names

    a4aout @ index <- FLQuants(index)

    tmpSumm <- with(out, c(nopar, nlogl, maxgrad, nrow(df.data)))
                             #logDetHess = logDetHess,           
    a4aout @ fitSumm <- array(tmpSumm, dimnames = list(c("nopar","nlogl","maxgrad","nobs")))                                  
                              
    
    if (fit %in% c("assessment", "debug", "Ext")) {

      # fill up an a4aFitSA object
      a4aout <- a4aFitSA(a4aout)
    
      # now the a4aFitSA bits
      #tmpSumm <- with(out, c(nopar, nlogl, maxgrad, nrow(df.data)))
                             #logDetHess = logDetHess, 
                             
      #a4aout @ fitSumm <- array(tmpSumm, dimnames = list(c("nopar","nlogl","maxgrad","nobs")))                                  
                                   
      a4aout @ pars <- new("SCAPars")
  
      #
      # fill up stkmodel
      a4aout @ pars @ stkmodel @ name      <- a4aout @ name
      a4aout @ pars @ stkmodel @ desc      <- a4aout @ desc
      a4aout @ pars @ stkmodel @ range     <- a4aout @ range
      a4aout @ pars @ stkmodel @ centering <- center.log[1]
      a4aout @ pars @ stkmodel @  fMod     <- fmodel
      a4aout @ pars @ stkmodel @ n1Mod     <- n1model 
      a4aout @ pars @ stkmodel @ srMod     <- srmodel
      a4aout @ pars @ stkmodel @ m         <- m(stock)
      
      pars <- out $ par.est
                                        
      active <- sapply(out $ par.std, length) > 0

      dimnames(out $ cov) <- dimnames(out $ prec) <- list(unlist(pnames), unlist(pnames))
      stkactive <- active
      stkactive[2:3] <- FALSE
      a4aout @ pars @ stkmodel @ params <- FLPar(structure(unlist(pars[stkactive]), names = unlist(pnames[stkactive])))
      units(a4aout @ pars @ stkmodel @ params) <- "NA"

      a4aout @ pars @ stkmodel @ distr <- "norm"
      whichcol <-  grep(paste("(^",c("f","n1","r","sra","srb"),"Mod:)",collapse="|",sep=""), unlist(pnames))
      # we can use the inverse of a subset of the precision matrix
      # if we want the vcov matrix conditional on the
      # other (qmodel etc.) parameter estimates.
      ##a4aout @ pars @ stkmodel @ vcov <- solve(out $ prec[whichcol, whichcol])
      # or just the full vcov matrix, unconditional on the other things...
      a4aout @ pars @ stkmodel @ vcov <- out $ cov[whichcol, whichcol]



      #
      # fill up qmodel
      qmodels <- 
        lapply(seq_along(indices), 
          function(i)
          {
             which <- grepl(fleet.names[i+1], pnames[[2]])
             submodel(Mod = qmodel[[i]],
                      params = FLPar(structure(pars[[2]][which], names = pnames[[2]][which])),
                      vcov = out $ cov[pnames[[2]][which],pnames[[2]][which], drop = FALSE],
                      distr = "norm",
                      centering = center.log[i+1],
                      name = fleet.names[i+1],
                      desc = indices[[i]] @ desc,
                      range = indices[[i]] @ range)
          })
      
      names(qmodels) <- fleet.names[-1]
      a4aout @ pars @ qmodel <- submodels(qmodels)
      #names(a4aout @ pars @ qmodel) <- fleet.names[-1]
            
      #
      # fill up vmodel
      vmodels <- 
        lapply(seq_along(fleet.names), 
          function(i)
          {
             which <- grepl(fleet.names[i], pnames[[3]])
             submodel(Mod = vmodel[[i]],
                      params = FLPar(structure(pars[[3]][which], names = pnames[[3]][which])),
                      vcov = out $ cov[pnames[[3]][which],pnames[[3]][which], drop = FALSE],
                      distr = "norm",
                      centering = 0,
                      name = fleet.names[i],
                      desc = "",
                      range = if (i==1) stock @ range else indices[[i-1]] @ range)
          })
      
      names(vmodels) <- fleet.names
      a4aout @ pars @ vmodel <- submodels(vmodels)        
      #names(a4aout @ pars @ vmodel) <- fleet.names
 
 
      if (fit == "Ext") {
  
        # just while code is developing
        a4aout <- a4aFitExt(a4aout)
        a4aout @ Sigma <- out $ cov
        a4aout @ L <- chol(out $ cov)
        # consider adding L to submodels in a4aFitSA class...
        pars <- unlist(out $ par.est)
        names(pars) <- unlist(pnames)
        which <- names(pars) %in% colnames(out $ cov)
        a4aout @ baseLvlPars <- pars[which]
        a4aout @ designMatrix <- list(Xf = Xf, Xq = Xq, Xn1 = Xny1, Xv = Xv, Xsra = Xsra, Xsrb = Xsrb, Xr = Xr)

        # info about where to split?        
      }
    }
  } else if (fit == "MCMC") {
  
    # no MCMC class stuff yet... could use FLQuantDistr?
    a4aout <- a4aFitMCMC()
    a4aout @ pars <- out    

  }

  # end time
  my.time.used[4] <- Sys.time() 
  a4aout @ clock <- c("Pre-processing"  = diff(my.time.used)[1], 
                      "Running a4a"     = diff(my.time.used)[2], 
                      "Post-processing" = diff(my.time.used)[3], 
                       Total            = my.time.used[4] - my.time.used[1])
  
  units(a4aout @ harvest) <- "f"  

  
  # remove temporary directory - keep only true when dir is not temp dir
  if (!keep) unlink(wkdir)

  #
  #attr(a4aout, "stuff") <- out
  return(a4aout)  
}
  
  
#  out @ index.name <- ind.names
#  out @ logq <- do.call(FLQuants, logq)

#  out @ catch.lvar <- getres(fit $ res, "sd", 1)^2 
#  out @ catch.lhat <- getres(fit $ res, "pred", 1) + center.log[1]
#  out @ catch.lres <- getres(fit $ res, "res", 1) 

#  out @ index <- lapply(getres2(fit $ res, "obs"), exp)
#  out @ index.lvar <- lapply(getres2(fit $ res, "sd"), "^", 2) 
#  out @ index.lhat <- getres2(fit $ res, "pred")
#  out @ index.lres <- getres2(fit $ res, "res") 
  
#  for (i in 1:length(out @ index)) {
#    out @ index[[i]]     <- out @ index[[i]] * exp(center.log[i+1])
#    out @ index.lhat[[i]] <- out @ index.lhat[[i]] + center.log[i+1]
#    out @ logq[[i]] <- out @ logq[[i]] + center.log[i+1] - center.log[1]
#  }  
  
#  out @ logQ <- out @ index
#  for (i in 1:length(out @ logQ)) out @ logQ[[i]] @ .Data <- fit $ logQ[[i]]


#  out @ coefficients <- fit $ est
#  out @ covariance <- fit $ cov
#  names(out @ coefficients) <- fit $ names
#  colnames(out @ covariance) <- rownames(out @ covariance) <- fit $ names





