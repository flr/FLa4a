
a4aInternal2 <- function(stock, indices, fmodel = defaultFmod(stock), qmodel = defaultQmod(indices), srmodel = defaultSRmod(stock), n1model = defaultN1mod(stock), vmodel = defaultVmod(stock, indices), covar=missing, wkdir=missing, verbose = FALSE, fit = "assessment", center = TRUE, mcmc=missing, new = FALSE)
{
  if (!new) {
    a4aInternal(
      stock,
      indices,
      fmodel = defaultFmod(stock),
      qmodel = defaultQmod(indices),
      srmodel = defaultSRmod(stock),
      n1model = defaultN1mod(stock),
      vmodel = defaultVmod(stock, indices),
      covar = missing,
      wkdir = missing,
      verbose = FALSE,
      fit = "assessment",
      center = TRUE,
      mcmc = missing
    )
  } else {
    stop("not yet")
    a4aInternal3()
  }
}

a4aInternal3 <- function(
    stock, indices, stk_model,
    wkdir = missing, verbose = FALSE, fit = "assessment",
    center = TRUE, mcmc = missing)
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

#  if (useADMB) { # fit using the ADMB code

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

#  } else { # fit using the TMB code
#    fitList <- fitTMB(fit, wkdir, df.data, stock, indices, full.df,
#                      fbar, plusgroup, surveytime, fleet.names,
#                      Xf, Xq, Xv, Xny1, srr, Xsra, Xsrb, Xr, Xvlist, Xqlist,
#                      my.time.used, mcmc, verbose)

#    if (fit == "setup") return(fitList)

#    out <- fitList$out
#    my.time.used <- fitList$my.time.used
#    wkdir <- fitList$wkdir
#    convergence <- fitList$convergence
#    pnames <- fitList$pnames
#    ages <- fitList$ages
#    years <- fitList$years
#  }

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
	a4aout@stock.n <- FLQuant(stk.n, units=units(catch.n(stock))) * exp(center.log[1])
	a4aout@harvest <- FLQuant(hvst, units = "f")
	Z <- a4aout@harvest + m(stock)
	a4aout@catch.n <- a4aout@harvest / Z * (1 - exp(-Z)) * a4aout@stock.n
  units(a4aout@catch.n) <- units(catch.n(stock))

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
	idna <- !is.na(catch.n(stock))
	cgcv <- length(a4aout@catch.n[idna, drop=TRUE]) * sum(c(log(catch.n(stock)/a4aout@catch.n))[idna, drop=TRUE]^2)/sum(1-flev)^2
	tmpSumm <- with(out, c(nopar, nlogl, maxgrad, nrow(df.data), cgcv, convergence, NA))

  # add in likelihood components here, using simple names for now.
  # comp1 is fleet 1 (catch), comp2 is fleet 2 (survey 1) etc.
  # if there is a SR relationship, it is the last comp
  nlogl_comps <- out$nlogl_comps
  nlogl_comps_names <- paste0("nlogl_comp", 1:length(nlogl_comps))

	a4aout@fitSumm <-
    array(c(tmpSumm, nlogl_comps),
          dimnames = list(c("nopar","nlogl","maxgrad","nobs","gcv", "convergence", "accrate", nlogl_comps_names)))

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
    a4aout@pars@stkmodel@units     <- units(catch.n(stock))
		pars <- out$par.est
		active <- sapply(out$par.std, length) > 0
		# add in iter dimension
		dim(out$cov) <- c(dim(out$cov), 1)
		dimnames(out$cov) <- list(unlist(pnames), unlist(pnames), 1)
		stkactive <- active
		if (convergence>0) stkactive[c(2,3,7)] <- FALSE else stkactive[c(2,3)] <- FALSE
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
				       #centering = FLPar(centering = center.log[i+1] - center.log[1]), # must also subtract catch scaling
				       centering = FLPar(centering = center.log[i+1]), # must not! subtract catch scaling
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

    #if (!useADMB) stop("not implemented in TMB version")

		# coerce into a4aFitMCMC
		a4aout <- a4aFitMCMC(a4aout, mcmc=mcmc)
		# fill parameters
		a4aout@pars@stkmodel@coefficients <- propagate(a4aout@pars@stkmodel@coefficients, out$mcmc$nit)
		a4aout@pars@stkmodel@coefficients[] <- t(out$mcmc$mcmcout[,dimnames(a4aout@pars@stkmodel@coefficients)[[1]]])
		dimnames(a4aout@pars@stkmodel@coefficients)[[2]] <- 1:out$mcmc$nit
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
			a4aout@pars@vmodel[[i]]@coefficients <- propagate(a4aout@pars@vmodel[[i]]@coefficients, out$mcmc$nit)
			a4aout@pars@vmodel[[i]]@coefficients[] <- t(out$mcmc$mcmcout[,dimnames(a4aout@pars@vmodel[[i]]@coefficients)[[1]]])
			a4aout@pars@vmodel[[i]]@vcov[] <- NA

			if(i!="catch"){
				a4aout@pars@qmodel[[i]]@coefficients <- propagate(a4aout@pars@qmodel[[i]]@coefficients, out$mcmc$nit)
				a4aout@pars@qmodel[[i]]@coefficients[] <- t(out$mcmc$mcmcout[,dimnames(a4aout@pars@qmodel[[i]]@coefficients)[[1]]])
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
