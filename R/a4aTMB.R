# a4aTMB.R - DESC
# /a4aTMB.R

# Copyright European Union, 2016
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@jrc.ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

a4aDM <- function(stock, indices, fmodel  = ~ s(age, k = 3) + factor(year), 
                qmodel  = lapply(seq(length(indices)), function(i) ~ 1), 
                srmodel = ~ factor(year),
                n1model = ~ factor(age), 
                vmodel  = lapply(seq(length(indices) + 1), function(i) ~ 1),
                covar=missing, verbose = FALSE, fit = "assessment", 
                center = TRUE, mcmc=missing)
{

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

  #
  dms <- dimnames(stock.n(stock))
  grid <- do.call(expand.grid, c(dimnames(catch.n(stock))[c(3,5)], list(iter = 1:max(dms$iter))))
  
  niters <- nrow(grid)
  res <- vector('list', length=niters)

  # LOOP over niters
  for (its in seq(niters)) {

    indices <- lapply(indices, function(x){
      idx <- x[,, grid$unit[its],, grid$area[its]]
		  iter(idx, min(grid$iter[its], dims(x)$iter))
    })
    
    indices <- FLIndices(indices)

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
  list.obs <- c(list(catch = FLa4a:::quant2mat(catch.n(stock)@.Data)), lapply(indices, function(x) FLa4a:::quant2mat(index(x)@.Data)))
	                   
	# convert the variances of catches and indices to a list of named arrays
	list.var <- c(list(catch = FLa4a:::quant2mat(catch.n(stock)@var)), lapply(indices, function(x) FLa4a:::quant2mat(index.var(x))))                
	                   
	# calculate appropriate centering for observations on log scale
	# a bit spaguetti ... if center is a numeric vector only those elements will be centered
	center.log <- sapply(list.obs, function(x) mean(log(x), na.rm = TRUE))
	if(is.numeric(center)) center.log[-center][] <- 0 else if(!isTRUE(center)) center.log[] <- 0 

	# convert to dataframe. NOTE: list2df also logs the observations and centers
	df.data <- do.call(rbind, lapply(1:length(list.obs), FLa4a:::list2df, list.obs=list.obs, list.var=list.var, center.log=center.log))
	
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
	full.df <- do.call(rbind, lapply(1:length(list.obs), function(i) cbind(fleet = i, FLa4a:::make.df(i, stock=stock, indices=indices))))

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
	Xqlist <- lapply(seq_along(indices), function(i) getX(qmodel[[i]], df=subset(full.df, fleet == fleet.names[i+1])))
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
	# NOTE: move to internal funs ?? map with FLSR !? 
	#-------------------------------------------------------------------------
	# process recruitment formula: 
	# builder functions - could be more hadley...
	bevholt <- function(a = ~ 1, b = ~ 1, CV = 0.5) {
		if (CV <= 0) stop ("CV in stock recruit relationship cannot be less than zero")
		list(srr = "bevholt", a = a, b = b, SPR0 = 1, srrCV = CV, ID = 1)
	}
	bevholtSV <- function(h = ~ 1, v = ~ 1, SPR0 = 1, CV = 0.5) {
		if (CV <= 0) stop ("CV in stock recruit relationship cannot be less than zero")
		list(srr = "bevholtSV", a = h, b = v, SPR0 = SPR0, srrCV = CV, ID = 5)
	}
	ricker <- function(a = ~ 1, b = ~ 1, CV = 0.5) {
		if (CV <= 0) stop ("CV in stock recruit relationship cannot be less than zero")
		list(srr = "ricker", a = a, b = b, SPR0 = 1, srrCV = CV, ID = 2)
	}
	hockey <- function(a = ~ 1, b = ~ 1, CV = 0.5) {
		if (CV <= 0) stop ("CV in stock recruit relationship cannot be less than zero")
		list(srr = "hockey", a = a, b = b, SPR0 = 1, srrCV = CV, ID = 3)
	}
	geomean <- function(a = ~ 1, CV = 0.5) {
		if (CV <= 0) stop ("CV in stock recruit relationship cannot be less than zero")
		list(srr = "geomean", a = a, b = ~ 1, SPR0 = 1, srrCV = CV, ID = 4)
	}
	none <- function() list(srr = "geomean", a = ~ 1, b = ~ 1, srrCV = -1, ID = 4)
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	#-------------------------------------------------------------------------
	# NOTE: how are covars being passed ? 
	#-------------------------------------------------------------------------
	# now separate model elements, find SR models, stop if more than one specified
	facs <- strsplit(as.character(srmodel)[length(srmodel)], "[+]")[[1]]
	facs <- gsub("(^ )|( $)", "", facs) # remove leading and trailing spaces
	a4as <- grepl(paste("(^",c("bevholt", "bevholtSV", "ricker","hockey","geomean"),"[(])", collapse = "|", sep = ""), facs)
	if (sum(a4as) > 1) stop("you can only specify one type of stock recruit relationship.")
	srrmod <- if (sum(a4as) == 0) "none()" else facs[a4as]
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

	# internal r model matrix (this if is doing nothing)
	if (sum(a4as) == 0) rmodel <- srmodel else rmodel <- ~ factor(year) 
	Xr <- getX(rmodel, subset(full.df, age == min(age) & fleet == "catch"))
	#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	#========================================================================
	# Write model matrices
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

	#-------------------------------------------------------------------------
	# write config list
	#-------------------------------------------------------------------------
  
	df.aux <- unique(full.df[c("year","age","m","m.spwn","harvest.spwn","mat.wt","stock.wt")])

Ldat<-list(
ageRange=range(full.df$age),
yearRange=range(full.df$year),
surveyMinAge=srvMinAge,
surveyMaxAge=srvMaxAge,
surveyTimes=surveytime,
fbarRange=fbar,
obs=as.matrix(df.data),
aux=as.matrix(df.aux),
designF=Xf,
designQ=Xq,
designV=Xv,
designNy1=Xny1,
designR=Xr,
designRa=Xsra,
designRb=Xsrb,
srCV=srr$srrCV,
spr0=ifelse(!is.null(srr$SPR0), srr$SPR0, 0),
Rmodel=srr$ID,
isPlusGrp=plusgroup
)

Ldat$locFleetVec=Ldat$obs[,1]
Ldat$locYearVec=Ldat$obs[,2]
Ldat$locAgeVec=Ldat$obs[,3]

Lpin=list(
	fpar=rep(0, ncol(Ldat$designF)),
	qpar=rep(0, ncol(Ldat$designQ)),
	vpar=rep(0, ncol(Ldat$designV)),
	ny1par=rep(0, ncol(Ldat$designNy1)),
	rpar=rep(0, ncol(Ldat$designR)),
	rapar=rep(0, ncol(Ldat$designRa)),
	rbpar=rep(0, ncol(Ldat$designRb))
	)

res[[its]] <- list(Ldat=Ldat, Lpin=Lpin)

  }
  names(res) <- paste0('it', seq(niters))

  return(res)
}
