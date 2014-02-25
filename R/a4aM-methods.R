#' Method to compute natural mortality 
#'
#' @param object a \code{a4aM} object
#' @param grMod a \code{a4aGr} object to get the K from 
#' @param ... placeolder for covariates of the models. The names must match formula's variables (not parameters), with the exception of the \code{a4aGr} individual growth model. To use a growth model it must be called \code{grMod} and be of class \code{a4aGr}, in which case the parameters will be matched. The main objective if to be able to use \code{K} from von Bertalanffy models in M. 
#' @details If the models use \code{age} and/or \code{year} as terms the method expects these to be included in the call (will be passed through the \ldots argument). If they're not, the method will use the range slot to work out the ages and/or years that should be predicted. If \code{age} and/or \code{year} are not model terms, the method will use the range slot to define the dimensions of the resulting \code{M} \code{FLQuant}.    
#' @return a \code{FLQuant} object
#' @aliases m,a4aM-method
#' @examples
#' age <- 0:15
#' k <- 0.4
#' shp <- eval(as.list(~exp(-age-0.5))[[2]], envir=list(age=age))
#' lvl <- eval(as.list(~1.5*k)[[2]], envir=list(k=k))
#' M <- shp*lvl/mean(shp)
#' mod1 <- FLModelSim(model=~exp(-age-0.5))
#' mod2 <- FLModelSim(model=~1.5*k, params=FLPar(k=0.4))
#' m1 <- a4aM(shape=mod1, level=mod2)
#' rngmbar(m1)<-c(0,15)
#' m(m1, age=0:15)
#' all.equal(M, c(m(m1, age=0:15)))
#' # another example m
#' rngmbar(m1) <- c(2,4)
#' m(m1, age=2:10)
#' # with iters
#' mod2 <- FLModelSim(model=~k^0.66*t^0.57, params=FLPar(matrix(c(0.4,10,0.5,11), ncol=2, dimnames=list(params=c("k","t"), iter=1:2))), vcov=array(c(0.004, 0.00,0.00, 0.001), dim=c(2,2,2)))
#' m2 <- a4aM(shape=mod1, level=mod2)
#' m(m2, age=0:10)
#' m3 <- a4aM(shape=mod1, level=mvrnorm(100, mod2))
#' m(m3, age=0:15)
#' mod1 <- FLModelSim(model=~exp(-age-0.5))
#' mod3 <- FLModelSim(model=~1+b*v, params=FLPar(b=0.05))
#' mObj <- a4aM(shape=mod1, level=mvrnorm(100, mod2), trend=mod3, range=c(min=0,max=15,minyear=2000,maxyear=2003,minmbar=0,maxmbar=0))
#' m(mObj, v=1:4)

setMethod("m", "a4aM", function(object, grMod="missing", ...){
	args <- list(...)
	#if(length(args)==0) args$nocovariateprovided <- TRUE
	#browser()
	# check variables in models for "age" or "year"
	allVars <- c(all.vars(model(shape(object))), all.vars(model(level(object))), all.vars(model(trend(object))))
	if("age" %in% allVars) if(!("age" %in% names(args))) args$age <- vecage(object)
	if("year" %in% allVars) if(!("year" %in% names(args))) args$year <- vecyear(object)
	
	#if(sum(c("age", "year") %in% allVars)>0){
	#	if(!("age" %in% names(args))) args$age <- vecage(object)
	#	# else rngmbar(object) <- c(min(args$age), max(args$age)) 
	#	if(!("year" %in% names(args))) args$year <- vecyear(object)
	#}	

	# check if there is a growth model to get K from
	if(!missing(grMod)){
		k <- c("k","K")[c("k","K") %in% dimnames(params(level(object)))$params]
		params(level(object))[k] <- getK(grMod)
	}

	args$object <- shape(object)
	shp <- do.call("predict", args)
	args$object <- level(object)
	lvl <- do.call("predict", args)
	args$object <- trend(object)
	trd <- do.call("predict", args)
	
	# build the FLQuant for output
	mat <- matrix(1, ncol=6, nrow=3, dimnames=list(model=c("shp","lvl","trd"), dim=c("age","year","unit", "season","area","iter")))
	mat["shp",c(1,6)] <- dim(shp)
	mat["lvl",c(1,6)] <- dim(lvl)
	mat["trd",c(2,6)] <- dim(trd)

	dm <- apply(mat,2,max)
	nms12 <- list(age="all",year="1")
	
	# a bit spagethi ...
	if(!("age" %in% allVars)){
		v <- vecage(object)
		mat["shp",1] <- dm["age"] <- length(v)	
		nms12$age <- v
	} else {
		nms12$age <- args$age
	}

	if(!("year" %in% allVars)){
		v <- vecyear(object)
		dm["year"] <- length(v)	
		nms12$year <- v
	} else {
		nms12$year <- args$year
	}

	flq <- FLQuant(NA, dim=dm)
	dimnames(flq)[1:2] <- nms12

	flqs <- flq
	flqs[] <- shp

	flql <- aperm(flq, c(6,1,2,3,4,5))
	flql[] <- lvl
	flql <- aperm(flql,c(2,3,4,5,6,1))

	flqt <- aperm(flq, c(2,6,1,3,4,5))
	flqt[] <- trd
	flqt <- aperm(flqt, c(3,1,4,5,6,2))
	m <- flqs*flql*flqt/quantMeans(flqs[as.character(vecmbar(object))])[rep(1,dm[1])]
	FLQuant(m)
})

#' Method to simulate multivariate normal parameters
#'
#' @param n the number of simulations to be generated
#' @param mu a \code{a4aM} object
#' @return a \code{a4aM} object with n iterations 
#' @aliases mvrnorm,numeric,a4aM,missing,missing,missing,missing-method
#' @examples
#' mod1 <- FLModelSim(model=~exp(-age-0.5))
#' mod2 <- FLModelSim(model=~k^0.66*t^0.57, params=FLPar(matrix(c(0.4,10,0.5,11), ncol=2, dimnames=list(params=c("k","t"), iter=1:2))), vcov=array(c(0.004, 0.00,0.00, 0.001), dim=c(2,2,2)))
#' mod3 <- FLModelSim(model=~1+b*v, params=FLPar(b=0.05))
#' mObj <- a4aM(shape=mod1, level=mod2, trend=mod3, range=c(min=0,max=15,minyear=2000,maxyear=2003,minmbar=0,maxmbar=0))
#' mObj <- mvrnorm(100, mObj)
#' m(mObj, v=c(1,1,1,1))

setMethod("mvrnorm", c(n="numeric", mu="a4aM", Sigma="missing",
	tol="missing", empirical="missing", EISPACK="missing"), function(n=1, mu) {
	args <- list()
	args$n <- n
	args$mu <- shape(mu)
	if(!is.empty(vcov(args$mu))) params(shape(mu)) <- params(do.call("mvrnorm", args))	

	args$mu <- level(mu)
	if(!is.empty(vcov(args$mu))) params(level(mu)) <- params(do.call("mvrnorm", args))	

	args$mu <- trend(mu)
	if(!is.empty(vcov(args$mu))) params(trend(mu)) <- params(do.call("mvrnorm", args))	

	mu	
})

