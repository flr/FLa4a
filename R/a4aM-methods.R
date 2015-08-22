#' @title natural mortality 
#' @description Method to compute natural mortality 
#' @name M 
#' @rdname M 
#' @param object a \code{a4aM} object
#' @param grMod a \code{a4aGr} object to get the K from 
#' @param ... placeholder for covariates of the models. The names must match formula's variables (not parameters), with the exception of the \code{a4aGr} individual growth model. To use a growth model it must be called \code{grMod} and be of class \code{a4aGr}, in which case the parameters will be matched. The main objective if to be able to use \code{K} from von Bertalanffy models in M. 
#' @details The method uses the range slot to define the quant and year dimensions of the resulting \code{M} \code{FLQuant}. The name fo the quant dimension is taken as the name of a variable that is present in the \code{shape} formula but is not present in the shape parameters (the \code{params} slot of the \code{shape} model). If more than one such than variable exists then there is a problem with the \code{shape} model definition.
#' @return a \code{FLQuant} object
#' @aliases m,a4aM-method
#' @template Example-a4aM
setMethod("m", "a4aM", function(object, grMod="missing", ...){
	args <- list(...)

    # Pick out the quant name from the shape model - the only one that goes over quant dimension
    # The range of this is set by vecquant(object)
    quantVars <- all.vars(model(shape(object)))
    quantParams <- dimnames(params(shape(object)))$params
    qname <- quantVars[!(quantVars %in% quantParams)]
    if(length(qname) < 1){
        qname <- "quant"
    }
    if(length(qname) > 1){
        stop("More than one variable in the shape formula with no params value.")
    }

	allVars <- c(all.vars(model(shape(object))), all.vars(model(level(object))), all.vars(model(trend(object))))

	args[[qname]] <- vecquant(object)
	args$year <- vecyear(object)

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

    # Use predict, formula argument is called 'object'
	args$object <- shape(object)
	shp <- do.call("predict", args)
	args$object <- level(object)
	lvl <- do.call("predict", args)
	args$object <- trend(object)
	trd <- do.call("predict", args)
	
	# build the FLQuant for output
	mat <- matrix(1, ncol=6, nrow=3, dimnames=list(model=c("shp","lvl","trd"), dim=c(qname,"year","unit", "season","area","iter")))
	mat["shp",c(1,6)] <- dim(shp)
	mat["lvl",c(1,6)] <- dim(lvl)
	mat["trd",c(2,6)] <- dim(trd)

	dm <- apply(mat,2,max)
	nms12 <- list(qname="all",year="1")
    names(nms12)[1] <- qname
	
	# a bit spagethi ...
    # if qname is not in shape model, the dimensions of the predicted shape values are not determined by range and we need to blow up dims to match rng
	if(!(qname %in% allVars)){  
		v <- vecquant(object)
		mat["shp",1] <- dm[qname] <- length(v)	
		nms12[[qname]] <- v
	} else {
	   nms12[[qname]] <- args[[qname]]
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
    names(dimnames(flq))[1] <- names(nms12)[1]

	flqs <- flq
	flqs[] <- shp

	flql <- aperm(flq, c(6,1,2,3,4,5))
	flql[] <- lvl
	flql <- aperm(flql,c(2,3,4,5,6,1))

	flqt <- aperm(flq, c(2,6,1,3,4,5))
	flqt[] <- trd
	flqt <- aperm(flqt, c(3,1,4,5,6,2))
	m <- flqs*flql*flqt/quantMeans(flqs[as.character(vecmbar(object))])[rep(1,dm[1])]
	return(FLQuant(m))
})

#' @title natural mortality 
#' @description Method to simulate multivariate normal parameters for a \code{a4aM} object. 
#' @name mvnorm 
#' @rdname mvnorm-a4aM 
#' @param n the number of simulations to be generated
#' @param mu a \code{a4aM} object
#' @return a \code{a4aM} object with n iterations 
#' @aliases mvrnorm,numeric,a4aM,missing,missing,missing,missing-method
#' @examples
#' mod1 <- FLModelSim(model=~exp(-age-0.5))
#' mod2 <- FLModelSim(model=~k^0.66*t^0.57, 
#'                    params=FLPar(matrix(c(0.4,10,0.5,11), ncol=2, 
#'                    	           dimnames=list(params=c("k","t"), iter=1:2))), 
#'                    vcov=array(c(0.004, 0.00,0.00, 0.001), dim=c(2,2,2)))
#' mod3 <- FLModelSim(model=~1+b*v, params=FLPar(b=0.05))
#' mObj <- a4aM(shape=mod1, level=mod2, trend=mod3, 
#'              range=c(min=0,max=15,minyear=2000,maxyear=2003,minmbar=0,maxmbar=0))
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

