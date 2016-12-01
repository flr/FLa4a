#' @title Get K
#' @name getK
#' @rdname getK
#' @description Method to get values of the growth parameter K
#' @param object an \code{a4aGr} object
#' @template dots
#' @return a \code{vector} with K values
#' @aliases getK getK-methods
#' @examples
#' mm <- matrix(NA, ncol=3, nrow=3)
#' diag(mm) <- c(50, 0.001,0.001)
#' mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
#' vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","yr^-1","yr")), vcov=mm, distr="norm")
#' vbObj <- mvrnorm(100,vbObj)
#' getK(vbObj)
setGeneric("getK", function(object, ...) standardGeneric("getK"))
#' @rdname getK
setMethod("getK", "a4aGr", function(object){
	k <- c("k","K")[c("k","K") %in% dimnames(params(object))$params]
	if(length(k)==1) params(object)[k]
}) 

#' @title mvrnorm 
#' @rdname mvrnorm-a4aGr
#' @name mvrnorm for a4aGR
#' @aliases mvrnorm,numeric,a4aGr-method
#' @description Method to generate multivariate normal parameters for \code{a4aGr} objects.
#' @param n the number of simulations to be generated
#' @param mu an \code{a4aGr} object
#' @return an \code{a4aGr} object with n iterations
#' @examples
#' mm <- matrix(NA, ncol=3, nrow=3)
#' diag(mm) <- c(50, 0.001,0.001)
#' mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
#' vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","yr^-1","yr")), vcov=mm, distr="norm")
#' vbObj <- mvrnorm(100,vbObj)
setMethod("mvrnorm", signature("numeric", "a4aGr"), function(n=1, mu) {
	args <- list()
	args$mu <- FLModelSim(model=grInvMod(mu), params=params(mu), vcov=vcov(mu), distr=distr(mu))
	args$n <- n
	res <- do.call("mvrnorm", args)	
	params(mu) <- params(res)
	mu	
})

#' @title mvrtriangle 
#' @name mvrtriangle for a4aGr 
#' @rdname mvrtriangle-a4aGr
#' @aliases mvrtriangle,numeric,a4aGr-method
#' @description Method to generate multivariate parameters with elliptical copulas and triangular marginals for \code{a4aGr} objects.
#' @param n the number of iterations
#' @template bothargs
#' @details The method is essentially a special case of \code{mvrcop}, where the copula is of type "ellipCopula" and family "t", and where the marginals are triangular.
#' @return an \code{a4aGr} object with n iterations
#' @examples
#' # Set up the a4aGr object and parameters for the marginals
#' mm <- matrix(NA, ncol=3, nrow=3)
#' diag(mm) <- c(50, 0.001,0.001)
#' mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
#' vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","yr^-1","yr")), vcov=mm, distr="norm")
#' pars <- list(list(a=50, b=100, c=58.5), list(a=0.06, b=0.2, c=0.086), list(a=0, b=0.005, c=0.001))
#'
#' # Apply mvrtriangle...
#'   vbObj1 <- mvrtriangle(10000, vbObj, paramMargins=pars)
#'   splom(data.frame(t(params(vbObj1)@@.Data)), pch=".")
#' #...and compare with mvrcop
#'   vbObj2 <- mvrcop(10000, vbObj, copula="ellipCopula", family="t", param=0, margins="triangle", paramMargins=pars)
#'   splom(data.frame(t(params(vbObj2)@@.Data)), pch=".")
setMethod("mvrtriangle", signature("numeric", "a4aGr"), function(n=1, object, ...) {
	args <- list(...)
	args$object <- FLModelSim(model=grInvMod(object), params=params(object), vcov=vcov(object), distr=distr(object))
	args$n <- n
	res <- do.call("mvrtriangle", args)	
	params(object) <- params(res)
	object	
})

#' @title mvrcop 
#' @name mvrcop for a4aGr
#' @rdname mvrcop-a4aGr
#' @aliases mvrcop,numeric,a4aGr-method
#' @description Method to generate multivariate parameters with user-defined copulas and marginals for \code{a4aGr} objects.
#' @param n the number of iterations
#' @param mvdc the \code{a4aGr} object
#' @param ... arguments to be passed to the rMvdc and copula methods
#' @return an \code{FLModelSim} object with n groups of parameters
#' @examples
#' mm <- matrix(NA, ncol=3, nrow=3)
#' diag(mm) <- c(50, 0.001,0.001)
#' mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.01,0.00004)
#' vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","yr^-1","yr")), vcov=mm, distr="norm")
#' pars <- list(list(a=50, b=100, c=58.5), list(a=0.06, b=0.2, c=0.086), list(a=0, b=0.005, c=0.001))
#' #In the following, the third, fourth and fifth arguments refer to the copula,
#' #  while the final two arguments refer to the marginal distributions:
#' vbObj <- mvrcop(10000, vbObj, copula="archmCopula", family="clayton", param=2, margins="triangle", paramMargins=pars)
#' splom(data.frame(t(params(vbObj)@@.Data)), pch=".")
setMethod("mvrcop", signature("numeric", "a4aGr"), function(n=1, mvdc, ...) {
	object <- mvdc
	rm(mvdc)
	args <- list(...)
	args$mvdc <- FLModelSim(model=grInvMod(object), params=params(object), vcov=vcov(object), distr=distr(object))
	args$n <- n
	res <- do.call("mvrcop", args)	
	params(object) <- params(res)
	object	
})

##' @title predict for \code{a4aGr} 
#' @name predict for a4aGr
#' @rdname predict-a4aGr
#' @description Predicts ages or lengths using a growth class
#' @param object the \code{a4aGr} object
#' @param ... arguments to be passed to the rMvdc and copula methods
#' @return a \code{matrix} object with lengths or ages
#' @examples
#' vbObj <- a4aGr(grMod=~linf*(1-exp(-k*(t-t0))), grInvMod=~t0-1/k*log(1-len/linf), params=FLPar(linf=58.5, k=0.086, t0=0.001, units=c("cm","yr^-1","yr")))
#' predict(vbObj, len=1:50+0.5)
#' predict(vbObj, t=1:20+0.5)
setMethod("predict", "a4aGr", function(object, ...){
	model <- object
	args <- list(...)
	ll <- all.vars(grInvMod(object))[!all.vars(grMod(object)) %in% all.vars(grInvMod(object))]
	if(ll %in% names(args)){
		args$object <- FLModelSim(model=grInvMod(model), params=params(model), vcov=vcov(model), distr=distr(model))
	} else {
		args$object <- FLModelSim(model=grMod(model), params=params(model), vcov=vcov(model), distr=distr(model))
	}
	do.call("predict", args)
})



