#' @title Simulation with a copula model and triangular distributions
#' @description Simulates model parameters using elliptical copulas and triangular marginals.
#' @param n the number of iterations
#' @param object the \code{FLModelSim} object
#' @param ... arguments to be passed to the rMvdc and copula methods
#' @return an \code{FLModelSim} object with n sets of parameters
#' @rdname mvrtriangle
#' @examples
#' # Set up the FLModelSim object
#' mm <- matrix(NA, ncol=3, nrow=3)
#' diag(mm) <- c(100, 0.001,0.001)
#' mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.1,0.0003)
#' md <- ~linf*(1-exp(-k*(t-t0)))
#' prs <- FLPar(linf=120, k=0.3, t0=0.1, units=c("cm","yr^-1","yr"))
#' vb <- FLModelSim(model=md, params=prs, vcov=mm, distr="norm")
#'
#' # Simulate from a multivariate normal distribution...
#'   set.seed(1)
#'   vbSim <- mvrnorm(10000, vb)
#'   mm <- predict(vbSim, t=0:20+0.5)
#' #...from a multivariate triangular distribution with default ranges (0.01 and
#' #   0.99 quantiles for min and max using a normal distribution with mean from
#' #   params and sigma from vcov, and with the apex located at params)...
#'   set.seed(1)
#'   vbSim1 <- mvrtriangle(10000, vb)
#'   mm1 <- predict(vbSim1, t=0:20+0.5)
#' #...and from a multivariate triangular distribution with specified ranges 
#' #   (note if "c" is missing, it will take the average of "a" and "b")
#'   set.seed(1)
#'   pars <- list(list(a=90, b=125, c=120), list(a=0.2, b=0.4), list(a=0, b=0.4, c=0.1))
#'   vbSim2 <- mvrtriangle(10000, vb, paramMargins=pars)
#'   mm2 <- predict(vbSim2, t=0:20+0.5)
#'
#' # Plot the results
#' par(mfrow=c(3,1))
#' boxplot(t(mm), main="normal")
#' boxplot(t(mm1), main="triangular")
#' boxplot(t(mm2), main="triangular2")
#' splom(data.frame(t(params(vbSim)@@.Data)), pch=".")
#' splom(data.frame(t(params(vbSim1)@@.Data)), pch=".")
#' splom(data.frame(t(params(vbSim2)@@.Data)), pch=".")


setGeneric("mvrtriangle", function(n, object, ...) standardGeneric("mvrtriangle"))
#' @rdname mvrtriangle
setMethod("mvrtriangle", signature("numeric", "FLModelSim"), function(n=1, object, ...) {
		args <- list(...)	
		model <- object
		pars <- params(model)
		dm <- dim(pars)
		dnm <- dimnames(pars)
		# check that params second dim is "iter"
		if(names(dnm[2])!="iter") stop("To apply this method params must have 2 dimensions only and the second has to be \"iter\".")	
	
		#--------------------------------------------------
		# if we have "n" => flatten iterations
		#--------------------------------------------------
		mu <- iterMedians(pars)
		sig2 = apply(vcov(model),c(1,2),median)

		#--------------------------------------------------
		# set the copula
		#--------------------------------------------------
		rho <- cov2cor(sig2)
		args$dim <- dm[1]
		# if dispstr not set use "un"
		if(!("dispstr" %in% names(args))){
			args$dispstr <- "un"
			args$param <- rho[upper.tri(rho)]
		}
		# if family not set use "t" copula
		if(!("family" %in% names(args))) args$family <- "t"
		# go
		lst <- args[names(args) %in% names(formals(ellipCopula))]
		args$copula <- do.call("ellipCopula", lst)

		#--------------------------------------------------
		# set and run mvdc
		#--------------------------------------------------
		args$margins <- rep("triangle", dm[1])
		if(!("paramMargins" %in% names(args))){
			pmrg <- data.frame(mean=c(mu), sd=sqrt(diag(sig2)))
			pmrg <- t(apply(pmrg, 1, function(x) qnorm(c(0.01,0.99,0.5), x["mean"], x["sd"])))
			dimnames(pmrg) <- list(dimnames(mu)[[1]], c("a", "b", "c"))
			pmrg <- split(data.frame(pmrg), 1:nrow(pmrg))
			args$paramMargins <- lapply(pmrg, as.list)
		}

		lst <- args[names(args) %in% names(formals(mvdc))]
		mvobj <- do.call("mvdc", lst)	
		res <- rMvdc(n, mvobj)		

		#--------------------------------------------------
		# output
		#--------------------------------------------------
		if(n>1) res <- t(res) else res <- matrix(res, ncol=1)
		dnm$iter <- 1:n
		dimnames(res) <- dnm
		res <- FLPar(res)
		units(res) <- units(mu)
		return(FLModelSim(model=model(model), params=res, vcov=vcov(model), distr=unique(with(args, paste(dispstr, copula@fullname, margins)))))
	}
)

#' @title Simulation using copula models
#' @description Simulates model parameters with user-defined copulas and marginals.
#' @param n the number of iterations
#' @param mvdc an \code{FLModelSim} object
#' @param copula the name of the copula to be used
#' @param ... arguments to be passed to the copula methods
#' @return an \code{FLModelSim} object with n groups of parameters
#' @rdname mvrcop
#' @examples
#' mm <- matrix(NA, ncol=3, nrow=3)
#' diag(mm) <- c(100, 0.001,0.001)
#' mm[upper.tri(mm)] <- mm[lower.tri(mm)] <- c(0.1,0.1,0.0003)
#' md <- ~linf*(1-exp(-k*(t-t0)))
#' prs <- FLPar(linf=120, k=0.3, t0=0.1, units=c("cm","yr^-1","yr"))
#' vb <- FLModelSim(model=md, params=prs, vcov=mm, distr="norm")
#' pars <- list(list(a=90, b=125, c=120), list(a=0.2, b=0.4), list(a=0, b=0.4, c=0.1))
#' vbSim <- mvrcop(10000, vb, copula="archmCopula", family="clayton", param=2, 
#'    margins="triangle", paramMargins=pars)
#' boxplot(t(predict(vbSim, t=0:20+0.5)))
#' splom(data.frame(t(params(vbSim)@@.Data)), pch=".")
setGeneric("mvrcop", function(n, mvdc, ...) standardGeneric("mvrcop"))
#' @rdname mvrcop
setMethod("mvrcop", signature("numeric", "FLModelSim"), function(n, mvdc, copula, ...) {

		args <- list(...)	
		model <- mvdc
		rm(mvdc)
		pars <- params(model)
		dm <- dim(pars)
		dnm <- dimnames(pars)
		# check that params second dim is "iter"
		if(names(dnm[2])!="iter") stop("To apply this method params must have 2 dimensions only and the second has to be \"iter\".")	
	
		#--------------------------------------------------
		# if we have "n" => flatten iterations
		#--------------------------------------------------
		mu <- iterMedians(pars)
		args$dim <- dm[1]

		#--------------------------------------------------
		# set the copula
		#--------------------------------------------------
		lst <- args[names(args) %in% names(formals(copula))]
		args$copula <- do.call(copula, lst)

		#--------------------------------------------------
		# set and run mvdc
		#--------------------------------------------------
		lst <- args[names(args) %in% names(formals(mvdc))]
		if(length(lst$margins)==1) lst$margins <- rep(lst$margins, args$dim)
		mvobj <- do.call("mvdc", lst)	
		res <- rMvdc(n, mvobj)		

		#--------------------------------------------------
		# output
		#--------------------------------------------------
		if(n>1) res <- t(res) else res <- matrix(res, ncol=1)
		dnm$iter <- 1:n
		dimnames(res) <- dnm
		res <- FLPar(res)
		units(res) <- units(mu)
		return(FLModelSim(model=model(model), params=res, vcov=vcov(model), distr=unique(with(args, paste(copula@fullname, margins)))))
	
	}
)

#' @title Check that the second dimension in params is "iter"
#' @name pars2dim
#' @rdname pars2dim-methods
#' @template object
#' @aliases pars2dim pars2dim-methods
#' @description Checks that the name of the second dimension in params is "iter". For internal use, not very interesting for users. It takes a \code{FLModelSim} object and returns a \code{logical}.
#' @examples
#' pars2dim(FLModelSim())
setGeneric("pars2dim", function(object) standardGeneric("pars2dim"))
#' @rdname pars2dim-methods
#' @aliases pars2dim,FLModelSim-method
setMethod("pars2dim", "FLModelSim", function(object) {

	pars <- params(object)
	dnm <- dimnames(pars)
	names(dnm)[2]=="iter"


})

