###############################################################################
# EJ&CM(2012821)
# Methods to add uncertainty to FLQuant objects
# NOTE #1:
# NOTE #2:
###############################################################################

#' Methods to generate FLStock objects
#' @description This method computes the \code{FLStock} slots consistently with the information provided by the \code{FLQuant}. It requires two of the triplet R/C/F to compute the third consistent with Baranov and survival's equations.
#' @param object an FLStock
#' @param R an FLQuant with iterations or missing
#' @param C an FLQuant with iterations or missing
#' @param F an FLQuant with iterations or missing
#' @param ... additional argument list that might not ever be used.
#' @return an FLStock
#' @docType methods
#' @rdname genFLStock-methods
#' @aliases genFLStock genFLStock-methods

setGeneric("genFLStock", function(object, R, C, F, ...) standardGeneric("genFLStock"))

#' @rdname genFLStock-methods
setMethod("genFLStock", c("FLStock", "FLQuant", "FLQuant", "missing"), function(object, R, C, F, ...){
	cat("Not implemented yet\n")
})

#' @rdname genFLStock-methods
setMethod("genFLStock", c("FLStock", "missing", "FLQuant", "FLQuant"), function(object, R, C, F, ...){
	cat("Not implemented yet\n")
})

#' @rdname genFLStock-methods
setMethod("genFLStock", c("FLStock", "FLQuant", "missing", "FLQuant"), function(object, R, C, F, ...){
	# requires checking dimensions
	if(!identical(dim(catch.n(object))[-c(1,6)], dim(R)[-c(1,6)])) stop("Recruitment vector must have consistent dimensions with the stock object")
	if(!identical(dim(catch.n(object))[-6]    , dim(F)[-6])) stop("Harvest matrix must have consistent dimensions with the stock object")
	if(dim(R)[6]!=dim(R)[6]) stop("R and F must have the same number of iterations")

	# get dims and set flq
	dms <- dims(object)
	nages <- dms$age
	nyrs <- dms$year
	niters <- dim(R)[6]
	flq <- FLQuant(dimnames=dimnames(F))
	
	# compute cumulative Z
	Z <- FLCohort(F + m(object))
	Z[] <- apply(Z, c(2:6), cumsum)

	# expand variability into [N] by R*[F] 
	Ns <- FLCohort(R[rep(1,nages)])
	Ns <- Ns*exp(-Z)
	Ns <- as(Ns, "FLQuant")

	# Update object
	stock.n(object) <- flq
	# [R]
	stock.n(object)[1] <- R
	# [N]
	stock.n(object)[-1,-1] <- Ns[-nages,-nyrs] 
	# plus group
	stock.n(object)[nages,-1] <- Ns[nages,-1] + stock.n(object)[nages,-1]
	stock(object) <- computeStock(object)
	# [F]
	harvest(object) <- F
	# [C]
	Z <- harvest(object) + m(object)
	Cs <- harvest(object)/Z*(1-exp(-Z))*stock.n(object) 
	catch.n(object) <- Cs
	catch(object) <- computeCatch(object)
	# [L] & [D] rebuilt from C
	# Ds=D/(D+L)*Cs where Cs is the simulated catch
	D <- discards.n(object)
	L <- landings.n(object)
	discards.n(object) <- D/(D+L)*Cs
	discards(object) <- computeDiscards(object)
	landings.n(object) <- Cs - discards.n(object)
	landings(object) <- computeLandings(object)
	# out
	object
})

#' @name getAcor
#' @rdname getAcor-methods
#' @title compute log-correlation matrix
#' @description Method to compute the log-correlation matrix for the first dimension (\code{quant}) of the \code{FLQuant} object.
#' @template bothargs
#' @return an \code{FLQuant} object with a quant log-correlation matrix
#' @aliases getAcor getAcor-methods
#' @examples
#' data(ple4)
#' getAcor(harvest(ple4))

setGeneric("getAcor", function(object, ...) standardGeneric("getAcor"))
#' @rdname getAcor-methods
setMethod("getAcor", c("FLQuant"), function(object, ...) {
		mu <- log(object)
		Rho <- cor(t(mu[drop = TRUE]))
		return(Rho)
})

#' @title Methods to generate FLQuant objects
#' @description This method uses the quant log-correlation matrix of the \code{FLQuant} object and generates a new \code{FLQuant} using a lognormal multivariate distribution.
#' @param object an FLQuant
#' @param cv the coefficient of variation
#' @param method the method used to compute the correlation matrix; for now only "ac" (autocorrelation) is implemented
#' @param niter the number of iterations to be generated
#' @param ... additional argument list that might not ever be used.
#' @return an FLQuant
#' @docType methods
#' @rdname genFLQuant-methods
#' @aliases genFLQuant genFLQuant-methods
#' @examples
#' data(ple4)
#' sim.F <- genFLQuant(harvest(ple4))
setGeneric("genFLQuant", function(object, ...) standardGeneric("genFLQuant"))

#' @rdname genFLQuant-methods
setMethod("genFLQuant", "FLQuant", 
	function(object, cv = 0.2, method = "ac", niter = 250) {
  # use log transform, to be expanded on later versions
	mu <- log(object)
	if(method == "ac") {
		Rho <- cor(t(mu[drop = TRUE]))
		flq <- mvrnorm(niter * dim(mu)[2], rep(0, nrow(Rho)), log(cv^2+1) * Rho)
		mu <- propagate(mu, niter)
		flq <- FLQuant(c(t(flq)), dimnames = dimnames(mu))
		flq <- exp(mu + flq)
	}
	units(flq) <- units(object)
	return(flq)
})

#' @rdname genFLQuant-methods
# if nsim > 0 the simulate nsim times
setMethod("genFLQuant", "submodel", 
	function(object, type = c("link", "response"), nsim = 0) {
	  type <- match.arg(type)
      # make empty FLQuant
      flq <- FLQuant(
                  matrix(NA,
                     nrow = object@range["max"] - object@range["min"] + 1,
                     ncol = object@range["maxyear"] - object@range["minyear"] + 1),
                     dimnames = list(age = object@range["min"]:object@range["max"],
                                     year = object@range["minyear"]:object@range["maxyear"])
                  )
      df <- as.data.frame(flq)
      # simulate from submodel?
      if (nsim > 0) {
        object <- simulate(object, nsim = nsim)
      }
      # this should have 2 dimensions!
      b <- as(coef(object), "matrix")
      # get design matrix
      X <- model.matrix(formula(object), df)
      # predict accross all iters (if dimensions don't match then coefs are the wrong length!)
      pred <- X %*% b      
      # add into flq
      flq <- propagate(flq, max(1, nsim))
      flq[] <- as(pred, "vector")
      # transform if asked
      if (type == "response") {
      	object@linkinv(flq)
      } else {
        flq
      }
    }
 )



#' Methods to generate FLIndex objects
#' @description This method produces an \code{FLIndex} object by using the \code{genFLQuant} method.
#' @param object an \code{FLIndex} object
#' @param cv the coefficient of variation
#' @param niter the number of iterations to be generated
#' @param ... additional argument list that might not ever be used.
#' @return an FLIndex
#' @docType methods
#' @rdname genFLIndex-methods
#' @aliases genFLIndex genFLIndex-methods

setGeneric("genFLIndex", function(object, ...) standardGeneric("genFLIndex"))

#' @rdname genFLIndex-methods
setMethod("genFLIndex", c("FLQuant"), function(object, cv = 0.2, niter = 250) {
      # use log transform, to be expanded on later versions
      mu <- log(object)
      
#      if(method == "ac") {
        Rho <- cor(t(mu[drop = TRUE]))
        flq <- mvrnorm(niter * dim(mu)[2], rep(0, nrow(Rho)), log(cv^2+1) * Rho)
        mu <- propagate(mu, niter)
        flq <- FLQuant(c(t(flq)), dimnames = dimnames(mu))
        flq <- exp(mu + flq)
#      }
      return(flq)
})

