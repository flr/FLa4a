###############################################################################
# EJ&CM(2012821)
# Methods to add uncertainty to FLQuant objects
# NOTE #1:
# NOTE #2:
###############################################################################

#' Methods to generate FLStock objects
#'
#' Some additional details about this S4 generic and its methods.
#' The extra blank line between this section and the title is
#' critical for roxygen2 to differentiate the title from the
#' description section.
#'
#' @param object an FLQuant
#' @param rec an FLQuant
#' @param catch.n an FLQuant
#' @param harvest an FLQuant
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname genFLStock-methods
#'
#' @examples
#' data(ple4)
#' sim.F <- genFLQuant(harvest(ple4))
setGeneric("genFLStock", function(object, R, C, F, ...) standardGeneric("genFLStock"))

#' @rdname genFLStock-methods
#' @aliases genFLStock,FLStock,FLQuant,FLQuant,missing-method
setMethod("genFLStock", c("FLStock", "FLQuant", "FLQuant", "missing"), function(object, R, C, F, ...){
	cat("Not implemented yet\n")
})

#' @rdname genFLStock-methods
#' @aliases genFLStock,FLStock,FLQuant,FLQuant,missing-method
setMethod("genFLStock", c("FLStock", "missing", "FLQuant", "FLQuant"), function(object, R, C, F, ...){
	cat("Not implemented yet\n")
})

#' @rdname genFLStock-methods
#' @aliases genFLStock,FLStock,FLQuant,FLQuant,missing-method
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

#' Methods to add uncertainty to FLQuant objects
#'
#' Some additional details about this S4 generic and its methods.
#' The extra blank line between this section and the title is
#' critical for roxygen2 to differentiate the title from the
#' description section.
#'
#' @param object an FLQuant
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname getAcor-methods
#'
#' @examples
#' data(ple4)
#' genFLQuant(harvest(ple4), method = "ac")
setGeneric("getAcor", function(object, ...) standardGeneric("getAcor"))

#' @rdname getAcor-methods
#' @aliases getAcor,FLQuant-method
setMethod("getAcor", c("FLQuant"), function(object, tf=log, ...) {
		mu <- log(object)
		Rho <- cor(t(mu[drop = TRUE]))
		return(Rho)
})

#' Methods to genetate FLStock objects
#'
#' Some additional details about this S4 generic and its methods.
#' The extra blank line between this section and the title is
#' critical for roxygen2 to differentiate the title from the
#' description section.
#'
#' @param object an FLQuant
#' @param rec an FLQuant
#' @param catch.n an FLQuant
#' @param harvest an FLQuant
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname genFLQuant-methods
#'
#' @examples
#' data(ple4)
#' sim.F <- genFLQuant(harvest(ple4), method = "ac")
setGeneric("genFLQuant", function(object, ...) standardGeneric("genFLQuant"))

#' @rdname genFLQuant-methods
#' @aliases genFLQuant,FLQuant-method
setMethod("genFLQuant", c("FLQuant"), function(object, cv = 0.2, method = "ac", niter = 250) {
  # use log transform, to be expanded on later versions
	mu <- log(object)
	if(method == "ac") {
		Rho <- cor(t(mu[drop = TRUE]))
		flq <- mvrnorm(niter * dim(mu)[2], rep(0, nrow(Rho)), log(cv^2+1) * Rho)
		mu <- propagate(mu, niter)
		flq <- FLQuant(c(t(flq)), dimnames = dimnames(mu))
		flq <- exp(mu + flq)
	}
#	if(method == "rw") {
#		n.ages  <- dim(mu)[1]
#	  	n.years <- dim(mu)[2]
#		n <- n.ages * n.years
#		# set up lcs to extract posterior means
#		B = diag(n)
#		B[B==0] <- NA
#		lcs = inla.make.lincombs(Predictor = B)
#		# treat mu as a GMRF model - 
#		# an independant random walk for each age (same variance accross ages)
#		form <- x ~ f(years, model = 'rw1', replicate = ages)
#		data <- list(x = c(mu), years = rep(1:n.years, each = n.ages), ages  = rep(1:n.ages, n.years))
#		result <- inla(form, data = data, control.predictor = list(compute = TRUE), 
#                      lincomb = lcs, control.inla = list(lincomb.derived.correlation.matrix = TRUE))
#		# the covariance of the fitted RW
#		RW.cov <- result $ misc $ lincomb.derived.correlation.matrix
#		# two options for the mean:
#		#  1) use the mean estimate of RW process from INLA
#		#     - this is potentially very smooth and lacking in strucure
#		#	mu.hat <- result $ summary.linear.predictor $ mean
#		#	flq <- mvrnorm(niter, mu.hat, cv^2 * RW.cov)
#		#  2) use the original data and add the noise to that
#		#  2 is more consistent with ac method and always maintains data structure
#		flq <- exp(mvrnorm(niter, c(mu), cv^2 * RW.cov))
#		flq <- FLQuant(c(t(flq)), dimnames = dimnames(propagate(mu, niter)))
#	}
	units(flq) <- units(object)
	return(flq)
})

#' Methods to create index FLQuants from stock.n
#'
#' Some additional details about this S4 generic and its methods.
#' The extra blank line between this section and the title is
#' critical for roxygen2 to differentiate the title from the
#' description section.
#'
#' @param object an FLQuant containing stock numbers at age
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant?
#' 
#' @seealso \code{\link{stock.n}} and \code{\link{FLQuant}}
#' 
#' @export
#' @docType methods
#' @rdname genFLIndex-methods
#'
#' @examples
#' data(ple4)
#' genFLQuant(harvest(ple4), method = "ac")
setGeneric("genFLIndex", function(object, ...) standardGeneric("genFLIndex"))

#' @rdname genFLIndex-methods
#' @aliases genFLIndex,FLQuant-method
setMethod("genFLIndex", c("FLQuant"), 
    function(object, cv = 0.2, 
                     catchability = "flat", 
                     niter = 250) {
      # use log transform, to be expanded on later versions
      mu <- log(object)
      
      if(method == "ac") {
        Rho <- cor(t(mu[drop = TRUE]))
        flq <- mvrnorm(niter * dim(mu)[2], rep(0, nrow(Rho)), cv^2 * Rho)
        mu <- propagate(mu, niter)
        flq <- FLQuant(c(t(flq)), dimnames = dimnames(mu))
        flq <- exp(mu + flq)
      }

      
      return(flq)
})

