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

