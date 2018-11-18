
#--------------------------------------------------------------------
setGeneric("rz", function(object, ...) standardGeneric("rz"))
setMethod("rz", "FLIndex", function(object){
	flq <- index(object)
	flq[flq==0] <- min(flq[flq!=0])
	index(object) <- flq
	object
})

setMethod("rz", "FLIndices", function(object){
	lapply(object, rz)
})

#--------------------------------------------------------------------
yearCumsum <- function(x, ...){
	x[] <- apply(x, c(1,3:6), cumsum)
	return(x)
}

#--------------------------------------------------------------------
yearDiffPerc <- function(x, ...){
	#x[,-1] <- x[,-1]/x[,-ncol(x)]-1
	x[,-1] <- x[,-1]/c(x[,1])-1
	x[,1] <- 0
	return(x)
}

#--------------------------------------------------------------------
as.table.FLQuants <- function(x){
	x <- mcf(x)
	df0 <- as.data.frame(do.call("cbind", lapply(x, c)))
	row.names(df0) <- dimnames(x)[[2]]
	df0
}

#--------------------------------------------------------------------
iterQuantiles <- function(x, prob=0.5, ...){
	return(apply(x, c(1:5), quantile, prob=prob, na.rm = FALSE))
}

#--------------------------------------------------------------------
an <- function(x, ...) as.numeric(x, ...)

cbind.FLQuant <- function(x, y){
	lst <- mcf(list(x, y))
	res <- lst[[1]]
	res[,dimnames(y)[[2]]] <- y
	res
} 

#--------------------------------------------------------------------
getFmod <- function(stk, model="interaction", dfm=c(2/3, 1/2)){
	dis <- dims(stk)
	KY=floor(dfm[1] * dis$year)
	KA=ceiling(dfm[2] *dis$age)
	if(model=="interaction"){
		frm <- substitute(~ te(age, year, k=c(KA, KY)), list(KA=KA, KY=KY))
	} else if (model=="separable"){
		frm <- substitute(~ s(age, k = KA) + s(year, k = KY), list(KA=KA, KY=KY))
	}
	as.formula(frm)
}

#--------------------------------------------------------------------
getQmod <- function(idx){
	lds <- lapply(idx, dims)
	lds <- lapply(lds, function(x){
		if(x$age<=3){
			frm <- ~factor(age)
		} else {
			frm <- substitute(~s(age, k=KA), list(KA=ceiling(3/4 * x$age)))
		} 
		as.formula(frm)
	})
	lds		
}

#--------------------------------------------------------------------
getCtrl <- function(values, quantity, years, it){
	dnms <- list(iter=1:it, year=years, c("min", "val", "max"))
	arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
	arr0[,,"val"] <- unlist(values)
	arr0 <- aperm(arr0, c(2,3,1))
	ctrl <- fwdControl(data.frame(year=years, quantity=quantity, val=NA))
	ctrl@trgtArray <- arr0
	ctrl
}

#--------------------------------------------------------------------
iem <- function(iem, ctrl){
	iem0 <- iem
  # strip out fun and multiplicative from list of parameters, 
  # to use do.call() later
	iem0$fun <- iem0$multiplicative <- NULL
  # eh? number of non NA values in target. But if there are NAs then we need to
  # know their position for the *ctrl@trgtArray later
	iem0$n <- sum(!is.na(ctrl@trgtArray))
	if(iem$multiplicative){
		ctrl@trgtArray <- do.call(iem$fun, iem0)*ctrl@trgtArray
	} else {
		ctrl@trgtArray <- do.call(iem$fun, iem0) + ctrl@trgtArray
	}
	ctrl
}

#--------------------------------------------------------------------
getF <- function(trgt, fc, mxy, ay, object, multiplicative=TRUE, correction=TRUE){
    # F trajectory to Ftrg
    # note f is estimated for ay-1
    # ay is assessment year, also the intermediate year. For this year the
    # decisions were takenon ay-1 but we don't have observations of its
    # implementation
    # ry is the decrease needed in ay
    # ry1 is the decrease needed in ay+1
    # object has to be a FLQuant

    # work in medians and add variability in the end
    object0 <- object
    object <- iterMedians(object)
    fc <- iterMedians(fc)

    if(multiplicative){
        ry1 <- (trgt/object[,ac(ay)])^(1/(mxy-ay))
        object[,ac(ay+1)] <- ifelse((ay+1) <= mxy, object[,ac(ay)]*ry1, trgt)
        if(correction){
            cny <- mxy-ay+1
            ry <- ifelse((ay+1) < mxy, (trgt/fc)^(1/(cny)), 1)
            object[,ac(ay+1)] <- object[,ac(ay+1)]*(object[,ac(ay+1)]/(ry*fc))
        }
    } else {
        stop("Only additive implemented\n")
      #  ry <- (trgt-fc)/(mxy-ay+1)
      #  ry1 <- 2*ry
      #  object[,ac(ay+1)] <- fc+ry1
      #  if(correction){
      #      object[,ac(ay+1)] <- object[,ac(ay+1)]-(object[,ac(ay)]-(fc+ry))
      #  }
        }
    object0[,ac(ay+1)] <- object[,ac(ay+1)]/object[,ac(ay)]*object0[,ac(ay)]

    object0
}

#--------------------------------------------------------------------
# effort to f to model hyperstability or hyperdepletion
e2f <- function(object, alpha=missing, beta, maxF=2.0){
    # object must be fwdControl
    if(missing(alpha)) alpha <- maxF^(1-beta) # linear meets curve at maxF
    object@trgtArray[object@target[,"quantity"]=="f",,] <- alpha *
      object@trgtArray[object@target[,"quantity"]=="f",,]^beta
    object
}

#--------------------------------------------------------------------
# ar1rlnorm {{{
ar1rlnorm <- function(rho, years, iters=1, margSD=0.6) {
  n <- length(years)
  rhosq <- rho ^ 2

  res <- matrix(rnorm(n*iters, mean=0, sd=margSD), nrow=n, ncol=iters)
  res <- apply(res, 2, function(x) {
    for(i in 2:n)
    x[i] <- sqrt(rhosq) * x[i-1] + sqrt(1-rhosq) * x[i]
    return(exp(x))
    }
  )
  return(FLQuant(array(res, dim=c(1,n,1,1,1,iters)),
    dimnames=list(year=years, iter=seq(1, iters))))
  } # }}}
  
  
  # add uncertainty using RW or correlation matrix
  setMethod("genFLQuant", c("FLQuant"),
    function(object, cv = 0.2, method = "rw", niter = 250) {
  # use log transform, to be expanded on later versions
	mu <- log(object)
	if(method == "ac") {
		Rho <- cor(t(mu[drop = TRUE]))
		flq <- mvrnorm(niter * dim(mu)[2], rep(0, nrow(Rho)), cv^2 * Rho)
		mu <- propagate(mu, niter)
		flq <- FLQuant(c(t(flq)), dimnames = dimnames(mu))
		flq <- exp(mu + flq)
	}
	if(method == "rw") {
		n.ages  <- dim(mu)[1]
	  	n.years <- dim(mu)[2]
		n <- n.ages * n.years
		# set up lcs to extract posterior means
		B = diag(n)
		B[B==0] <- NA
		lcs = inla.make.lincombs(Predictor = B)
		# treat mu as a GMRF model - 
		# an independant random walk for each age (same variance accross ages)
		form <- x ~ f(years, model = 'rw1', replicate = ages)
		data <- list(x = c(mu), years = rep(1:n.years, each = n.ages), ages  = rep(1:n.ages, n.years))
		result <- inla(form, data = data, control.predictor = list(compute = TRUE), 
                       lincomb = lcs, control.inla = list(lincomb.derived.correlation.matrix = TRUE))
		# the covariance of the fitted RW
		RW.cov <- result $ misc $ lincomb.derived.correlation.matrix
		# two options for the mean:
		#  1) use the mean estimate of RW process from INLA
		#     - this is potentially very smooth and lacking in strucure
		#	mu.hat <- result $ summary.linear.predictor $ mean
		#	flq <- mvrnorm(niter, mu.hat, cv^2 * RW.cov)
		#  2) use the original data and add the noise to that
		#  2 is more consistent with ac method and always maintains data structure
		flq <- exp(mvrnorm(niter, c(mu), cv^2 * RW.cov))
		flq <- FLQuant(c(t(flq)), dimnames = dimnames(propagate(mu, niter)))
	}
	units(flq) <- units(object)
	return(flq)
})


#--------------------------------------------------------------------
# add uncertainty to FLStock using quants
setGeneric("au", function(object, R, C, F, ...) standardGeneric("au"))

setMethod("au", c("FLStock", "FLQuant", "missing", "FLQuant"), 
  function(object, R, C, F, ...){
	# requires checking dimensions
	if(!identical(dim(catch.n(object))[-c(1,6)], dim(R)[-c(1,6)]))
    stop("Recruitment vector must have consistent dimensions with the stock object")
	if(!identical(dim(catch.n(object))[-6]    , dim(F)[-6])) 
    stop("Harvest matrix must have consistent dimensions with the stock object")
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
	#Ns <- FLCohort(R[rep(1,nages)])
	#Ns <- Ns*exp(-Z)
	Ns <- R[rep(1,nages)]
	# stupid hack to get the right dimensions
	Z[,dimnames(Ns)[[2]]] <- Ns*exp(-Z[,dimnames(Ns)[[2]]])
	Ns <- as(Z, "FLQuant")

	# Update object
	stock.n(object) <- flq
	# [R]
	stock.n(object)[1] <- R
	# [N]
	stock.n(object)[-1,-1] <- Ns[-nages,-nyrs] 
	# plus group
	stock.n(object)[nages,-1] <- Ns[nages,-nyrs] + stock.n(object)[nages,-1]
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

