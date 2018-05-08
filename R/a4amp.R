mp <- function(opModel, obsModel, impModel, ctrl.mp, mpPars, scenario="test", tracking){

	#============================================================
	# prepare the om
	stk.om <- stock(opModel)	
	name(stk.om) <- scenario
	sr.om <- sr(opModel)
	sr.om.res <- residuals(sr.om)
	sr.om.res.mult <- sr.om@logerror
	fy <- mpPars$fy # final year
	y0 <- mpPars$y0 # initial data year
	dy <- mpPars$dy # final data year
	iy <- mpPars$iy # initial year of projection (also intermediate)
	nsqy <- mpPars$nsqy # number of years to compute status quo metrics
	ny <- fy - iy + 1 # number of years to project from intial year
	vy <- ac(iy:fy) # vector of years to be projected

	# start tracking
	if (missing(tracking)) 
		tracking <- FLQuant(NA, dimnames=list(metric=c("Fperc", "convergence", "Fhcr", "Implementation", "IEM", "FleetDyn","OM.f", "OM.ssb", "OM.catch"), year=c(iy-1,vy), iter=1:it))
	tracking["Implementation", ac(iy)] <- catch(stk.om)[,ac(iy)]
	if (!is.null(mpPars$seed)) set.seed(mpPars$seed)
  
	#============================================================
	# go fish
	for(i in vy[-length(vy)]){
		gc()
		ay <- an(i)
		cat(i, " > ")
		vy0 <- 1:(ay-y0) # data years (positions vector) - one less than current year
		sqy <- ac((ay-1):(ay-nsqy)) # years for status quo computations 
		
		tracking["OM.f", ac(ay-1)] <- fbar(stk.om)[,ac(ay-1)]    
		tracking["OM.ssb", ac(ay-1)] <- ssb(stk.om)[,ac(ay-1)]    
		tracking["OM.catch", ac(ay-1)] <- catch(stk.om)[,ac(ay-1)]    
		
		#==========================================================
		# OEM
		#----------------------------------------------------------
		# function o()
		ctrl.oem <- args(obsModel)
		ctrl.oem$method <- method(obsModel)
		ctrl.oem$deviances <- deviances(obsModel)
		ctrl.oem$observations <- observations(obsModel)
		ctrl.oem$stk <- stk.om
		ctrl.oem$vy0 <- vy0
		ctrl.oem$ay <- ay
		ctrl.oem$tracking <- tracking
		ctrl.oem$ioval <- list(iv=list(t1=flsval), ov=list(t1=flsval, t2=flival))
#		if (!is.null(ctrl.mp$ctrl.oem)){
#			o.out <- do.call("a4ampDispatch", ctrl.oem)
#		} else {
#      		ctrl.oem$method <- 'perfectInfo.oem'
#			o.out <- do.call("a4ampDispatch", ctrl.oem)
#		}
		o.out <- do.call("a4ampDispatch", ctrl.oem)
		stk0 <- o.out$stk
		idx0 <- o.out$idx
		observations(obsModel) <- o.out$observations
		tracking <- o.out$tracking

		#==========================================================
		# MP
		#----------------------------------------------------------
		# Estimator of stock statistics
		# function f()
		if (!is.null(ctrl.mp$ctrl.f)){
			ctrl.f <- ctrl.mp$ctrl.f
			ctrl.f$stk <- stk0
			ctrl.f$idx <- idx0
			ctrl.f$tracking <- tracking
			ctrl.f$ioval <- list(iv=list(t1=flsval, t2=flival), ov=list(t1=flsval))
			out.assess <- do.call("a4ampDispatch", ctrl.f)
			stk0 <- out.assess$stk
			tracking <- out.assess$tracking
		}
		tracking["Fperc",ac(ay)] <- fdy <- fbar(stk0)[,ac(ay-1)]

		#----------------------------------------------------------
		# HCR parametrization
		# function x()
		if (!is.null(ctrl.mp$ctrl.x)){
			ctrl.x <- ctrl.mp$ctrl.x
			ctrl.x$stk <- stk0
			ctrl.x$ay <- ay
			ctrl.x$iy <- iy
			ctrl.x$tracking <- tracking
			if(exists("hcrpars")) ctrl.x$hcrpars <- hcrpars
			ctrl.x$ioval <- list(iv=list(t1=flsval), ov=list(t1=flpval))
			out <- do.call("a4ampDispatch", ctrl.x)
			hcrpars <- out$hcrpars
			tracking <- out$tracking
		}

		#----------------------------------------------------------
		# HCR
		# function h()
		if (!is.null(ctrl.mp$ctrl.h)){
			ctrl.h <- ctrl.mp$ctrl.h
			ctrl.h$stk <- stk0
			ctrl.h$ay <- ay
			ctrl.h$tracking <- tracking
			if(exists("hcrpars")) ctrl.h$hcrpars <- hcrpars
			ctrl.h$ioval <- list(iv=list(t1=flsval), ov=list(t1=flfval))
			out <- do.call("a4ampDispatch", ctrl.h)
			ctrl <- out$ctrl
			tracking <- out$tracking
		} else {
			ctrl <- getCtrl(yearMeans(fbar(stk0)[,sqy]), "f", ay+1, it)
		}
		tracking["Fhcr", ac(ay)] <- ctrl@trgtArray[ac(ay+1),"val",]
		
		#----------------------------------------------------------
		# Management Implementation
		# function k()
		if (!is.null(ctrl.mp$ctrl.k)){
			ctrl.k <- ctrl.mp$ctrl.k
			ctrl.k$ctrl <- ctrl
			ctrl.k$stk <- stk0
			ctrl.k$ay <- ay
			ctrl.k$tracking <- tracking
			ctrl.k$ioval <- list(iv=list(t1=flsval, t2=flfval), ov=list(t1=flfval))
			out <- do.call("a4ampDispatch", ctrl.k)
			ctrl <- out$ctrl
			tracking <- out$tracking
			tracking["Implementation", ac(ay)] <- ctrl@trgtArray[ac(ay+1),"val",]
		} else {
			tracking["Implementation", ac(ay)] <- tracking["Fhcr", ac(ay+1)]
		}

		#----------------------------------------------------------
		# Technical measures
		# function w()
		if (!is.null(ctrl.mp$ctrl.w)){
			ctrl.w <- ctrl.mp$ctrl.w
			ctrl.w$stk <- stk0
			ctrl.w$sqy <- sqy
			ctrl.w$tracking <- tracking
			ctrl.w$ioval <- list(iv=list(t1=flsval), ov=list(t1=flqval))
			out <- do.call("a4ampDispatch", ctrl.w)
			attr(ctrl, "snew") <- out$flq
			tracking <- out$tracking
		}

		#==========================================================
		# IEM
		#----------------------------------------------------------
		# function l()
		if (!is.null(ctrl.mp$ctrl.l)){
			ctrl.l <- ctrl.mp$ctrl.l
			ctrl.l$ctrl <- ctrl
			ctrl.l$tracking <- tracking
			ctrl.l$ioval <- list(iv=list(t1=flfval), ov=list(t1=flfval))
			out <- do.call("a4ampDispatch", ctrl.l)
			ctrl <- out$ctrl
			tracking <- out$tracking
		}
		tracking["IEM",ac(ay)] <- ctrl@trgtArray[,"val",]

		#==========================================================
		# OM
		#----------------------------------------------------------
		# fleet dynamics/behaviour
		# function j()
		if (exists(fleetBehaviour(opModel))){
			ctrl.j <- args(fleetBehaviour(opModel))
			ctrl.j$method <- method(fleetBehaviour(opModel))
			ctrl.j$tracking <- tracking
			ctrl.j$ctrl <- ctrl
			ctrl.j$ioval <- list(iv=list(t1=flfval), ov=list(t1=flfval))
			out <- do.call("a4ampDispatch", ctrl.j)
			ctrl <- out$ctrl
			tracking <- out$tracking
		}
		tracking["FleetDyn",ac(ay)] <- ctrl@trgtArray[,"val",]

		#----------------------------------------------------------
		# stock dynamics and OM projections
		# function g()
		if(!is.null(attr(ctrl, "snew"))) harvest(stk.om)[,ac(ay+1)] <- attr(ctrl, "snew")
		stk.om <- fwd(stk.om, ctrl=ctrl, sr=sr.om, sr.residuals = sr.om.res, sr.residuals.mult = sr.om.res.mult, maxF=2)

	}
    cat("\n")

	#============================================================
    attr(stk.om, "tracking") <- tracking
    list(opModel=stk.om, obsModel=obsModel)

}

# dispatching
flsval <- list(object="stk", test="!is(object, \"FLS\")", msg="\"stk must be of class FLStock\"")
flival <- list(object="idx", test= "!is(object, \"FLIndices\")", msg="\"idx must be of class FLIndices\"")
flpval <- list(object="hcrpars", test= "!is(object, \"FLPar\")", msg="\"hcrpars must be of class FLPar\"")
flfval <- list(object="ctrl", test= "!is(object, \"fwdControl\")", msg="\"ctrl must be of class fwdControl\"")
flqval <- list(object="flq", test= "!is(object, \"FLQuant\")", msg="\"flq must be of class FLQuant\"")

a4ampDispatch <- function(ioval, ...){
	args <- list(...)
	method <- args$method
	args$method <- NULL
	# checks in
	for(i in ioval$iv){
		object <- args[i$object]
		str <- paste("if(", i$test, ")", i$msg, sep=" ")
		eval(parse(text=str))
	}
	# dispatch
	out <- do.call(method, args)
	# checks out
	for(i in ioval$ov){
		object <- out[i$object]
		str <- paste("if(", i$test, ")", i$msg, sep=" ")
		eval(parse(text=str))
	}
	out
}

