#==================================================================== 
#    predict methods
#==================================================================== 

#' @title Predict methods for SCA
#' @description Predict methods for a4a stock assessment fits.
#' @name predict
#' @rdname predict-methods
#' @aliases predict,a4aFitSA-method
#' @examples
#' data(ple4)
#' data(ple4.index)
#' fmodel <- ~factor(age) + factor(year)
#' qmodel <- list(~factor(age))
#' fit1 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#' flqs <- predict(fit1)

setMethod("predict", signature(object = "a4aFitSA"),
  function(object) {
  	obj <- pars(object)
	# need to tag biomass indices
  	for(i in seq_along(index(object))) 
  		attr(obj@vmodel[[i+1]], "FLIndexBiomass") <- attr(obj@qmodel[[i]], "FLIndexBiomass") <- attr(index(object)[[i]], "FLIndexBiomass")
  	predict(obj)
  })

#' @rdname predict-methods
#' @aliases predict,SCAPars-method
setMethod("predict", signature(object = "SCAPars"),
  function(object) {

	sm <- stkmodel(object)
	qm <- qmodel(object)
	vm <- vmodel(object)
	# need to update centering to include stock centering
	for(i in 1:length(qm)) qm[[i]]@centering <- qm[[i]]@centering-sm@centering
	# run predict 
    lst <- list(
      stkmodel = predict.stkpars(sm),
      qmodel   = predict.submods(qm),
      vmodel   = predict.submods(vm)
    )
	# rec is being estimated from the srmodel so the n1model doesn't
	# have a recruitment estimate, need to update from rec predictions
	lst$stkmodel$ny1[1,1] <- lst$stkmodel$rec[,1]
	lst
})

predict.stkpars <- function(object) {
      ages <- range(object)["min"]:range(object)["max"]
      years <- range(object)["minyear"]:range(object)["maxyear"]
      cnames <- rownames(coef(object))
      df <- expand.grid(age = ages,
                        year = years)
      X <- getX(object @ fMod, df)
      b <- coef(object)[grep("fMod", cnames)]
      niter <- dim(b)[2] # reuse this for the others
      fit <- exp(c(X %*% b))
      harvest <- 
        FLQuant(array(fit, dim = c(length(ages), length(years), 1, 1, 1, niter), 
                      dimnames = list(age = ages, year = years, 
                                         unit = "unique", season = "all", area = "unique",
                                         iter = seq(niter))),
                units = "f")


      X <- getX(object @ n1Mod, data.frame(age = ages[-1]))
      b <- coef(object)[grep("n1Mod", cnames)]
      fit <- exp(c(rbind(NA, X %*% b)) + object @ centering)
      ny1 <-       
        FLQuant(array(fit, dim = c(length(ages), 1, 1, 1, 1, niter), 
                      dimnames = list(age = ages, year = years[1], 
                                         unit = "unique", season = "all", area = "unique",
                                         iter = seq(niter))))

	  # check if internal S/R was used
	  facs <- strsplit(as.character(object @ srMod)[length(object @ srMod)], "[+]")[[1]]
	  facs <- gsub("(^ )|( $)", "", facs) # remove leading and trailing spaces
	  a4as <- grepl(paste("(^",c("bevholt", "bevholtSV", "ricker","hockey","geomean"),"[(])", collapse = "|", sep = ""), facs)
	  if(a4as){
	      X <- getX(~factor(year), data.frame(year = years))
	      b <- coef(object)[grep("rMod", cnames)]
	      fit <- c(exp(c(X %*% b) + object @ centering))
	      rec <- FLQuant(array(fit, dim = c(1, length(years), 1, 1, 1, niter),
	      		dimnames = list(age = ages[1], year = years,
	      		unit = "unique", season = "all", area = "unique",
	      		iter = seq(niter))))

	  # FFFFFFFFFUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUCCCCCCCCCCCCCCCCCCCCCCCKKKKKKKKKKKKKKKKKKKKK
#		  bevholt <- function(a = ~ 1, b = ~ 1, CV = 0.5) {
#			if (CV <= 0) stop ("CV in stock recruit relationship cannot be less than zero")
#			list(srr = "bevholt", a = a, b = b, SPR0 = 1, srrCV = CV, ID = 1)
#		  }
#		  bevholtSV <- function(h = ~ 1, v = ~ 1, SPR0 = 1, CV = 0.5) {
#			if (CV <= 0) stop ("CV in stock recruit relationship cannot be less than zero")
#			list(srr = "bevholtSV", a = h, b = v, SPR0 = SPR0, srrCV = CV, ID = 5)
#		  }
#		  ricker <- function(a = ~ 1, b = ~ 1, CV = 0.5) {
#			if (CV <= 0) stop ("CV in stock recruit relationship cannot be less than zero")
#			list(srr = "ricker", a = a, b = b, SPR0 = 1, srrCV = CV, ID = 2)
#		  }
#		  hockey <- function(a = ~ 1, b = ~ 1, CV = 0.5) {
#			if (CV <= 0) stop ("CV in stock recruit relationship cannot be less than zero")
#			list(srr = "hockey", a = a, b = b, SPR0 = 1, srrCV = CV, ID = 3)
#		  }
#		  geomean <- function(a = ~ 1, CV = 0.5) {
#			if (CV <= 0) stop ("CV in stock recruit relationship cannot be less than zero")
#			list(srr = "geomean", a = a, b = ~ 1, SPR0 = 1, srrCV = CV, ID = 4)
#		  }
#		  # get that stuff
#	  	  srr <- eval(parse(text = facs))
#		  # a
#	      X <- getX(srr$a, data.frame(year = years))
#	      b <- coef(object)[grep("sraMod", cnames)]
#	      #fita <- c(exp(c(X %*% b) + object @ centering))
#	      fita <- X %*% b
#	      # b
#	      X <- getX(srr$b, data.frame(year = years))
#	      b <- coef(object)[grep("srbMod", cnames)]
#	      fitb <- X %*% b
#		  # need to compute predictions from S/R and weight with factor(year) ?
	  # END OF FUCK
	  
	  } else {
	      X <- getX(object @ srMod, data.frame(year = years))
	      b <- coef(object)[grep("rMod", cnames)]
	      fit <- c(exp(c(X %*% b) + object @ centering))
	      rec <- FLQuant(array(fit, dim = c(1, length(years), 1, 1, 1, niter),
	      		dimnames = list(age = ages[1], year = years,
	      		unit = "unique", season = "all", area = "unique",
	      		iter = seq(niter))))
	  }
      FLQuants(harvest = harvest, rec = rec, ny1 = ny1)
}

predict.submods <- function(object, ...) {
	lst <- lapply(object, function(x){
		if(isTRUE(attr(x, "FLIndexBiomass"))) ages <- "all" else ages <- range(x)["min"]:range(x)["max"]
		years <- try(range(x)["minyear"]:range(x)["maxyear"], silent=TRUE)
		if(is(years, "try-error")) years <- 1
		df <- expand.grid(age = ages, year = years)
		X <- getX(x @ Mod, df)
		b <- coef(x)
		niter <- dim(b)[2]
		fit <- exp(c(X %*% b) + x @ centering)
		FLQuant(array(fit, dim = c(length(ages), length(years), 1, 1, 1, niter), dimnames = list(age = ages, year = years, unit = "unique", season = "all", area = "unique", iter = seq(niter))))
	})
	lst
}

