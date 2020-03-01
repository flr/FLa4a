#' @title deprecated
#' @name deprecated
#' @docType methods
#' @rdname deprecated
#' @template dots
#' @description Deprecated methods.
#' @aliases a4aSCA
a4aSCA <- function(...){
	stop("The method \"a4aSCA\" was removed, please use \"sca\", which now gives the user assess to the same arguments \"a4aSCA\".")
}

#' @title Default sub-models
#' @name defaultSubModels
#' @docType methods
#' @rdname defaultsubmodels
#' @description Methods to create formulas for sub-models. The sub-models are set automagically using defaults.
#' @param stock an FLStock object
#' @param indices an FLIndices object
#' @param dfm numeric vector with the data points fraction to be used to set the spline ks.
#' @return a FLStock object
#' @aliases defaultFmod
defaultFmod <- function(stock, dfm=c(0.5, 0.7)){
	dis <- dims(stock)
	KY=floor(dfm[1] * dis$year)
	KA=ceiling(dfm[2] *dis$age)
	if (KA >= 3) {
		KA <- min(max(3, KA), 6)
		KB <- min(max(3, KA), 10)
	    fmodel <- formula(paste("~ te(age, year, k = c(", KA,",", KY,"), bs = 'tp') + s(age, k=", KB, ")"))
	  } else {
		fmodel <- formula(paste("~ age + s(year, k = ", KY,")"))
	  }
	fmodel
}

#' @rdname defaultsubmodels
#' @aliases defaultQmod
defaultQmod <- function(indices, dfm=0.6){
	lds <- lapply(indices, dims)
	lds <- lapply(lds, function(x){
		if(x$age==1){
			frm <- ~1
		} else if(x$age>1 & x$age<=3){
			frm <- ~factor(age)
		} else {
			frm <- substitute(~s(age, k=KA), list(KA=min(ceiling(dfm * x$age), 6)))
		}
		as.formula(frm)
	})
	lds
}

#' @rdname defaultsubmodels
#' @aliases defaultN1mod
defaultN1mod <- function(stock){
  dis <- dims(stock)
  if(dis$age==1){
	frm <- ~1
  } else if(dis$age>1 & dis$age<=3){
	frm <- ~factor(age)
  } else {
	frm <- ~ s(age, k = 3)
  }
  as.formula(frm)
}

#' @rdname defaultsubmodels
#' @aliases defaultVmod
defaultVmod <- function(stock, indices){
  vmodel  <- lapply(seq(length(indices) + 1), function(i) ~ 1)
  dis <- dims(stock)
  if(dis$age==1){
	frm <- ~1
  } else if(dis$age>1 & dis$age<=3){
	frm <- ~factor(age)
  } else {
	frm <- ~ s(age, k = 3)
  }
  vmodel[[1]] <- frm
  vmodel
}

#' @rdname defaultsubmodels
#' @aliases defaultSRmod
defaultSRmod <- function(stock){~factor(year)}

#' @title Breakpoints
#' @name breakpts
#' @rdname breakpts
#' @description Method to set breakpoints in submodels
#' @param var a \code{numeric} object that defines the variable to be "broken"
#' @param breaks a \code{numeric} object that defines the breakpoints
#' @template dots
#' @return a \code{factor} with levels according to the defined breaks
#' @aliases breakpts breakpts-methods
setGeneric("breakpts", function(var, ...) standardGeneric("breakpts"))
#' @rdname breakpts
setMethod("breakpts", "numeric", function(var, breaks, ...) {
  if (min(var, na.rm = TRUE) < min(breaks)) breaks <- c(min(var, na.rm = TRUE) - 1, breaks)
  if (max(var, na.rm = TRUE) > max(breaks)) breaks <- c(breaks, max(var, na.rm = TRUE))
  label <- paste0("(",breaks[-length(breaks)], ",", breaks[-1], "]")
  cut(var, breaks = breaks, label = label)
})
