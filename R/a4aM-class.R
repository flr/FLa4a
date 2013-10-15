# validity: iters must be 1 or n in all models
valida4aM <- function(object){

	lvl <- object@level	
	shp <- object@shape
	trd <- object@trend

	v1 <- c(level=pars2dim(lvl), shape=pars2dim(shp), trend=pars2dim(trd))
	v2 <- c(level=length(dimnames(lvl@params)[["iter"]]), shape=length(dimnames(shp@params)[["iter"]]), trend=length(dimnames(trd@params)[["iter"]]))
	v2 <- v2==1 | v2==max(v2)
	df0 <- data.frame(dim2iter=v1,iter1orn=v2)
	if(sum(!df0)>0) return("Object is not valid. Check that all params have 2 dimensions with the second being named \"iter\" and that the number of iters in all models are 1 or N")
	return(TRUE)
}

#' Natural mortality class
#'
#' @section Slot: 
#' \itemize{
#'	\item \code{shape} the shape of M by age.
#'	\item \code{level} the mean level of M in a range of ages, which will be used to scale the \code{shape}.
#'	\item \code{trend} the yearly trend in M.
#' '}
#' 
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export

setClass("a4aM",
        representation(
				"FLComp",
                shape = "FLModelSim",
                level = "FLModelSim",
				trend = "FLModelSim"),
        prototype = prototype(
				name	= character(0),
				desc	= character(0),
				range	= c(min=0, max=0, plusgroup=0, minyear=0, maxyear=0, minmbar=0, maxmbar=0),
                shape = FLModelSim(),
                level = FLModelSim(),
                trend = FLModelSim()),
		validity=valida4aM
)

#' @title Constructor for the \code{a4aM} class 
#'
#' @description Constructor for the natural mortality class using \code{FLModelSim} objects.
#' @param shape a \code{FLModelSim} object
#' @param level a \code{FLModelSim} object
#' @param trend a \code{FLModelSim} object
#' @return a \code{a4aM} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
#' @examples
#' mod1 <- FLModelSim(model=~exp(-age-0.5))
#' mod2 <- FLModelSim(model=~1.5*k, params=FLPar(k=0.4))
#' m1 <- a4aM(shape=mod1, level=mod2)
setGeneric("a4aM", function(object, ...) standardGeneric("a4aM"))
setMethod("a4aM", signature(object="missing"),
  function(...) {
    # empty
  	if(missing(...)){
	  	new("a4aM")
    # or not
  	} else {
      args <- list(...)
	  args$Class <- 'a4aM'
      do.call("new", args)
	  }
  }
)

#' @title shape accessor
#' @description accessor method for \code{a4aM} object's slot \code{shape}.
#' @param object a \code{a4aM} object
#' @return a \code{FLModelSim} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export

setGeneric("shape", function(object, ...) standardGeneric("shape"))
setMethod("shape", "a4aM", function(object) object@shape)

#' @title shape replacement
#' @description replacement method for \code{a4aM} object's slot \code{shape}.
#' @param object a \code{a4aM} object
#' @param value a \code{FLModelSim} object
#' @return a \code{a4aM} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export

setGeneric("shape<-", function(object,value) standardGeneric("shape<-"))
setReplaceMethod("shape", "a4aM", function(object, value){
	if(all.equal(is(value), is(object@shape))) object@shape <- value
	object
})

#' @title level accessor
#' @description accessor method for \code{a4aM} object's slot \code{level}.
#' @param object a \code{a4aM} object
#' @return a \code{FLModelSim} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export

setGeneric("level", function(object, ...) standardGeneric("level"))
setMethod("level", "a4aM", function(object) object@level)

#' @title level replacement
#' @description replacement method for \code{a4aM} object's slot \code{level}.
#' @param object a \code{a4aM} object
#' @param value a \code{FLModelSim} object
#' @return a \code{a4aM} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export

setGeneric("level<-", function(object,value) standardGeneric("level<-"))
setReplaceMethod("level", "a4aM", function(object, value){
	if(all.equal(is(value), is(object@level))) object@level <- value
	object
})

#' @title trend accessor
#' @description accessor method for \code{a4aM} object's slot \code{trend}.
#' @param object a \code{a4aM} object
#' @return a \code{FLModelSim} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export

setGeneric("trend", function(object, ...) standardGeneric("trend"))
setMethod("trend", "a4aM", function(object) object@trend)

#' @title trend replacement
#' @description replacement method for \code{a4aM} object's slot \code{trend}.
#' @param object a \code{a4aM} object
#' @param value a \code{FLModelSim} object
#' @return a \code{a4aM} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export

setGeneric("trend<-", function(object,value) standardGeneric("trend<-"))
setReplaceMethod("trend", "a4aM", function(object, value){
	if(all.equal(is(value), is(object@trend))) object@trend <- value
	object
})

#' @title M range
#' @description method to extract from \code{a4aM} objects the age range for M.
#' @param object a \code{a4aM} object
#' @return a \code{vector} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export

setGeneric("rngmbar", function(object, ...) standardGeneric("rngmbar"))
setMethod("rngmbar", "a4aM", function(object){
	object@range[c("minmbar","maxmbar")]
})

#' @title M range replacement
#' @description method to replace \code{a4aM} object's age range for M.
#' @param object a \code{a4aM} object
#' @param value a \code{vector} with max and min age range 
#' @return a \code{a4aM} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("rngmbar<-", function(object,value) standardGeneric("rngmbar<-"))
setReplaceMethod("rngmbar", "a4aM", function(object, value){
	object@range[c("minmbar","maxmbar")] <- value
	object
})


#' @title M age vector
#' @description method to extract from \code{a4aM} objects the vector of ages for the mean M.
#' @param object a \code{a4aM} object
#' @return a \code{vector} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("vecmbar", function(object, ...) standardGeneric("vecmbar"))
setMethod("vecmbar", "a4aM", function(object){
	rng <- object@range[c("minmbar","maxmbar")]
	rng[1]:rng[2]
})

