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

#' @title Natural mortality class
#' @name a4aM
#' @rdname a4aM-class
#' @template ClassDescription
#' @section Slot: 
#' \describe{
#'	\item{\code{shape}}{the shape of M by age}
#'	\item{\code{level}}{the mean level of M over a range of ages, which will be used to scale the \code{shape}}
#'	\item{\code{trend}}{the yearly trend in M}
#' }
#' @aliases a4aM-class

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

#' @rdname a4aM-class
#' @aliases a4aM a4aM-methods a4aM,missing-method
#' @template Accessors
#' @template Constructors
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

#' @rdname a4aM-class
#' @aliases show,a4aM-method
setMethod("show", "a4aM", function(object){
	cat(paste("a4aM object", object@name, ":\n", sep=""))
	cat(paste("  shape: ", deparse(object@shape@model), "\n", sep=""))
	cat(paste("  level: ", deparse(object@level@model), "\n", sep=""))
	cat(paste("  trend: ", deparse(object@trend@model), "\n", sep=""))
})

#' @rdname a4aM-class
#' @aliases shape shape-methods shape,a4aM-method
setGeneric("shape", function(object, ...) standardGeneric("shape"))
setMethod("shape", "a4aM", function(object) object@shape)

#' @rdname a4aM-class
#' @aliases shape<- shape<--methods shape<-,a4aM-method
setGeneric("shape<-", function(object,value) standardGeneric("shape<-"))
setReplaceMethod("shape", "a4aM", function(object, value){
	if(all.equal(is(value), is(object@shape))) object@shape <- value
	object
})

#' @rdname a4aM-class
#' @aliases level level-methods level,a4aM-method
setGeneric("level", function(object, ...) standardGeneric("level"))
setMethod("level", "a4aM", function(object) object@level)

#' @rdname a4aM-class
#' @aliases level<- level<--methods level<-,a4aM-method
setGeneric("level<-", function(object,value) standardGeneric("level<-"))
setReplaceMethod("level", "a4aM", function(object, value){
	if(all.equal(is(value), is(object@level))) object@level <- value
	object
})

#' @rdname a4aM-class
#' @aliases trend trend-methods trend,a4aM-method
setGeneric("trend", function(object, ...) standardGeneric("trend"))
setMethod("trend", "a4aM", function(object) object@trend)

#' @rdname a4aM-class
#' @aliases trend<- trend<--methods trend<-,a4aM-method
setGeneric("trend<-", function(object,value) standardGeneric("trend<-"))
setReplaceMethod("trend", "a4aM", function(object, value){
	if(all.equal(is(value), is(object@trend))) object@trend <- value
	object
})

#' @name rngmbar
#' @rdname rngmbar
#' @title M bar range extract and replacement
#' @description Methods to extract from \code{a4aM} objects the M bar age range, or replace its value.
#' @param object an \code{a4aM} object
#' @param value a \code{vector} with max and min age range to replace the object info
#' @return a \code{vector} object when extracting or an \code{a4aM} object when replacing
#' @aliases rngmbar rngmbar-methods rngmbar,a4aM-method
#' @aliases rngmbar<- rngmbar<--methods rngmbar<-,a4aM-method
#' @examples
#' mod1 <- FLModelSim(model=~exp(-age-0.5))
#' mod2 <- FLModelSim(model=~1.5*k, params=FLPar(k=0.4))
#' m1 <- a4aM(shape=mod1, level=mod2)
#' rngmbar(m1)<-c(1,5)
#' rngmbar(m1)

setGeneric("rngmbar", function(object, ...) standardGeneric("rngmbar"))
setMethod("rngmbar", "a4aM", function(object){
	object@range[c("minmbar","maxmbar")]
})

setGeneric("rngmbar<-", function(object,value) standardGeneric("rngmbar<-"))
setReplaceMethod("rngmbar", "a4aM", function(object, value){
	object@range[c("minmbar","maxmbar")] <- value
	object
})

#' @name vecmbar
#' @rdname vecmbar
#' @title M age vector
#' @description Method to extract from \code{a4aM} objects the vector of ages for the mean M.
#' @param object a \code{a4aM} object
#' @return an \code{vector} object
#' @aliases vecmbar vecmbar-methods vecmbar,a4aM-method
#' @examples
#' mod1 <- FLModelSim(model=~exp(-age-0.5))
#' mod2 <- FLModelSim(model=~1.5*k, params=FLPar(k=0.4))
#' m1 <- a4aM(shape=mod1, level=mod2)
#' rngmbar(m1)<-c(1,5)
#' vecmbar(m1)

setGeneric("vecmbar", function(object, ...) standardGeneric("vecmbar"))
setMethod("vecmbar", "a4aM", function(object){
	rng <- object@range[c("minmbar","maxmbar")]
	rng[1]:rng[2]
})


