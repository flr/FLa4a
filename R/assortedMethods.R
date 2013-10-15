#' Assorted methods needed by FLa4a
#' @import FLCore
#' @include FLModelSimMethods.R

#' Check if object is empty
#'
#' @param object an object
#' @return \code{TRUE} if object is of length 0 
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
#' @examples
#' is.empty(list())
#' is.empty(list(a=2))

is.empty <- function(object) {
	length(object) == 0
}

#' Checks that the name of the second dimension in params is "iter". 
#' For internal use, not very interesting for users 
#' @param object a \code{FLPar} object
#' @return logical
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
#' @examples
#' pars2dim(FLPar())
#' pars2dim(FLPar(array(dim=c(1,1,1))))
setMethod("pars2dim", "FLPar", function(object) {

	dnm <- dimnames(object)
	names(dnm)[2]=="iter"

})

#' Gets the FLQuant's numeric id for a vector of "years".
#' For internal use, not very interesting for users 
#' @param object a \code{FLQuant} object
#' @param year a \code{vector} of years
#' @return \code{numeric vector} that can be used to subset the \code{FLQuant}
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
#' @examples
#' data(ple4)
#' flq <- catch(ple4)
#' getYidx(flq, 2000:2004)
#' flq[, getYidx(flq, 2000:2004)]

setGeneric("getYidx", function(object, ...) standardGeneric("getYidx"))
setMethod("getYidx", "FLQuant", function(object, year) {
	yrs <- dimnames(object)[[2]]
	if(sum(year>0)>0){
		idx <- match(as.character(year), yrs, nomatch=0)
		if(sum(idx)>0){
			idx
		} else {
			year
		}
	} else {
		length(yrs)+year+1
	} 

})


setMethod("iterMedians", "FLPar", function(x){
	apply(x, -match("iter", names(dimnames(x))), median)
})

#' Compute number of iterations
#'
#' @param object object a FL* class object
#' @return numeric
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("niters", function(object, ...) standardGeneric("niters"))
setMethod("niters", "FLModelSim", function(object){
	dim(params(object))[2]
})

##' Generates accessores and replacements for classes
##'
##' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
##' @param object a class object
##' @return generates the accessors and replacement methods
##' @export

#setGeneric("genAR", function(object, ...) standardGeneric("genAR"))
#setMethod("genAR", "ANY", function(object, ...) {
#	slts <- slotNames(object)

#	for(i in slts){
#		if(!isGeneric(i)){
#			do.call("setGeneric", list(name=i, def=eval(parse(text=paste("function(object, ...) standardGeneric(\"", i, "\")", sep="")))))
#		}
#		do.call("setMethod", list(f=i, signature=is(object)[1], definition=eval(parse(text=paste("function(object) object@", i, sep="")))))
#		
#		if(!isGeneric(paste(i, "<-", sep=""))){
#			do.call("setGeneric", list(name=paste(i, "<-", sep=""), def=eval(parse(text=paste("function(object, value) standardGeneric(\"", paste(i, "<-", sep=""), "\")", sep="")))))
#		}

#		do.call("setReplaceMethod", list(f=i, signature=is(object)[1], def=eval(parse(text=paste("function(object, value){if(all.equal(is(value), is(object@", i, "))) object@", i, "<- value; object}", sep="")))))
#	}

#})


