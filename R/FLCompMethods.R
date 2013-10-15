#' @title year range 
#' @description method to extract from \code{FLComp} objects the year range.
#' @param object a \code{FLComp} object
#' @return a \code{vector} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("rngyear", function(object, ...) standardGeneric("rngyear"))
setMethod("rngyear", "FLComp", function(object){
	object@range[c("minyear","maxyear")]
})

#' @title year range replacement
#' @description method to replace \code{FLComp} object's year range.
#' @param object a \code{FLComp} object
#' @param value a \code{vector} with max and min year range 
#' @return a \code{a4aM} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("rngyear<-", function(object,value) standardGeneric("rngyear<-"))
setReplaceMethod("rngyear", "FLComp", function(object, value){
	object@range[c("minyear","maxyear")] <- value
	object
})

#' @title age range
#' @description method to extract from \code{FLComp} objects the age range.
#' @param object a \code{FLComp} object
#' @return a \code{vector} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("rngage", function(object, ...) standardGeneric("rngage"))
setMethod("rngage", "FLComp", function(object){
	object@range[c("min","max")]
})

#' @title age range replacement
#' @description method to replace \code{FLComp} object's age range.
#' @param object a \code{FLComp} object
#' @param value a \code{vector} with max and min age range 
#' @return a \code{a4aM} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("rngage<-", function(object,value) standardGeneric("rngage<-"))
setReplaceMethod("rngage", "FLComp", function(object, value){
	object@range[c("min","max")] <- value
	object
})

#' @title year vector
#' @description method to extract from \code{FLComp} objects the vector of years.
#' @param object a \code{a4aM} object
#' @return a \code{vector} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("vecyear", function(object, ...) standardGeneric("vecyear"))
setMethod("vecyear", "FLComp", function(object){
	rng <- object@range[c("minyear","maxyear")]
	rng[1]:rng[2]
})

#' @title age vector
#' @description method to extract from \code{FLComp} objects the vector of ages.
#' @param object a \code{a4aM} object
#' @return a \code{vector} object
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#' @export
setGeneric("vecage", function(object, ...) standardGeneric("vecage"))
setMethod("vecage", "FLComp", function(object){
	rng <- object@range[c("min","max")]
	rng[1]:rng[2]
})
