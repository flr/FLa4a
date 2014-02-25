#' @name rngyear
#' @rdname rngyear
#' @title year range extract and replacement
#' @description Methods to extract from \code{FLComp} objects the year range or replace its value.
#' @param object a \code{FLComp} object
#' @param value a \code{vector} with max and min year range to replace the object info. 
#' @return a \code{vector} object when extracting or a \code{FLComp} object when replacing
#' @aliases rngyear rngyear-methods rngyear,FLComp-method rngyear<- rngyear<--methods rngyear<-,FLComp,numeric-method
setGeneric("rngyear", function(object, ...) standardGeneric("rngyear"))
setMethod("rngyear", "FLComp", function(object){
	object@range[c("minyear","maxyear")]
})

setGeneric("rngyear<-", function(object,value) standardGeneric("rngyear<-"))
setReplaceMethod("rngyear", signature("FLComp","numeric"), function(object, value){
	object@range[c("minyear","maxyear")] <- value
	object
})

#' @name rngage
#' @rdname rngage
#' @title age range extract and replacement
#' @description Methods to extract from \code{FLComp} objects the age range or replace its value.
#' @param object a \code{FLComp} object
#' @param value a \code{vector} with max and min age range to replace the object info. 
#' @return a \code{vector} object when extracting or a \code{FLComp} object when replacing
#' @aliases rngage rngage-methods rngage,FLComp-method rngage<- rngage<--methods rngage<-,FLComp,numeric-method
setGeneric("rngage", function(object, ...) standardGeneric("rngage"))
setMethod("rngage", "FLComp", function(object){
	object@range[c("min","max")]
})

setGeneric("rngage<-", function(object,value) standardGeneric("rngage<-"))
setReplaceMethod("rngage", signature("FLComp","numeric"), function(object, value){
	object@range[c("min","max")] <- value
	object
})

#' @name vecyear
#' @rdname vecyear
#' @title year vector extract
#' @description method to extract from \code{FLComp} objects the vector of years.
#' @param object a \code{FLComp} object
#' @return a \code{vector} object
#' @aliases vecyear vecyear-methods vecyear,FLComp-method
setGeneric("vecyear", function(object, ...) standardGeneric("vecyear"))
setMethod("vecyear", "FLComp", function(object){
	rng <- object@range[c("minyear","maxyear")]
	rng[1]:rng[2]
})

#' @name vecage
#' @rdname vecage
#' @title age vector extract
#' @description method to extract from \code{FLComp} objects the vector of ages
#' @param object a \code{FLComp} object
#' @return a \code{vector} object
#' @aliases vecage vecage-methods vecage,FLComp-method
setGeneric("vecage", function(object, ...) standardGeneric("vecage"))
setMethod("vecage", "FLComp", function(object){
	rng <- object@range[c("min","max")]
	rng[1]:rng[2]
})

