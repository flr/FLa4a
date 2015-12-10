#' @name rngyear
#' @rdname rngyear
#' @title year range extract and replacement
#' @description Methods to extract from \code{FLComp} objects the year range, or to replace its value.
#' @param object an \code{FLComp} object
#' @param value a \code{vector} with max and min year range to replace the object info
#' @return a \code{vector} object when extracting, or an \code{FLComp} object when replacing
#' @aliases rngyear rngyear-methods rngyear,FLComp-method rngyear<- rngyear<--methods rngyear<-,FLComp,numeric-method
#' @examples
#' data(ple4)
#' data(ple4.index)
#' rngyear(ple4)
#' rngyear(ple4.index)

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
#' @description Methods to extract from \code{FLComp} objects the age range, or to replace its value.
#' @param object an \code{FLComp} object
#' @param value a \code{vector} with max and min age range to replace the object info
#' @return a \code{vector} object when extracting, or an \code{FLComp} object when replacing
#' @aliases rngage rngage-methods rngage,FLComp-method rngage<- rngage<--methods rngage<-,FLComp,numeric-method
#' @examples
#' data(ple4)
#' data(ple4.index)
#' rngage(ple4)
#' rngage(ple4.index)

setGeneric("rngage", function(object, ...) standardGeneric("rngage"))
setMethod("rngage", "FLComp", function(object){
	object@range[c("min","max")]
})

setGeneric("rngage<-", function(object,value) standardGeneric("rngage<-"))
setReplaceMethod("rngage", signature("FLComp","numeric"), function(object, value){
	object@range[c("min","max")] <- sort(value)
	object
})

#' @name rnglen
#' @rdname rnglen
#' @title length range extract and replacement
#' @description Methods to extract from \code{FLComp} objects the length range, or to replace its value.
#' @param object an \code{FLComp} object
#' @param value a \code{vector} with max and min length range to replace the object info
#' @return a \code{vector} object when extracting, or an \code{FLComp} object when replacing
#' @aliases rnglen rnglen-methods rnglen,FLComp-method rnglen<- rnglen<--methods rnglen<-,FLComp,numeric-method
#' @examples
#' data(rfLen)
#' rnglen(rfTrawl.idx)

setGeneric("rnglen", function(object, ...) standardGeneric("rnglen"))
setMethod("rnglen", "FLComp", function(object){
	object@range[c("min","max")]
})

setGeneric("rnglen<-", function(object,value) standardGeneric("rnglen<-"))
setReplaceMethod("rnglen", signature("FLComp","numeric"), function(object, value){
	object@range[c("min","max")] <- sort(value)
	object
})

#' @name rngquant
#' @rdname rngquant
#' @title quant range extract and replacement
#' @description Methods to extract from \code{FLComp} objects the quant range, or to replace its value.
#' @param object an \code{FLComp} object
#' @param value a \code{vector} with max and min quant range to replace the object info
#' @return a \code{vector} object when extracting, or an \code{FLComp} object when replacing
#' @aliases rngquant rngquant-methods rngquant,FLComp-method rngquant<- rngquant<--methods rngquant<-,FLComp,numeric-method
#' @examples
#' data(ple4.index)
#' data(rfLen)
#' rngquant(ple4.index)
#' rngquant(rfTrawl.idx)

setGeneric("rngquant", function(object, ...) standardGeneric("rngquant"))
setMethod("rngquant", "FLComp", function(object){
	object@range[c("min","max")]
})

setGeneric("rngquant<-", function(object,value) standardGeneric("rngquant<-"))
setReplaceMethod("rngquant", signature("FLComp","numeric"), function(object, value){
	object@range[c("min","max")] <- sort(value)
	object
})

#' @name vecyear
#' @rdname vecyear
#' @title year vector extract
#' @description Method to extract from \code{FLComp} objects the vector of years.
#' @param object an \code{FLComp} object
#' @return a \code{vector} object
#' @aliases vecyear vecyear-methods vecyear,FLComp-method
#' @examples
#' data(ple4)
#' data(ple4.index)
#' vecyear(ple4)
#' vecyear(ple4.index)

setGeneric("vecyear", function(object, ...) standardGeneric("vecyear"))
setMethod("vecyear", "FLComp", function(object){
	rng <- object@range[c("minyear","maxyear")]
	rng[1]:rng[2]
})

#' @name vecage
#' @rdname vecage
#' @title age vector extract
#' @description Method to extract from \code{FLComp} objects the vector of ages.
#' @param object an \code{FLComp} object
#' @return a \code{vector} object
#' @aliases vecage vecage-methods vecage,FLComp-method
#' @examples
#' data(ple4)
#' data(ple4.index)
#' vecage(ple4)
#' vecage(ple4.index)

setGeneric("vecage", function(object, ...) standardGeneric("vecage"))
setMethod("vecage", "FLComp", function(object){
	rng <- object@range[c("min","max")]
	rng[1]:rng[2]
})

#' @name veclen
#' @rdname veclen
#' @title length vector
#' @description Method to extract from \code{FLComp} objects the vector of lengths.
#' @param object an \code{FLComp} object
#' @return a \code{vector} object
#' @aliases veclen veclen-methods veclen,FLComp-method
#' @examples
#' data(rfLen)
#' veclen(rfTrawl.idx)

setGeneric("veclen", function(object, ...) standardGeneric("veclen"))
setMethod("veclen", "FLComp", function(object){
	rng <- object@range[c("min","max")]
	rng[1]:rng[2]
})

#' @name vecquant
#' @rdname vecquant
#' @title length vector
#' @description method to extract from \code{FLComp} objects the vector of lengths.
#' @param object a \code{FLComp} object
#' @return a \code{vector} object
#' @aliases vecquant vecquant-methods vecquant,FLComp-method
#' @examples
#' data(ple4.index)
#' data(rfLen)
#' vecquant(ple4.index)
#' vecquant(rfTrawl.idx)

setGeneric("vecquant", function(object, ...) standardGeneric("vecquant"))
setMethod("vecquant", "FLComp", function(object){
	rng <- object@range[c("min","max")]
	rng[1]:rng[2]
})

