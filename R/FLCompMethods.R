##' @name rngyear
##' @rdname rngyear
##' @title year range extract and replacement
##' @description Methods to extract from \code{FLComp} objects the year range, or to replace its value.
##' @param object an \code{FLComp} object
##' @param value a \code{vector} with max and min year range to replace the object info
##' @return a \code{vector} object when extracting, or an \code{FLComp} object when replacing
##' @template dots
##' @aliases rngyear rngyear-methods
##' @examples
##' data(ple4)
##' data(ple4.index)
##' rngyear(ple4)
##' rngyear(ple4.index)
#setGeneric("rngyear", function(object, ...) standardGeneric("rngyear"))
##' @rdname rngyear
#setMethod("rngyear", "FLComp", function(object){
#	object@range[c("minyear","maxyear")]
#})

##' @rdname rngyear
##' @aliases rngyear<- rngyear<--methods
#setGeneric("rngyear<-", function(object,value) standardGeneric("rngyear<-"))
##' @rdname rngyear
#setReplaceMethod("rngyear", signature("FLComp","numeric"), function(object, value){
#	object@range[c("minyear","maxyear")] <- value
#	object
#})

##' @name rngage
##' @rdname rngage
##' @title age range extract and replacement
##' @description Methods to extract from \code{FLComp} objects the age range, or to replace its value.
##' @template bothargs
##' @param value a \code{vector} with max and min age range to replace the object info
##' @return a \code{vector} object when extracting, or an \code{FLComp} object when replacing
##' @aliases rngage rngage-methods
##' @examples
##' data(ple4)
##' data(ple4.index)
##' rngage(ple4)
##' rngage(ple4.index)
#setGeneric("rngage", function(object, ...) standardGeneric("rngage"))
##' @rdname rngage
#setMethod("rngage", "FLComp", function(object){
#	object@range[c("min","max")]
#})

##' @rdname rngage
##' @aliases rngage<- rngage<--methods
#setGeneric("rngage<-", function(object,value) standardGeneric("rngage<-"))
##' @rdname rngage
#setReplaceMethod("rngage", signature("FLComp","numeric"), function(object, value){
#	object@range[c("min","max")] <- sort(value)
#	object
#})

##' @name rnglen
##' @rdname rnglen
##' @title length range extract and replacement
##' @description Methods to extract from \code{FLComp} objects the length range, or to replace its value.
##' @template bothargs
##' @param value a \code{vector} with max and min length range to replace the object info
##' @return a \code{vector} object when extracting, or an \code{FLComp} object when replacing
##' @aliases rnglen rnglen-methods
##' @examples
##' data(rfLen)
##' rnglen(rfTrawl.idx)
#setGeneric("rnglen", function(object, ...) standardGeneric("rnglen"))
##' @rdname rnglen
#setMethod("rnglen", "FLComp", function(object){
#	object@range[c("min","max")]
#})

##' @rdname rnglen
##' @aliases rnglen<- rnglen<--methods
#setGeneric("rnglen<-", function(object,value) standardGeneric("rnglen<-"))
##' @rdname rnglen
#setReplaceMethod("rnglen", signature("FLComp","numeric"), function(object, value){
#	object@range[c("min","max")] <- sort(value)
#	object
#})

##' @name rngquant
##' @rdname rngquant
##' @title quant range extract and replacement
##' @description Methods to extract from \code{FLComp} objects the quant range, or to replace its value.
##' @template bothargs
##' @param value a \code{vector} with max and min quant range to replace the object info
##' @return a \code{vector} object when extracting, or an \code{FLComp} object when replacing
##' @aliases rngquant rngquant-methods
##' @examples
##' data(ple4.index)
##' data(rfLen)
##' rngquant(ple4.index)
##' rngquant(rfTrawl.idx)
#setGeneric("rngquant", function(object, ...) standardGeneric("rngquant"))
##' @rdname rngquant
#setMethod("rngquant", "FLComp", function(object){
#	object@range[c("min","max")]
#})

##' @rdname rngquant
##' @aliases rngquant<- rngquant<--methods
#setGeneric("rngquant<-", function(object,value) standardGeneric("rngquant<-"))
##' @rdname rngquant
#setReplaceMethod("rngquant", signature("FLComp","numeric"), function(object, value){
#	object@range[c("min","max")] <- sort(value)
#	object
#})

##' @name vecyear
##' @rdname vecyear
##' @title year vector extract
##' @description Method to extract from \code{FLComp} objects the vector of years.
##' @template bothargs
##' @return a \code{vector} object
##' @aliases vecyear vecyear-methods
##' @examples
##' data(ple4)
##' data(ple4.index)
##' vecyear(ple4)
##' vecyear(ple4.index)
#setGeneric("vecyear", function(object, ...) standardGeneric("vecyear"))
##' @rdname vecyear
#setMethod("vecyear", "FLComp", function(object){
#	rng <- object@range[c("minyear","maxyear")]
#	rng[1]:rng[2]
#})

##' @name vecage
##' @rdname vecage
##' @title age vector extract
##' @description Method to extract from \code{FLComp} objects the vector of ages.
##' @template bothargs
##' @return a \code{vector} object
##' @aliases vecage vecage-methods
##' @examples
##' data(ple4)
##' data(ple4.index)
##' vecage(ple4)
##' vecage(ple4.index)
#setGeneric("vecage", function(object, ...) standardGeneric("vecage"))
##' @rdname vecage
#setMethod("vecage", "FLComp", function(object){
#	rng <- object@range[c("min","max")]
#	rng[1]:rng[2]
#})

##' @name veclen
##' @rdname veclen
##' @title length vector
##' @description Method to extract from \code{FLComp} objects the vector of lengths.
##' @template bothargs
##' @return a \code{vector} object
##' @aliases veclen veclen-methods
##' @examples
##' data(rfLen)
##' veclen(rfTrawl.idx)
#setGeneric("veclen", function(object, ...) standardGeneric("veclen"))
##' @rdname veclen
#setMethod("veclen", "FLComp", function(object){
#	rng <- object@range[c("min","max")]
#	rng[1]:rng[2]
#})

##' @name vecquant
##' @rdname vecquant
##' @title length vector
##' @description method to extract from \code{FLComp} objects the vector of lengths.
##' @template bothargs
##' @return a \code{vector} object
##' @aliases vecquant vecquant-methods
##' @examples
##' data(ple4.index)
##' data(rfLen)
##' vecquant(ple4.index)
##' vecquant(rfTrawl.idx)
#setGeneric("vecquant", function(object, ...) standardGeneric("vecquant"))
##' @rdname vecquant
#setMethod("vecquant", "FLComp", function(object){
#	rng <- object@range[c("min","max")]
#	rng[1]:rng[2]
#})

