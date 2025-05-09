#' @title S4 class \code{a4aFit}
#'
#' @description The \code{a4aFit} class was built to store the a4a stock assessment fits.
#'
#' @section Slots:
#' \describe{
#'    \item{call}{The function call}
#'    \item{clock}{Information on call duration}
#'    \item{fitSumm}{Fit summary}
#'    \item{stock.n}{Estimates of stock numbers-at-age}
#'    \item{harvest}{Estimates of fishing mortality at age}
#'    \item{catch.n}{Estimates of catch numbers-at-age}
#'    \item{index}{Estimates of survey or CPUE indices-at-age}
#'  }
#' @template Accessors
#' @template Constructors
#' @docType class
#' @name a4aFit-class
#' @rdname a4aFit-class
#' @aliases a4aFit-class
#' @examples
#' data(ple4)
#' data(ple4.index)
#'
#' obj <- sca(stock=ple4, indices=FLIndices(ple4.index))
#' obj
#'
#' slotNames(obj)
#' clock(obj)
#' fitSumm(obj)
#'
#' flq <- stock.n(obj)
#' is(flq)
#' flq <- index(obj)
#' is(flq)
#'
#' logLik(obj)
#' AIC(obj)
#' BIC(obj)

a4aFit <-
  setClass("a4aFit",
           contains = "FLComp",
           slots = c(call    = "call",
                     clock   = "numeric",
                     fitSumm = "array",
                     stock.n = "FLQuant",
                     harvest = "FLQuant",
                     catch.n = "FLQuant",
                     index   = "FLQuants"))

#' @rdname a4aFit-class
#' @template bothargs
#' @aliases a4aFit a4aFit-methods
setGeneric("a4aFit")
 
setMethod("initialize", "a4aFit",
    function(.Object,
             ...,
             call, clock, fitSumm,
             stock.n, harvest, catch.n, index) {
      if (!missing(call)) .Object@call <- call
      if (!missing(clock)) .Object@clock   <- clock
      if (!missing(fitSumm)) .Object@fitSumm <- fitSumm
      if (!missing(stock.n)) .Object@stock.n <- stock.n
      if (!missing(harvest)) .Object@harvest <- harvest
      if (!missing(catch.n)) .Object@catch.n <- catch.n
      if (!missing(index)) .Object@index   <- index
      .Object <- callNextMethod(.Object, ...)
      .Object
    })

setValidity("a4aFit",
  function(object) {
    # All FLQuant objects must have same dimensions
    if (any(dim(object@harvest) != dim(object@stock.n)) ||
        any(dim(object@catch.n) != dim(object@stock.n)))
      "stock.n, catch.n and harvest slots must have same dimensions"
    else # Everything is fine
      TRUE
})


#
#  accessor methods
#

#' @rdname a4aFit-class
#' @aliases clock clock-methods
setGeneric("clock", function(object, ...) standardGeneric("clock"))
#' @rdname a4aFit-class
setMethod("clock", "a4aFit", function(object) object@clock)

#' @rdname a4aFit-class
#' @aliases fitSumm fitSumm-methods
setGeneric("fitSumm", function(object, ...) standardGeneric("fitSumm"))
#' @rdname a4aFit-class
setMethod("fitSumm", "a4aFit", function(object) object@fitSumm)

#' @rdname a4aFit-class
setMethod("stock.n", "a4aFit", function(object) object@stock.n)

#' @rdname a4aFit-class
setMethod("harvest", "a4aFit", function(object) object@harvest)

#' @rdname a4aFit-class
setMethod("catch.n", "a4aFit", function(object) object@catch.n)

#' @rdname a4aFit-class
setMethod("index", "a4aFit", function(object) object@index)


#' @rdname a4aFit-class
setMethod("show", signature(object = "a4aFit"),
  function(object)
  {
    cat("a4a model fit for:", object@name, "\n")
    cat("\nCall:\n", paste(deparse(object @ call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")

    cat("Time used:\n")
    print(object @ clock)

 })


#
# Other methods
#

#' @rdname a4aFit-class
setMethod("logLik", signature(object = "a4aFit"),
  function(object, ...)
  {
    dim2 <- length(dim(object @ fitSumm))
    if (dim2 == 1) {
      val <- -1 * unname(object @ fitSumm["nlogl"])
      attr(val, "nobs") <- unname(object @ fitSumm["nobs"])
      attr(val, "df") <- unname(object@fitSumm["nopar"])
    } else if (dim2 == 2) {
      val <- -1 * unname(object @ fitSumm["nlogl",])
      attr(val, "nobs") <- unname(object @ fitSumm["nobs",])
      attr(val, "df") <- unname(object@fitSumm["nopar",])
    }
    class(val) <- "logLik"
    val
 })

#' @rdname a4aFit-class
#' @param obj the object to be subset
#' @param it iteration to be extracted
setMethod("iter", "a4aFit", function(obj, it){
	if(dim(obj@fitSumm)[2]>1) obj@fitSumm <- obj@fitSumm[,it, drop=FALSE]
	obj@harvest <- iter(obj@harvest, it)
	obj@stock.n <- iter(obj@stock.n, it)
	obj@catch.n <- iter(obj@catch.n, it)
	obj@index <- iter(obj@index, it)
	obj
})


setMethod("coerce", signature(from = "a4aFit", to = "FLSR"),
  function (from, to, strict = TRUE)
  {
    res <- as(stkmodel(pars(from)), "FLSR")

    return(res)
  }
)

#' @rdname a4aFit-class
#' @param x the object to be subset
#' @param start initial year
#' @param end final year
#' @param extend if object is shorter than end-start extend to cover year range
#' @param frequency interval between years if extended
setMethod("window", signature(x="a4aFit"),
  function(x, start=dims(x)$minyear, end=dims(x)$maxyear, extend=TRUE, frequency=1)
  {
    x <- qapply(x, window, start=start, end=end, extend=extend, frequency=frequency)
    x@index <- window(x@index, start=start, end=end, extend=extend, frequency=frequency)
    x@range["minyear"] <- as.numeric(dims(x)$minyear)
    x@range["maxyear"] <- as.numeric(dims(x)$maxyear)
    return(x)
  }
)



