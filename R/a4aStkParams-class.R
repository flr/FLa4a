#' @title Stock parameters class
#' @docType class
#' @name a4aStkParams
#' @rdname a4aStkParams-class
#' @template ClassDescription
#' @section Slot:
#' \describe{
#'  \item{\code{fMod}}{F submodel \code{formula}}
#'  \item{\code{n1Mod}}{first year N \code{formula}}
#'  \item{\code{srMod}}{stock-recruitment submodel \code{formula}}
#'  \item{\code{params}}{\code{FLPar} with parameters}
#'  \item{\code{vcov}}{\code{array} with variance-covariance}
#'  \item{\code{centering}}{centering values \code{numeric}}
#'  \item{\code{distr}}{statistical distribution \code{character}}
#'  \item{\code{m}}{natural mortality \code{FLQuant}}
#'  \item{\code{units}}{data units \code{character}}
#' }
#' @aliases a4aStkParams-class

setClass("a4aStkParams",
  contains = "FLComp",
  slots =
    c(
      fMod         = "formula",
      n1Mod        = "formula",
      srMod        = "formula",
      coefficients = "FLPar",
      vcov         = "array",
      centering    = "FLPar",
      distr        = "character",
      m            = "FLQuant",
      wt           = "FLQuant",
      mat          = "FLQuant",
      units        = "character",
      link         = "function",
      linkinv      = "function"
    )
)

setValidity("a4aStkParams",
  function(object) {
    # check dimensions of m and wt are the same
    if (!identical(unlist(dimnames(object@m)[2:5]),
                   unlist(dimnames(object@wt)[2:5]))) {
      "m and wt elements must share dimensions 2 to 5"
      # can also check against range slot here...
    } else {
      TRUE
    }
})

setMethod("initialize", "a4aStkParams",
    function(.Object,
              fMod         = ~ 1,
              n1Mod        = ~ 1,
              srMod        = ~ 1,
              coefficients = FLPar(),
              vcov         = array(),
              centering    = FLPar(centering = 0),
              distr        = "norm",
              m            = FLQuant(),
              wt           = FLQuant(),
              mat          = FLQuant(),
              units        = "NA",
              ink          = log,
              linkinv      = exp,
              ...) {
      # initialize FLComp slots
      .Object <- callNextMethod(.Object, ...)
      # initialize remaining slots
      .Object@fMod  <- fMod
      .Object@n1Mod <- n1Mod
      .Object@srMod <- srMod
      .Object@coefficients <- coefficients
      .Object@vcov <- vcov
      if (length(dim(vcov)) == 2) dim(.Object@vcov) <- c(dim(vcov), 1)
      .Object@centering <- centering
      .Object@distr <- distr
      # if missing set dimensions of of m and wt based on range
      if (missing(m) || missing(wt))
        flq <- flq_from_range(.Object)
      .Object@m <- if (missing(m)) flq else m
      .Object@wt <- if (missing(wt)) flq else wt
      .Object@wt <- if (missing(mat)) flq else mat
      # throw error if range from FLComp doesn't match FLQuants
      # (can't check this in setValidity due to callNextMethod resulting in an
      # invalid a4aStkParams object when range is supplied)
      if (abs(as.numeric(dimnames(.Object@m)$year[1]) -
              .Object@range["minyear"]) > 1e-9 ||
          abs(as.numeric(dimnames(.Object@m)$year[dim(.Object@m)[2]]) -
              .Object@range["maxyear"]) > 1e-9) {
            stop("range does not match supplied m and wt dimensions")
      }
      .Object@units <- units
      .Object
    })






#' @rdname a4aStkParams-class
#' @template bothargs
#' @aliases a4aStkParams a4aStkParams-methods
#' @template Accessors
#' @template Constructors
setGeneric("a4aStkParams",
  function(object, ...) standardGeneric("a4aStkParams"))

#' @rdname a4aStkParams-class
setMethod("a4aStkParams", signature(object = "missing"),
  function(...) {
    # empty
    if(missing(...)){
        new("a4aStkParams")
    # or not
    } else {
      args <- list(...)
      args$Class <- "a4aStkParams"
      do.call("new", args)
      }
  }
)


#' @rdname a4aStkParams-class
setMethod("m", signature(object = "a4aStkParams"),
  function(object) object@m)

#' @rdname a4aStkParams-class
setMethod("wt", signature(object = "a4aStkParams"),
  function(object) object@wt)

#' @rdname a4aStkParams-class
setMethod("mat", signature(object = "a4aStkParams"),
  function(object) object@mat)

#' @rdname a4aStkParams-class
#' @aliases fMod fMod-methods
setGeneric("fMod", function(object, ...) standardGeneric("fMod"))
#' @rdname a4aStkParams-class
setMethod("fMod", "a4aStkParams", function(object) object@fMod)

#' @rdname a4aStkParams-class
#' @aliases fmodel fmodel-methods
setGeneric("fmodel", function(object, ...) standardGeneric("fmodel"))
#' @rdname a4aStkParams-class
setMethod("fmodel", "a4aStkParams",
  function(object) {
    stk_submodel <- as(object, "submodels")
    stk_submodel$fmodel
  }
)


#' @rdname a4aStkParams-class
#' @param value the new object
#' @aliases fMod<- fMod<--methods
setGeneric("fMod<-", function(object, value) standardGeneric("fMod<-"))
#' @rdname a4aStkParams-class
setReplaceMethod("fMod", signature("a4aStkParams", "formula"),
  function(object, value) {
    if (all.equal(is(value), is(object@fMod)))
      object@fMod <- value
    object
  }
)

#' @rdname a4aStkParams-class
#' @aliases n1Mod n1Mod-methods
setGeneric("n1Mod", function(object, ...) standardGeneric("n1Mod"))
#' @rdname a4aStkParams-class
setMethod("n1Mod", "a4aStkParams", function(object) object@n1Mod)

#' @rdname a4aStkParams-class
#' @aliases n1Mod<- n1Mod<--methods
setGeneric("n1Mod<-", function(object, value) standardGeneric("n1Mod<-"))
#' @rdname a4aStkParams-class
setReplaceMethod("n1Mod", signature("a4aStkParams", "formula"),
  function(object, value) {
    if (all.equal(is(value), is(object@n1Mod)))
      object@n1Mod <- value
    object
  }
)

#' @rdname a4aStkParams-class
#' @aliases srMod rMod-methods
setGeneric("srMod", function(object, ...) standardGeneric("srMod"))
#' @rdname a4aStkParams-class
setMethod("srMod", "a4aStkParams", function(object) object@srMod)

#' @rdname a4aStkParams-class
#' @aliases srMod<- srMod<--methods
setGeneric("srMod<-", function(object, value) standardGeneric("srMod<-"))
#' @rdname a4aStkParams-class
setReplaceMethod("srMod", signature("a4aStkParams", "formula"),
  function(object, value){
    if (all.equal(is(value), is(object@srMod)))
      object@srMod <- value
    object
  })

#' @rdname a4aStkParams-class
setMethod("params", "a4aStkParams", function(object) object@coefficients)

#' @rdname a4aStkParams-class
setReplaceMethod("params", signature("a4aStkParams", "FLPar"),
  function(object, value) {
    if (all.equal(is(value), is(params(object))))
      object@coefficients <- value
    object
  }
)

#' @rdname a4aStkParams-class
#' @aliases coefficients coefficients-methods
setGeneric("coefficients",
  function(object, ...) standardGeneric("coefficients"))
#' @rdname a4aStkParams-class
setMethod("coefficients", "a4aStkParams", function(object) object@coefficients)

#' @rdname a4aStkParams-class
#' @aliases coefficients<- coefficients<--methods
setGeneric("coefficients<-",
  function(object, value) standardGeneric("coefficients<-"))
#' @rdname a4aStkParams-class
setReplaceMethod("coefficients", signature("a4aStkParams", "FLPar"),
  function(object, value){
    if (all.equal(is(value), is(coefficients(object))))
      object@coefficients <- value
    object
  }
)

#' @rdname a4aStkParams-class
setMethod("distr", "a4aStkParams", function(object) object@distr)

#' @rdname a4aStkParams-class
setReplaceMethod("distr", signature("a4aStkParams", "character"),
  function(object, value){
    if (all.equal(is(value), is(object@distr)))
      object@distr <- value
    object
  }
)

#' @rdname a4aStkParams-class
setMethod("vcov", "a4aStkParams", function(object) object@vcov)

#' @rdname a4aStkParams-class
setReplaceMethod("vcov", signature("a4aStkParams", "array"),
  function(object, value){
    if(all.equal(is(value), is(object@vcov)))
      object@vcov <- value
    object
  }
)



#
#  show methods
#

setMethod("show", "a4aStkParams",
  function(object) {
    cat("stkmodel:\n")
    if (length(object) == 0) {
      cat("empty object\n")
    } else {
      fmt <- paste0("\t %2smodel: ")
      cat(sprintf(fmt, "f")); print(fMod(object), showEnv = FALSE)
      cat(sprintf(fmt, "n1")); print(n1Mod(object), showEnv = FALSE)
      cat(sprintf(fmt, "sr")); print(srMod(object), showEnv = FALSE)
    }
 })





#
# Coersion methods
#

# method.skeleton("coerce", "a4aStkParams",  file = stdout())

setMethod("coerce", signature(from = "a4aStkParams", to = "FLSR"),
  function (from, to, strict = TRUE) {
    # get SR model formula
    srmodel <- geta4aSRmodel(srMod(from))
    # get FLSR definition
    expr_model <- a4aSRmodelDefinitions(srmodel)

    # build skeleton FLSR
    flsr <-
      FLSR(
        formula(
          paste("rec ~ (",
                deparse(expr_model, width.cutoff = 500),
                ") * ",
                exp(from@centering))
          ))

    # get SR pars
    cnames <- rownames(coef(from))
    params(flsr) <- FLPar(a = exp(coef(from)[grep("sraMod", cnames)]),
                          b = exp(coef(from)[grep("srbMod", cnames)]))

    which <- c(grep("sraMod", cnames), grep("srbMod", cnames))
    vcov(flsr) <- vcov(from)[which, which,, drop = FALSE]
    dimnames(vcov(flsr)) <- list(c("a", "b"), c("a", "b"))

    flqs <- genFLQuant(from)

    rec(flsr) <- flqs$stock.n[1, ]
    ssb <- quantSums(flqs$stock.n * mat(from) * wt(from))
    ssb(flsr) <- ssb
    fitted(flsr) <-
      ssb / ssb * eval(expr_model, c(as(params(flsr), "list"), ssb = ssb))
    residuals(flsr) <- log(rec(flsr)) - log(fitted(flsr))

    range(flsr) <- range(from)[c("min", "max", "minyear", "maxyear")]
    name(flsr) <- name(from)
    desc(flsr) <- desc(from)

    flsr
  }
)

# method.skeleton("coerce", "a4aStkParams",  file = stdout())

setMethod("coerce", signature(from = "a4aStkParams", to = "submodels"),
  function (from, to, strict = TRUE) {

    # extract come bits from stkmodel
    b <- coef(from)
    vmat <- vcov(from)
    f_which <- which(grepl("fMod", rownames(b)))
    r_which <- which(grepl("rMod", rownames(b)))
    n1_which <- which(grepl("n1Mod", rownames(b)))

    fmodel <-
      submodel(formula = fMod(from),
               coefficients = b[f_which],
               vcov = vmat[f_which, f_which,, drop = FALSE],
               centering = from@centering,
               distr = from@distr,
               link = from@link,
               linkinv = from@link,
               range = range(from)
           )
    n1model <-
      submodel(formula = n1Mod(from),
               coefficients = b[n1_which],
               vcov = vmat[n1_which, n1_which,, drop = FALSE],
               centering = from@centering,
               distr = from@distr,
               link = from@link,
               linkinv = from@link,
               range = range(from)
           )
    rmodel <-
      submodel(formula = srMod(from),
               coefficients = b[r_which],
               vcov = vmat[r_which, r_which,, drop = FALSE],
               centering = from@centering,
               distr = from@distr,
               link = from@link,
               linkinv = from@link,
               range = range(from)
           )


    # construct submodels
    stk_submodel <-
      submodels(list(fmodel, n1model, rmodel),
                names = c("fmodel", "n1model", "rmodel"))

    # calculate correlation matrix for each iter
    corrmat <- vmat
    for (i in 1:dim(vmat)[3]) {
      corrmat[,, i] <- cov2cor(vmat[,, i])
    }
    # f and r correlation
    stk_submodel@corBlocks$fmodel.n1model[] <-
      corrmat[f_which, n1_which,, drop = FALSE]
    # f and n1 model correlation
    stk_submodel@corBlocks$fmodel.rmodel[] <-
      corrmat[f_which, r_which,, drop = FALSE]
    # r and n1 model correlation
    stk_submodel@corBlocks$n1model.rmodel[] <-
      corrmat[n1_which, r_which,, drop = FALSE]

    stk_submodel
  }
)




#' @rdname a4aStkParams-class
#' @param iter the number of iterations to create
#' @param fill.iter should the new iterations be filled with values (TRUE) or NAs (FALSE)
setMethod("propagate", signature(object = "a4aStkParams"),
  function(object, iter, fill.iter = TRUE) {

    # propagate coefs, centering, m and wt
    object@coefficients <-
      propagate(object@coefficients, iter, fill.iter = fill.iter)
    object@centering <- propagate(object@centering, iter, fill.iter = fill.iter)
    object@m <- propagate(object@m, iter, fill.iter = fill.iter)
    object@wt <- propagate(object@wt, iter, fill.iter = fill.iter)
    object@mat <- propagate(object@mat, iter, fill.iter = fill.iter)

    # now propagate vcov
    vcov.iter <- vcov(object)
    dob <- dim(vcov.iter)

    if (iter != dob[3]) {
      # CHECK no iters in object
      if (dob[3] > 1) stop("propagate can only extend objects with no iters")

      object@vcov <-
        array(NA, dim = c(dob[1:2], iter),
              dimnames = c(dimnames(vcov.iter)[1:2], list(1:iter)))
      if (fill.iter) {
        object@vcov[] <- as.vector(vcov.iter)
      } else {
        object@vcov[,, 1] <- as.vector(vcov.iter)
      }
    }

    object
  }
)


#' @rdname a4aStkParams-class
#' @param obj the object to be subset
#' @param it iteration to be extracted
setMethod("iter", "a4aStkParams",
  function(obj, it) {
    obj@vcov <- obj@vcov[,, it, drop = FALSE]
    obj@coefficients <- iter(obj@coefficients, it)
    obj@m <- iter(obj@m, it)
    obj@wt <- iter(obj@wt, it)
    obj@mat <- iter(obj@mat, it)
    obj@centering <- iter(obj@centering, it)
    obj
  }
)
