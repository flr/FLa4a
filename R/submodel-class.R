#' @title Submodel class
#' @docType class
#' @name submodel
#' @rdname submodel-class
#' @template ClassDescription
#' @section Slot:
#' \describe{
#'
#'  \item{\code{Mod}}{\code{formula} describing the model}
#'
#'  \item{\code{params}}{\code{FLPar} with model parameters}
#'
#'  \item{\code{vcov}}{\code{array} with variance covariance paramaters
#'                                  related to the variance model}
#'
#'  \item{\code{centering}}{\code{numeric} value used for centering the data}
#'
#'  \item{\code{distr}}{a character with the parameters' statistical
#'                      distribution; it must match a known distribution for
#'                      R (\emph{e.g.} "norm" for gaussian) so that \code{rnorm}
#'                      can be called}
#' }
#' @aliases submodel-class
setClass("submodel",
  contains = "FLComp",
  slots = c(formula      = "formula",
            coefficients = "FLPar",
            vcov         = "array",
            centering    = "FLPar",
            distr        = "character",
            link         = "function",
            linkinv      = "function")
)

setValidity("submodel",
  function(object) {
    # no unit, area, season fits
    if (length(dim(coef(object))) > 2 | length(dim(object@vcov)) > 3) {
      "Params or vcov have unit, area or season. Can't work with that!"
    } else
    if (FALSE) {
      "coefficients do not match formula"
    } else {
      # Everything is fine
      TRUE
    }
})

setMethod("initialize", "submodel",
  function(.Object,
           ...,
           formula = ~ 1,
           coefficients,
           vcov,
           centering = FLPar(centering = 0),
           distr = "norm",
           link = log,
           linkinv = exp
           ) {
      # initialize FLComp slots
      .Object <- callNextMethod(.Object, ...)
       # initialize remaining slots
      formula(.Object) <- formula
      if (!missing(coefficients)) {
        coef(.Object) <- coefficients
      }
      # need hard assignment first time round
      npar <- length(coef(.Object))
      parnames <- rownames(coef(.Object))
      .Object@vcov <- array(NA, dim = c(npar, npar, 1),
                            dimnames = list(parnames, parnames, 1))
      .Object@vcov[] <- diag(npar)

      # check dims in the following?
      if (!missing(vcov)) vcov(.Object) <- vcov
      .Object@centering <- centering
      .Object@distr <- distr
      .Object@link <- link
      .Object@linkinv <- linkinv
      .Object
})




#' @rdname submodel-class
#' @template Accessors
#' @template Constructors
#' @template bothargs
#' @aliases submodel submodel-methods
setGeneric("submodel", function(object, ...)
  standardGeneric("submodel"))
#' @rdname submodel-class
setMethod("submodel", signature(object = "missing"),
  function(...) {
    # empty
    if (missing(...)) {
      new("submodel")
    # or not
    } else {
      args <- list(...)
      args$Class <- "submodel"
      do.call("new", args)
    }
  }
)

#' @rdname submodel-class
setMethod("params", "submodel", function(object) object@coefficients)

#' @rdname submodel-class
#' @aliases sMod sMod-methods
setGeneric("sMod", function(object, ...) standardGeneric("sMod"))
#' @rdname submodel-class
setMethod("sMod", "submodel", function(object) object@formula)



#
#  show methods
#

setMethod("show", "submodel",
  function(object) {
    cat(paste0("\t", name(object),": ", sep = ""))
    print(formula(object), showEnv = FALSE)
  }
)


#
# Other methods
#


#' @rdname submodel-class
#' @param obj the object to be subset
#' @param it iteration to be extracted
setMethod("iter", "submodel",
  function(obj, it) {
    niters <- dim(vcov(obj))[3]
    # follow behaviour of iter see
    # showMethods(iter, classes = "FLPar", includeDefs = TRUE)
    if (niters == 1) {
      obj@vcov <- obj@vcov
    } else {
      obj@vcov <- obj@vcov[,, it, drop = FALSE]
    }
    obj@coefficients <- iter(obj@coefficients, it)
    obj@centering <- iter(obj@centering, it)
    obj
  }
)


#' @rdname submodel-class
#' @param iter the number of iterations to create
#' @param fill.iter should the new iterations be filled with values (TRUE) or NAs (FALSE)
setMethod("propagate", signature(object = "submodel"),
  function(object, iter, fill.iter = TRUE) {
    # propagate coefs and centering
    object@coefficients <-
      propagate(object@coefficients, iter, fill.iter = fill.iter)
    object@centering <-
      propagate(object@centering, iter, fill.iter = fill.iter)

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

#
#  duplicate accessors to help link with existing functions
#

#' @rdname submodel-class
#' @param x the submodel object that is to be modified
setMethod("formula", "submodel", function(x) x@formula)

#' @rdname coef-methods
#' @aliases coef<-,a4aFitSA-methods
setGeneric("formula<-", function(object, value) standardGeneric("formula<-"))

#' @rdname coef-methods
setMethod("formula<-", c("submodel", "formula"),
  function(object, value) {
    object@formula <- value

    # cannot deal with factors at the moment
    stopifnot(!formula_has_covariate_factors(value))

    # recalc coefficients
    df <- as.data.frame(object)

    vars <- all.vars(value)
    vars <- setdiff(vars, names(df))
    args <- sapply(vars, function(x) rnorm(nrow(df)), simplify = FALSE)
    if (length(vars)) {
      df <- cbind.data.frame(df, args)
    }
    Xmat <- getX(value, df, check = FALSE)

    coef(object) <- 
      FLPar(structure(rep(0, ncol(Xmat)),names = colnames(Xmat)))

    object
  }
)




setMethod("as.data.frame",
  signature(x = "submodel", row.names = "missing", optional = "missing"),
  function (x, ...) 
  {
    flq <- flq_from_range(x)
    df <- as.data.frame(flq)

    df
  }
)
