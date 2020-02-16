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
#'  \item{\code{vcov}}{\code{array} with variance covariance paramaters related to the variance model}
#'
#'  \item{\code{centering}}{\code{numeric} value used for centering the data}
#'
#'  \item{\code{distr}}{a character with the parameters' statistical distribution; it must match a known distribution for R (\emph{e.g.} "norm" for gaussian) so that \code{rnorm} can be called}
#' }
#' @aliases submodel-class
setClass(
  "submodel",
  contains = "FLComp",
  slots =
    c(
      formula = "formula",
      coefficients = "FLPar",
      vcov = "array",
      centering = "FLPar",
      distr = "character",
      link = "function",
      linkinv = "function"
    )
)

setValidity(
  "submodel",
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
  }
)

# show
setMethod(
  "show", "submodel",
  function(object) {
    tab <- "  "
    cat(name(object), ":\n", sep = "")

    cat(tab, "formula: ", sep = "")
    print(formula(object), showEnv = FALSE)

    rng <- range(object)
    cat(
      tab,
      "range: (", rng["min"], " - ", rng["max"], "), ",
      "years: (", rng["minyear"], " - ", rng["maxyear"], ").\n",
      sep = ""
    )

    coefs <- coef(object)
    dimnames <- dimnames(coefs)[1]
    names(dimnames) <- paste0(tab, "coef, iters = ", dim(coefs)[2], ":")
    if (dim(coefs)[2] > 1) {
      v1 <- apply(coefs@.Data, 1, median, na.rm = TRUE)
      v2 <- apply(coefs@.Data, 1, mad, na.rm = TRUE)
      v3 <-
        paste0(
          format(v1, digits = 5),
          "(", format(v2, digits = 3), ")"
        )
    } else {
      v3 <- format(coefs@.Data, digits = 3)
    }
    print(
      array(
        v3,
        dim = dim(coefs)[1],
        dimnames = dimnames
      ),
      quote = FALSE
    )
  }
)

# initialize
setMethod(
  "initialize",
  "submodel",
  function(.Object, ..., formula, coefficients, vcov, centering) {
    # initialize FLComp slots
    .Object <- callNextMethod(.Object, ...)

    # fill in missing FLComp slots
    if (length(name(.Object)) == 0 || name(.Object) == "") {
      name(.Object) <- "submodel"
    }

    # initialise submodel
    if (!missing(formula)) {
      formula(.Object) <- formula
    } else {
      formula(.Object) <- ~1
    }

    # set up coefficients and vcov structure
    if (!missing(coefficients)) {
      # warn if coefficients are wrong length
      tryCatch(
        .Object@coefficients[] <- coefficients,
        error = function(e) {
          warning(
            "Error assigning coefficients, check length of vector.",
            call. = FALSE
          )
        }
      )
    }
    if (!missing(vcov)) {
      tryCatch(
        .Object@vcov[] <- vcov,
        error = function(e) {
          warning(
            "Error assigning vcov, check dimensions.",
            call. = FALSE
          )
        }
      )
    }
    .Object@centering <-
      FLPar(structure(0, names = "centering"))
    if (!missing(centering)) {
      tryCatch(
        .Object@centering[] <- centering,
        error = function(e) {
          warning(
            "Error assigning centering.",
            call. = FALSE
          )
        }
      )
    }

    # process left over dots
    dots <- list(...)
    if.null <- function(x, value) if (is.null(x)) value else x
    .Object@distr <- if.null(dots$distr, "norm")
    .Object@link <- if.null(dots$link, log)
    .Object@linkinv <- if.null(dots$linkinv, exp)

    .Object
  }
)

#' @rdname submodel-class
#' @template Accessors
#' @template Constructors
#' @template bothargs
#' @aliases submodel submodel-methods
setGeneric("submodel", function(object, ...) {
  standardGeneric("submodel")
})
#' @rdname submodel-class
setMethod(
  "submodel", signature(object = "missing"),
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
#' @param obj the object to be subset
#' @param it iteration to be extracted
setMethod("iter", "submodel", function(obj, it) {
  niters <- dim(vcov(obj))[3]
  # follow behaviour of iter see
  # showMethods(iter, classes = "FLPar", includeDefs = TRUE)
  if (niters == 1) {
    obj@vcov <- obj@vcov
  } else {
    obj@vcov <- obj@vcov[, , it, drop = FALSE]
  }
  obj@coefficients <- iter(obj@coefficients, it)
  obj@centering <- iter(obj@centering, it)
  obj
})

#' @rdname submodel-class
#' @param iter the number of iterations to create
#' @param fill.iter should the new iterations be filled with values (TRUE) or NAs (FALSE)
setMethod(
  "propagate", signature(object = "submodel"),
  function(object, iter, fill.iter = TRUE) {

    # propagate coefs and centering
    object@coefficients <- propagate(object@coefficients, iter, fill.iter = fill.iter)
    object@centering <- propagate(object@centering, iter, fill.iter = fill.iter)

    # now propagate vcov
    vcov.iter <- vcov(object)
    dob <- dim(vcov.iter)

    if (iter != dob[3]) {
      # CHECK no iters in object
      if (dob[3] > 1) stop("propagate can only extend objects with no iters")

      object@vcov <- array(NA, dim = c(dob[1:2], iter), dimnames = c(dimnames(vcov.iter)[1:2], list(1:iter)))
      if (fill.iter) {
        object@vcov[] <- as.vector(vcov.iter)
      } else {
        object@vcov[, , 1] <- as.vector(vcov.iter)
      }
    }

    object
  }
)


#' @rdname submodel-class
#' @param x the submodel object that is to be modified
setMethod("formula", "submodel", function(x) x@formula)

#' @rdname formula-methods
#' @aliases formula<-,a4aFitSA-methods
setGeneric(
  "formula<-",
  function(object, value) standardGeneric("formula<-")
)

#' @rdname formula-methods
setMethod(
  "formula<-", c("submodel", "formula"),
  function(object, value) {
    object@formula <- value

    # cannot deal with factors at the moment
    stopifnot(!formula_has_covariate_factors(object@formula))

    # get data.frame with filled out columns for covars
    df <- as.data.frame(object, fill = TRUE)

    # recalc coefficients
    Xmat <- getX(value, df, check = FALSE)

    object@coefficients <-
      FLPar(structure(rep(0, ncol(Xmat)), names = colnames(Xmat)))

    object@vcov <-
      array(
        NA,
        dim = c(ncol(Xmat), ncol(Xmat), 1),
        dimnames = list(colnames(Xmat), colnames(Xmat), 1)
      )
    # set as diagonal to begin with
    object@vcov[,,1] <- diag(ncol(Xmat))

    object
  }
)


# as data.frame
# fill = TRUE means extra columns added to cover covars in formula
setMethod(
  "as.data.frame",
  signature(
    x = "submodel", row.names = "missing", optional = "missing"
  ),
  function(x, drop = FALSE, fill = FALSE, ...) {
    flq <- flq_from_range(x)
    df <- as.data.frame(flq, drop = drop)

    if (fill) {
      # check formula for extra covars to include
      vars <- all.vars(formula(x))
      vars <- setdiff(vars, names(df))
      args <-
        sapply(
          vars,
          function(x) rnorm(nrow(df)),
          simplify = FALSE)
      if (length(vars)) {
        df <- cbind.data.frame(df, args)
      }
    }

    # add centering if present
    cdf <- as.data.frame(x@centering, drop = FALSE)
    iter_idx <- as.numeric(df$iter)
    if (is.null(iter_idx)) iter_idx <- rep(1, nrow(df))


    cbind.data.frame(
      df[names(df) != "data"],
      centering = cdf$data[iter_idx]
    )
  }
)


# coef methods

setMethod(
  "coef", "submodel",
  function(object, ...) object@coefficients
)

setGeneric("coef<-", function(object, value) standardGeneric("coef<-"))

setMethod(
  "coef<-", c("submodel", "FLPar"),
  function(object, value) {
    # retain structure
    object@coefficients[] <- value
    object
  }
)

setMethod(
  "coef<-", c("submodel", "numeric"),
  function(object, value) {
    # retain structure
    object@coefficients[] <- value
    object
  }
)


# range methods

setMethod(
  "range<-",
  c("submodel", "ANY", "numeric"),
  function(x, i, value) {
    warning(
      "changing the range may change the fitted values if smoothers",
      "are in the formula"
    )

    x@range[i] <- value
    x
  }
)


setMethod(
  "range<-",
  c("submodel", "missing", "numeric"),
  function(x, i, value) {
    warning(
      "changing the range may change the fitted values if smoothers",
      "are in the formula"
    )
    x@range[names(value)] <- value
    x
  }
)




# methods for back compatability

#' @rdname submodel-class
setMethod("params", "submodel", function(object) object@coefficients)

#' @rdname submodel-class
#' @aliases sMod sMod-methods
setGeneric("sMod", function(object, ...) standardGeneric("sMod"))
#' @rdname submodel-class
setMethod("sMod", "submodel", function(object) object@formula)
