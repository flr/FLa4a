#' @title Submodel_sim class
#' @docType class
#' @name submodel_sim
#' @rdname submodel_sim-class
#' @template ClassDescription
#' @section Slot:
#' \describe{
#'
#'  \item{\code{formula}}{\code{formula} describing the model}
#'
#'  \item{\code{coefficients}}{\code{FLPar} with model parameters}
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
#' 
#'  \item{\code{link}}{\code{function} value used for centering the data}
#' 
#'  \item{\code{linkinv}}{\code{function} value used for centering the data}
#' 
#'  \item{\code{covariates}}{\code{FLQuants} list of FLQuant objects holding}
#' 
#' }
#' @aliases submodel_sim-class
setClass("submodel_sim",
  contains = c("FLComp", "submodel"),
  slots = c(covariates   = "FLQuants")
)

setValidity("submodel_sim",
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

setMethod("initialize", "submodel_sim",
    function(.Object,
             ...,
            covariates = FLQuants()) {
      .Object <- callNextMethod(.Object, ...)
      .Object@covariates <- covariates
      .Object
})




#' @rdname submodel-class
#' @template Accessors
#' @template Constructors
#' @template bothargs
#' @aliases submodel submodel-methods
setGeneric("submodel_sim", function(object, covariates, ...)
  standardGeneric("submodel_sim"))
#' @rdname submodel-class
setMethod("submodel_sim", signature(object = "missing", covariates = "missing"),
  function(...) {
    # empty
    if (missing(...)) {
      new("submodel_sim")
    # or not
    } else {
      args <- list(...)
      args$Class <- "submodel_sim"
      do.call("new", args)
    }
  }
)




#
#  Convert to data.frame
#

setMethod("as.data.frame",
  signature(x = "submodel_sim", row.names = "missing", optional = "missing"),
  function (x, ...) 
  {
    flq <- flq_from_range(x)
    df <- as.data.frame(flq)

    # add in covariates
    if (length(x@covariates) > 0) {
      covar <- x@covariates
      # add in covariates to data.frame
      tmp <- lapply(seq_along(covar), function(i) {
        x <- as.data.frame(covar[[i]])[c("age", "year", "data")]
        if (length(unique(x$age)) == 1) x <- x[names(x) != "age"]
        if (length(unique(x$year)) == 1) x <- x[names(x) != "year"]
        names(x) <- gsub("data", names(covar)[i], names(x))
        x
      })
      covar.df <- tmp[[1]]
      for (i in seq_along(tmp[-1]))
        covar.df <- merge(covar.df, i, all = TRUE, sort = FALSE)

      df <- merge(df, covar.df, all.x = TRUE, all.y = FALSE)
    }

    df
  }
)
