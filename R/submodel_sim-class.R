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
setClass(
  "submodel_sim",
  contains = "submodel",
  slots = c(covariates = "data.frame")
)

setMethod(
  "initialize", "submodel_sim",
  function(.Object, ..., covariates) {
    .Object <- callNextMethod(.Object, ...)

    # do assignments
    if (!missing(covariates)) {
      .Object@covariates <- covariates
    } else if (identical(.Object@covariates, new("data.frame"))) {
      .Object@covariates <- as.data.frame(.Object, drop = TRUE)
    }

    .Object
  }
)


#
#  Convert to data.frame
#

setMethod("as.data.frame",
  signature(x = "submodel_sim", row.names = "missing", optional = "missing"),
  function (x, drop = FALSE, fill = FALSE, ...)
  {
    df <-
      callNextMethod(
        x, drop = drop, fill = nrow(x@covariates) == 0, ...
      )

    # add in covariates
    if (nrow(x@covariates) > 0) {
      df$`_id` <- 1:nrow(df)
      df <-
        merge(
          df,
          x@covariates,
          all.x = TRUE, all.y = FALSE,
          sort = FALSE
        )
      df <- df[order(df$`_id`), setdiff(names(df), "_id")]
    }

    df
  }
)
