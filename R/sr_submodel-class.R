#' @title Stock Recuirtment submodel class
#' @docType class
#' @name sr_submodel
#' @rdname sr_submodel-class
#' @template ClassDescription
#' @section Slot:
#' \describe{
#'  \item{\code{fmodel}}{F submodel \code{formula}}
#'  \item{\code{n1model}}{first year N \code{formula}}
#'  \item{\code{srmodel}}{stock-recruitment submodel \code{formula}}
#'  \item{\code{m}}{natural mortality \code{FLQuant}}
#'  \item{\code{wt}}{stock weights \code{FLQuant}}
#'  \item{\code{units}}{data units \code{character}}
#' }

setClass(
  "sr_submodel",
  contains = c("FLComp", "submodels"),
  slots =
    c(
      "formula" = "formula",
      "srr" = "character",
      "ID" = "integer",
      "CV" = "numeric",
      "SPR0" = "numeric"
    )
)

setMethod(
  "initialize",
  "sr_submodel",
  function(.Object, ..., formula) {
    # initialize FLComp slots
    .Object <- callNextMethod(.Object, ...)

    # initialise object
    if (!missing(formula)) {
      formula(.Object) <- formula
    } else if (identical(.Object@formula, new("formula"))) {
      formula(.Object) <- ~ none()
    } else {
      # seems odd, but resets any wonky settings
      # if a user has messed with the internals
      formula(.Object) <- formula(.Object)
    }

    # fill in missing FLComp slots
    if (length(name(.Object)) == 0 || name(.Object) == "") {
      name(.Object) <- "sr_submodel"
    }

    .Object
  }
)


#' @rdname stk_submodel-class
#' @aliases stk_submodel stk_submodel-methods
#' @template Accessors
#' @template Constructors
setGeneric("sr_submodel", function(object, ...) {
  standardGeneric("sr_submodel")
})

#' @rdname submodel-class
setMethod(
  "sr_submodel", signature(object = "missing"),
  function(...) {
    # empty
    if (missing(...)) {
      new("sr_submodel")
      # or not
    } else {
      args <- list(...)
      args$Class <- "sr_submodel"
      do.call("new", args)
    }
  }
)


#' @rdname formula-methods
setMethod(
  "formula<-", c("sr_submodel", "formula"),
  function(object, value) {

    object@formula <- value

    # process special sr model formula
    a4as <- isPresenta4aSRmodel(value)
    if (sum(a4as) > 1) {
      stop(
        "you can only specify one type of stock recruit relationship."
      )
    }
    if (sum(a4as) > 0 && any(!a4as)) {
      # ignore trailing elements
      message("Note: Trailing formula elements in the sr_submodel formula have been removed")
    }
    # what if no sr model specified?
    if (sum(a4as) == 0) {
      # put formula into a parameter of none()
    }
    srrmod <- geta4aSRmodel(value)
    srr_info <- eval(parse(text = srrmod))

    # fill in submodels for a and b parameters
    amod <-
      submodel(
        name = "a",
        range = range(object),
        formula = srr_info$a
      )
    bmod <-
      submodel(
        name = "b",
        range = range(object),
        formula = srr_info$b
      )

    object@.Data <- list(a = amod, b = bmod)

    # fill in corBlocks
    tmp <- new("submodels", object@.Data)
    object@corBlocks <- tmp@corBlocks

    # fill in other slots
    object@srr <- srr_info$srr
    object@SPR0 <- srr_info$SPR0
    object@CV <- srr_info$srrCV
    object@ID <- as.integer(srr_info$ID)

    object
  }
)


#setValidity


# coerce methods



#  accessor methods


#
#  show methods
#

#' @rdname stk_submodels-class
setMethod(
  "show",
  "sr_submodel",
  function(object)
  {
    cat("a4a stock recruitment model for:", object@name, "\n")
    cat("srr: ", object@srr)
    show(as(object, "submodels"))
  }
)


#
# Coercion
#

setMethod(
  "as.data.frame",
  signature(
    x = "sr_submodel", row.names = "missing", optional = "missing"
  ),
  function(x, drop = FALSE, fill = FALSE, centering = FALSE, ...) {
    as.data.frame(as(x, "submodels"))
  }
)


#
# Other methods
#

# propagate