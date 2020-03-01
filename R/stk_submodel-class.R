#' @title Stock submodels class
#' @docType class
#' @name stk_submodel
#' @rdname stk_submodel-class
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
#' @aliases a4astkmodel-class

setClass(
  "stk_submodel",
  contains = c("FLComp", "submodels"),
  slots =
    c(
      "srmod" = "sr_submodel",
      "m" = "FLQuant",
      "mat" = "FLQuant",
      "stock.wt" = "FLQuant",
      "catch.wt" = "FLQuant",
      "harvest.spwn" = "FLQuant",
      "m.spwn" = "FLQuant"
    )
)

setMethod(
  "initialize",
  "stk_submodel",
  function(.Object, ..., fmod, n1mod, srmod) {
    # initialize FLComp slots
    .Object <- callNextMethod(.Object, ...)

    # fill in missing FLComp slots
    if (length(name(.Object)) == 0 || name(.Object) == "") {
      name(.Object) <- "stk_submodel"
    }

    # there must be at least 2 years, to allow for an SR model
    rng <- range(.Object)
    if (rng["maxyear"] == rng["minyear"]) {
      warning("increasing maxyear by 1")
      .Object@range["maxyear"] <- .Object@range["maxyear"] + 1
    }

    # fill in submodels
    if (missing(fmod)) {
      fmod <-
        submodel(
          range = range(.Object)
        )
    }
    if (missing(n1mod)) {
      n1range <- range(.Object)[c("min", "max", "minyear", "maxyear")]
      n1range["maxyear"] <- n1range["minyear"]
      n1mod <-
        submodel(
          range = n1range
        )
    }
    if (missing(srmod)) {
      srrange <- range(.Object)[c("min", "max", "minyear", "maxyear")]
      srrange["max"] <- srrange["min"]
      srrange["minyear"] <- srrange["minyear"] + 1
      srmod <-
        sr_submodel(
          range = srrange
        )
    }

    # convert srmod to submodels, and add into main submodels
    # note rmod is a placeholder for the realised recruitments
    .Object@.Data <-
      list(
        fmod = submodel(fmod, name = "fmod"),
        n1mod = submodel(n1mod, name = "n1mod"),
        rmod =
          submodel(
            srmod$a, name = "rmod", formula = ~ factor(year) - 1
          ),
        sramod = submodel(srmod$a, name = "sramod"),
        srbmod = submodel(srmod$b, name = "srbmod")
      )

    tmp <- new("submodels", .Object@.Data)
    .Object@corBlocks <- tmp@corBlocks

    .Object@srmod <- srmod

    # fill in missing m, stock.wt or mat, etc.
    dots <- list(...)
    flqnames <-
      c(
        "m", "mat", "stock.wt", "catch.wt", "harvest.spwn", "m.spwn"
      )
    if (!all(flqnames %in% names(dots))) {
      range <- range(.Object)
      df <-
        expand.grid(
          age = range["min"]:range["max"],
          year = range["minyear"]:range["maxyear"],
          data = NA,
          stringsAsFactors = FALSE
        )
      emptyflq <- as.FLQuant(df)
    }
    for (flqname in flqnames) {
      if (!flqname %in% names(dots)) {
        slot(.Object, flqname) <- emptyflq
      }
    }

    .Object
  }
)

#' @rdname stk_submodel-class
#' @aliases stk_submodel stk_submodel-methods
#' @template Accessors
#' @template Constructors
setGeneric("stk_submodel", function(object, ...) {
  standardGeneric("stk_submodel")
})

#' @rdname submodel-class
setMethod(
  "stk_submodel", signature(object = "missing"),
  function(...) {
    # empty
    if (missing(...)) {
      new("stk_submodel")
      # or not
    } else {
      args <- list(...)
      args$Class <- "stk_submodel"
      do.call("new", args)
    }
  }
)


# setValidity(
#  "stk_submodels",
#  function(object) {
#    # the following is not throwing errors somehow...
#    #if (length(object) == 0) {
#    #  "There must be 3 submodels named f, n1, and sr"
#    #} else
#    # check dimensions of m and wt are the same as f submodel
#    if (length(.Object) < 3) {
#      stop("There must be 3 submodels named f, n1, and sr")
#    }
#    if (
#      !identical(
#        unlist(dimnames(object@m)[2:5]),
#        unlist(dimnames(object@stock.wt)[2:5])
#        )) {
#      "m and wt elements must share dimensions 2 to 5"
#    } else {
#      TRUE
#    }
#  }
# )

#
# coerce methods
#

# coerce from FLStock to stk_submodel

#
#  accessor methods
#

#' @rdname stk_submodels-class
setMethod(
  "m",
  "stk_submodel",
  function(object) object@m
)

#' @rdname stk_submodels-class
setMethod(
  "mat",
  "stk_submodel",
  function(object) object@mat
)

#' @rdname stk_submodels-class
setMethod(
  "stock.wt",
  "stk_submodel",
  function(object) object@stock.wt
)

#' @rdname stk_submodels-class
setMethod(
  "catch.wt",
  "stk_submodel",
  function(object) object@catch.wt
)

#' @rdname stk_submodels-class
setMethod(
  "harvest.spwn",
  "stk_submodel",
  function(object) object@harvest.spwn
)

#' @rdname stk_submodels-class
setMethod(
  "m.spwn",
  "stk_submodel",
  function(object) object@m.spwn
)


#
#  show methods
#

#' @rdname stk_submodels-class
setMethod(
  "show",
  "stk_submodel",
  function(object) {
    cat("a4a stock model for:", object@name, "\n")
    show(as(object, "submodels"))
  }
)


#
# Coercion
#

setMethod(
  "as.data.frame",
  signature(
    x = "stk_submodel", row.names = "missing", optional = "missing"
  ),
  function(x, drop = FALSE, fill = FALSE, centering = FALSE, ...) {
    as.data.frame(as(x, "submodels"))
  }
)


#
# Other methods
#

#' @rdname stk_submodels-class
#' @param obj the object to be subset
#' @param it iteration to be extracted
setMethod(
  "iter",
  "stk_submodel",
  function(obj, it) {
    message("todo - update ")
    obj@vcov <- obj@vcov[, , it, drop = FALSE]
    obj@m <- iter(obj@m, it)
    obj@wt <- iter(obj@wt, it)
    obj
  }
)
