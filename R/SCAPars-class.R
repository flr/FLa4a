
#' a4aFitSA extends \code{"a4aFit"} class.
#'
#' Some details about this class and my plans for it in the body.
#'
#' \describe{
#'    \item{myslot1}{A logical keeping track of something.}
#'
#'    \item{myslot2}{An integer specifying something else.}
#' 
#'    \item{myslot3}{A data.frame holding some data.}
#'  }
#' @name SCAPars-class
#' @rdname SCAPars-class
#' @exportClass SCAPars
setClass("SCAPars",
        representation(
               stkmodel    = "a4aStkParams",
                 qmodel    = "submodels",
                 vmodel    = "submodels"),
        prototype = prototype(
               stkmodel    = new("a4aStkParams"),
                 qmodel    = new("submodels"),
                 vmodel    = new("submodels"))
)


# constructor
#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("SCAPars", function(object, ...) standardGeneric("SCAPars"))

#' Title 
#' @name SCAPars
#' @docType methods
#' @rdname SCAPars-methods
#' @aliases SCAPars,SCAPars-method
setMethod("SCAPars", signature(object="missing"),
  function(...) {
    # empty
  	if(missing(...)){
	  	new("SCAPars")
    # or not
  	} else {
      args <- list(...)
	  args$Class <- 'SCAPars'
      do.call("new", args)
	  }
  }
)

# accessors


#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("stkmodel", function(object, ...) standardGeneric("stkmodel"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("stkmodel", "SCAPars", function(object) object @ stkmodel)


#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("n1model", function(object, ...) standardGeneric("n1model"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("n1model", "SCAPars", function(object) object @ stkmodel @ n1Mod)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("srmodel", function(object, ...) standardGeneric("srmodel"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("srmodel", "SCAPars", function(object) object @ stkmodel @ srMod)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("fmodel", function(object, ...) standardGeneric("fmodel"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("fmodel", "SCAPars", function(object) object @ stkmodel @ fMod)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("qmodel", function(object, ...) standardGeneric("qmodel"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("qmodel", "SCAPars", function(object) object@qmodel)


#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("vmodel", function(object, ...) standardGeneric("vmodel"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("vmodel", "SCAPars", function(object) object@vmodel)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("srPars", function(object, ...) standardGeneric("srPars"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("srPars", "SCAPars", function(object) object@srmodel@pars)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("srCovar", function(object, ...) standardGeneric("srCovar"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("srCovar", "SCAPars", function(object) object@srmodel@covar)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("srFrml", function(object, ...) standardGeneric("srFrml"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("srFrml", "SCAPars", function(object) object@srmodel@model)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("fPars", function(object, ...) standardGeneric("fPars"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("fPars", "SCAPars", function(object) object@fmodel@pars)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("fCovar", function(object, ...) standardGeneric("fCovar"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("fCovar", "SCAPars", function(object) object@fmodel@covar)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("fFrml", function(object, ...) standardGeneric("fFrml"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("fFrml", "SCAPars", function(object) object@fmodel@model)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("qPars", function(object, ...) standardGeneric("qPars"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("qPars", "SCAPars", function(object) object@qmodel@pars)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("qCovar", function(object, ...) standardGeneric("qCovar"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("qCovar", "SCAPars", function(object) object@qmodel@covar)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("qFrml", function(object, ...) standardGeneric("qFrml"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("qFrml", "SCAPars", function(object) object@qmodel@model)



#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("vPars", function(object, ...) standardGeneric("vPars"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("vPars", "SCAPars", function(object) object@vmodel@pars)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("vCovar", function(object, ...) standardGeneric("vCovar"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("vCovar", "SCAPars", function(object) object@vmodel@covar)

#' Calculate the median accross iterations
#'
#' @param object an FLQuant with iters
#'
#' @param ... Additional argument list that might not ever
#'  be used.
#'
#' @return an FLQuant
#' 
#' @seealso \code{\link{print}} and \code{\link{cat}}
#' 
#' @export
#' @docType methods
#' @rdname a4aFitSA-methods
#'
#' @examples
#' data(ple4)
setGeneric("vFrml", function(object, ...) standardGeneric("vFrml"))

#' Title 
#' @name a4aFitSA
#' @docType methods
#' @rdname a4aFitSA-methods
#' @aliases a4aFitSA,a4aFitSA-method
setMethod("vFrml", "SCAPars", function(object) object@vmodel@model)






