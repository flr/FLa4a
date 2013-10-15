
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
#' @rdname logq-methods
#'
#' @examples
#' data(ple4)
#' genFLQuant(harvest(ple4), method = "ac")
setGeneric("FLblockboot", function(object, ...) standardGeneric("FLblockboot"))

#' @rdname logq-methods
#' @aliases logq,FLa4aFit-method
setMethod("FLblockboot", signature(object = "FLQuant"),
  function(object, R = 999) {
    	
    dms <- dim(object)
    arr <- aperm(object @ .Data, c(2,1,3:6))
	
    make.ends <- function (a, n) 
    {
      mod <- function(i, n) 1 + (i - 1)%%n
      if (a[2L] == 0) {
        numeric()
      } else {
        mod(seq.int(a[1L], a[1L] + a[2L] - 1, length.out = a[2L]), n)
      }
    }

    # sample blocsk of years 
    n <- dms[2]
  
    # find block length:
    dim(arr) <- c(dms[2], prod(dms[-2]))
    out <- b.star(arr)[,1]
    dim(out) <- dms[-2]
    l <- c(apply(out, 3:5, max))

    # reshape for boot strapping
    dim(arr) <- c(dms[2], prod(dms[c(1,3,4,5)]), dms[6])

    # loop over iters

    doone <-
      function(arri, li) {

        len.tot <- rep(0, R)
        lens <- NULL; cont <- TRUE
        while (cont) {
          temp <- 1 + rgeom(R, 1/li)
          temp <- pmin(temp, n - len.tot)
          lens <- cbind(lens, temp)
          len.tot <- len.tot + temp
          cont <- any(len.tot < n)
        }
        dimnames(lens) <- NULL 
  
        i.a <- list(starts  = matrix(sample.int(n, ncol(lens) * R, replace = TRUE), R), 
                    lengths = lens)
       
        inds <- 
          sapply(seq_len(R),
               function(r) {
                 ends <- cbind(i.a$starts[r, ], i.a$lengths[r, ])
                 inds <- apply(ends, 1L, make.ends, n)
                 inds <- if (is.list(inds)) matrix(unlist(inds)[1L:n], n, 1L) else matrix(inds, n, 1L)
              })
    
        res <- propagate(object, R)
        # pre-whiten  and post-blacken ..
        arri <- matrix(arri, nrow = dms[2])
        barri <- apply(arri, 2, function(x) fitted(loess(x ~ y, data = list(y = 1:length(x)))))
        res[] <- as.vector(t(arri - barri)[,c(inds)] + c(t(barri)))
        res
      }
      

  
  
    if (dms[6] > 1) {
      out <- lapply(1:dms[6], function(i) doone(arr[, , i], l[i]) )
      names(out) <- paste0("iter", 1:dms[6])
      FLQuants(out)
    } else {
      doone(arr[, , 1], l[1])
    }
  }
)




#' Calculate the optimal block length for stationary block bootstrap
#'
#' @param object an FLQuant
#'
#' @return an FLQuant containing the optimal block lengths for each year
#' 
#' @export
#' @docType methods
#' @rdname FLBstar-methods
#'
#' @examples
#' data(ple4)
#' FLBstar(catch.n(ple4))
setGeneric("FLBstar", function(object, ...) standardGeneric("FLBstar"))

#' @rdname FLBstar-methods
#' @aliases FLBstar,FLQuant-method
setMethod("FLBstar", signature(object = "FLQuant"),
  function (object) {

    dms <- dim(object)
    arr <- aperm(object @ .Data, c(2,1,3:6))
    dim(arr) <- c(dms[2], prod(dms[-2]))

    out <- b.star(arr)[,1]
  
    dms[2] <- 1
    dmns <- dimnames(object)
    dmns[[2]] <- 1
    dim(out) <- dms
    dimnames(out) <- dmns
  
    FLQuant(out)
  })

