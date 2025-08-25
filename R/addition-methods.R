#==================================================================== 
# "+" methods
#==================================================================== 

#' + methods
#' @name addition
#' @description Update \code{FLStock} and \code{FLIndex} objects with stock assessment results.
#' @param e1 the original \code{FLStock} or \code{FLIndex} object
#' @param e2 a \code{a4aFit} object from where the new \code{FLStock} or \code{FLIndex} slots will be extracted.
#' @details If both objects have the same number of iterations, the \code{FLStock} slots will be replaced by the \code{a4aFit} slots, in the case of 1 iter, or \code{a4aFitSA} slots, in the case of n iters. If one of the objects has 1 iter and the other n, the method will simulate using the fit results from the \code{a4aFitSA} object to update the slots of the \code{FLStock} object.
#' @rdname addition-methods
#' @aliases +,FLStock,a4aFit-method +,FLIndices,a4aFit-method
setMethod("+", c("FLStock", "a4aFit"), function(e1, e2)
{

  nit1 <- dims(e1)$iter
  nit2 <- dims(e2)$iter
  if(nit1>nit2) {
	e2 <- propagate(e2, nit1)
  } else if(nit1<nit2){
  	e1 <- propagate(e1, nit2)
  }

  stock.n(e1) <- stock.n(e2)
  landings.n(e1) <- landings.n(e1) * (catch.n(e2)/catch.n(e1))
  discards.n(e1) <- discards.n(e1) * (catch.n(e2)/catch.n(e1))
  catch.n(e1) <- catch.n(e2)
  harvest(e1) <- harvest(e2)
  
  catch(e1) <- computeCatch(e1, na.rm=FALSE)
  stock(e1) <- computeStock(e1, na.rm=FALSE)
  landings(e1) <- computeLandings(e1, na.rm=FALSE)
  discards(e1) <- computeDiscards(e1, na.rm=FALSE)
  
  e1
})

setMethod("+", c("FLIndices", "a4aFit"), function(e1, e2) 
{

  #niters <- dims(e1) $ iter
  #if (niters > 1) stop("adding a basic a4aFit object only makes sence with 1 iteration")

  for (i in seq(FLIndices)) {
    index(e1[[i]]) <- index(e2)[[i]]
    #catch.n(e1[[i]]) <- index(e1[[i]]) * effort(e1[[1]])
    #index.q(e1[[1]])
    #sel.pattern(e1[[1]]) 
    #??index.var(e1[[1]])
  }
    
  e1
})


#' + methods
#' @name addition
#' @description Update \code{FLStocks} objects with multiple stock assessment results in a \code{a4aFits}.
#' @param e1 the original \code{FLStocks} object
#' @param e2 a \code{a4aFits} object from where the new \code{FLStock} slots will be extracted.
#' @details If both objects have the same number of iterations, the \code{FLStocks} slots will be replaced by the \code{a4aFits} slots, in the case of 1 iter, or \code{a4aFitSA} slots, in the case of n iters. If one of the objects has 1 iter and the other n, the method will simulate using the fit results from the \code{a4aFitSA} object to update the slots of the \code{FLStock} object.
#' @rdname addition-methods
#' @aliases +,FLStocks,a4aFits-method
setMethod("+", c("FLStocks", "a4aFits"), function(e1, e2)
{

	# checks 1 or n
	ns <- length(e1)
	nf <- length(e2)
	if(ns!=1 & nf!= 1 & ns!=nf) stop("objects must be of equal size or of size 1")

	# set same sizes
	n <- max(ns, nf)
	if(n>1 & ns==1){ 
		e1[1:n] <- e1[1]
		names(e1) <- rep(names(e1[1]), n)
		} else if(n>1 & nf==1){
		e2[1:n] <- e2[1]
		names(e2) <- rep(names(e2[1]), n)
		}
	
	# call +
	for(i in 1:n) e1[[i]] <- e1[[i]] + e2[[i]]

	# out	
	e1
})

#' + methods
#' @name addition
#' @description Add residual uncertainty to \code{FLStock} objects.
#' @param e1 the original \code{FLStocks} object
#' @param e2 a \code{a4aFitResiduals} object with residuals (suggest to use "deviances").
#' @details Random normal draws will be added to the log transformed stock.n and catch.n slots using the standard deviation computed from the residuals catch.n and index, for catch.n and stock.n respectively. Returns a \code{FLStock} object.
#' @rdname addition-methods
#' @aliases +,FLStocks,a4aFitResiduals-method
setMethod("+", c("FLStock", "a4aFitResiduals"), function(e1, e2)
{
  # if not deviances sd is 1
  if(e2@desc!="deviances") warning("Residuals are not deviances, in which case sd=1.")

  # number of iterations will come from stock object
  nit1 <- dims(e1)$iter
  # catch.n standard deviation
  flqsdc <- e1@catch.n
  flqsdc[] <- sqrt(yearVars(e2$catch.n, na.rm=TRUE))
  e1@catch.n <- rlnorm(nsim, log(e1@catch.n), flqsdc)
  # return
  e1
})


#==================================================================== 
# "*" methods
#==================================================================== 

#' * methods
#' @name *
#' @description Update \code{FLStock} and \code{FLIndex} objects with simulations from stock assessment fits.
#' @param e1 the original \code{FLStock} or \code{FLIndex} object
#' @param e2 a \code{a4aFit} object from where the new \code{FLStock} or \code{FLIndex} slots will be extracted.
#' @rdname multiplication-methods
#' @aliases *,FLStock,a4aFitSA-method *,FLIndices,a4aFitSA-method
setMethod("*", c("FLStock", "a4aFitSA"), function(e1, e2) 
{
  niters <- dims(e1)$iter
  niters2 <- dims(e2)$iter
  if (niters==niters2) {
    nsim = niters
  } else {
    nsim = max(c(niters, niters2))
  }
  e2 <- simulate(e2, nsim = nsim)  
  e1 + e2	
})

#' @rdname multiplication-methods
setMethod("*", c("FLStock", "SCAPars"), function(e1, e2) 
{
stop("method * FLStock, SCAPars deprecated !\n")
})

#' @rdname multiplication-methods
setMethod("*", c("FLIndices", "a4aFitSA"), function(e1, e2) 
{
  e1 * pars(e2)
})

#' @rdname multiplication-methods
setMethod("*", c("FLIndices", "SCAPars"), function(e1, e2) 
{

#  niters <- dims(e1) $ iter
#  niters2 <- dim(e2 @ stkmodel @ params)[2]
#  if (niters > 1 & niters2 == 1) {
#    nsim = niters
#  } else {
#    nsim = 1
#    if (niters > niters2) stop("oh oh")
#    if (niters == 1 & niters2 > 0) {
#      niters <- niters2
#      e1 <- propagate(e1, niters)
#    }
#  }

  stop("not implemented yet")

  #FLIndices(out)
})

