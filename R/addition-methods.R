#==================================================================== 
# "+" methods
#==================================================================== 

#' + methods
#' @name +
#' @description Update \code{FLStock} and \code{FLIndex} objects with stock assessment results
#' @details If both objects have the same number of iterations, the \code{FLStock} slots will be replaced by the \code{a4aFit} slots, in case of 1 iter, or \code{a4aFitSA} slots, in case of n iters. If one of the objects has 1 iter and the other n, the method will simulate using the fit results from the \code{a4aFitSA} object, to update the slots of the \code{FLStock} object.
#' @rdname addition-methods
#' @aliases +,FLStock,a4aFit-method
setMethod("+", c("FLStock", "a4aFit"), function(e1, e2)
{

  niters <- dims(e1) $ iter
  if (niters > 1) stop("adding a basic a4aFit object only makes sence with 1 iteration")

  stock.n(e1) <- stock.n(e2)
  catch.n(e1) <- catch.n(e2)
  harvest(e1) <- harvest(e2)
  
  catch(e1) <- computeCatch(e1)
  stock(e1) <- computeStock(e1)
  
  e1
})

#' @rdname addition-methods
#' @aliases +,FLStock,a4aFit-method
setMethod("+", c("FLStock", "a4aFitSA"), function(e1, e2)
{
  nit1 <- dims(e1) $ iter
  nit2 <- dim(qmodel(pars(e2))[[1]]@params)["iter"]
  v <- c(nit1, nit2)
  if(min(v)==max(v)) e1 <- e1 + as(e2, "a4aFit") else e1 <- e1 * e2
  e1
})

#' @rdname addition-methods
#' @aliases +,FLIndices,a4aFit-method
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

#==================================================================== 
# "*" methods
#==================================================================== 

#' * methods
#' @name *
#' @description Update \code{FLStock} and \code{FLIndex} objects with simulations from stock assessment fits
#' @rdname multiplication-methods
#' @aliases *,FLStock,a4aFitSA-method

setMethod("*", c("FLStock", "a4aFitSA"), function(e1, e2) 
{
  niters <- dims(e1) $ iter
  niters2 <- dim(pars(e2) @ stkmodel @ params)[2]
  if (niters > 1 & niters2 == 1) {
    nsim = niters
  } else {
    nsim = 1
    if (niters > niters2) stop("oh oh")
    if (niters == 1 & niters2 > 0) {
      niters <- niters2
      e1 <- propagate(e1, niters)
    }
  }

  # build up a4aFitSA to simulate from
  mod <- new("a4aFitSA")
  mod @ pars <- pars(e2)
  mod @ index <- index(e2)
  mod @ catch.n <- catch.n(e1)
  mod @ stock.n <- stock.n(e1)
  mod @ harvest <- harvest(e1)
  mod @ range <- range(e1)

  simstock <- simulate(mod, nsim = nsim)  

  catch.n(e1) <- catch.n(simstock)
  stock.n(e1) <- stock.n(simstock)
  harvest(e1) <- harvest(simstock)
    
  catch(e1) <- computeCatch(e1)
  stock(e1) <- computeStock(e1)
  
  e1
})

#' @rdname multiplication-methods
#' @aliases *,FLStock,SCAPars-method
setMethod("*", c("FLStock", "SCAPars"), function(e1, e2) 
{
stop("method * FLStock, SCAPars deprecated !\n")
})

#' @rdname multiplication-methods
#' @aliases *,FLIndices,a4aFitSA-method
setMethod("*", c("FLIndices", "a4aFitSA"), function(e1, e2) 
{
  e1 * pars(e2)
})

#' @rdname multiplication-methods
#' @aliases *,FLIndices,SCAPars-method
setMethod("*", c("FLIndices", "SCAPars"), function(e1, e2) 
{

  niters <- dims(e1) $ iter
  niters2 <- dim(e2 @ stkmodel @ params)[2]
  if (niters > 1 & niters2 == 1) {
    nsim = niters
  } else {
    nsim = 1
    if (niters > niters2) stop("oh oh")
    if (niters == 1 & niters2 > 0) {
      niters <- niters2
      e1 <- propagate(e1, niters)
    }
  }

  stop("not implemented yet")

  #FLIndices(out)
})

