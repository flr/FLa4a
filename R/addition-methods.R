#==================================================================== 
# "+" methods
#==================================================================== 

#' + methods
#' @name +
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
#' @aliases +,FLStock,a4aFitSA-method

setMethod("+", c("FLStock", "a4aFitSA"), function(e1, e2) 
{
  e1 + pars(e2)
})

#' @rdname addition-methods
#' @aliases +,FLStock,SCAPars-method
setMethod("+", c("FLStock", "SCAPars"), function(e1, e2) 
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

  # build up a4aFitSA to simulate from
  mod <- new("a4aFitSA")
  mod @ pars <- e2
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


#' @rdname addition-methods
#' @aliases +,FLIndices,a4aFitSA-method
setMethod("+", c("FLIndices", "a4aFitSA"), function(e1, e2) 
{
  e1 + pars(e2)
})


#' @rdname addition-methods
#' @aliases +,FLIndices,SCAPars-method
setMethod("+", c("FLIndices", "SCAPars"), function(e1, e2) 
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

