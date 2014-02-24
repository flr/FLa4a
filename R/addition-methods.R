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

  mod <- new("a4aFitSA")
  mod @ pars <- e2

  simstock <- simulate(mod, nsim = nsim)  

  catch.n(e1) <- catch.n(simstock)
  stock.n(e1) <- stock.n(simstock)
  harvest(e1) <- harvest(simstock)
    
  catch(e1) <- computeCatch(e1)
  stock(e1) <- computeStock(e1)
  
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


  out <- lapply(seq_along(e1), 
  function(i) 
  {

    iind <- e1[[i]]
    
    years <- range(iind)[c("minyear","maxyear")]
    ages <- range(iind)[c("min","max")]

  
    #
    # Build design matrix for catches only
    #
    full.df <- expand.grid(age  = ages[1]:ages[2],
                           year = years[1]:years[2])[2:1]
 
  #  if (!is.null(covar)) {
  #  # add in covariates to data.frame - it is easiest to provide covariates in one list
  #  tmp <- 
  #    lapply(seq_along(covar), 
  #      function(i) {
  #        x <- as.data.frame(covar[[i]])[c(1,2,7)]
  #        if (length(unique(x $ age)) == 1) x <- x[names(x) != "age"]
  #        if (length(unique(x $ year)) == 1) x <- x[names(x) != "year"]
  #        names(x) <- gsub("data", names(covar)[i], names(x))
  #        x
  #      })
  #  covar.df <- tmp[[1]]
  #  for (i in seq(length(covar) - 1)) covar.df <- merge(covar.df, tmp[[i + 1]], all = TRUE, sort = FALSE)
  #
  #  full.df <- merge(full.df, covar.df, all.x = TRUE, all.y = FALSE)
  #  } 

    # make sure contrasts are set to sumto zero to match fit
    opts <- options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly")) 
  
    # f model matrix
    Xf <- Matrix(getX(e2 @ stkmodel @ fMod, full.df))

    # initial age structure model matrix
    Xny1 <- getX(e2 @ stkmodel @ n1Mod, subset(full.df, year == min(year) & age > min(age)))
    
    # now separate the sr model element
    facs <- strsplit(as.character(e2 @ stkmodel @ srMod)[length(e2 @ stkmodel @ srMod)], "[+]")[[1]]
    facs <- gsub("(^ )|( $)", "", facs) # remove leading and trailing spaces
    a4as <- grepl(paste("(^",c("bevholt", "ricker","hockey","geomean"),"[(])", collapse = "|", sep = ""), facs)

    # internal r model matrix
    if (sum(a4as) == 0) rmodel <- e2 @ stkmodel @ srMod else rmodel <- ~ factor(year) 
    Xr <- getX(rmodel, subset(full.df, age == min(age)))

    # reset options
    options(opts)

    # always simulate from b distribution for SA class.  If you want fitted values do FLStock + a4aFit(a4aFitSA)
    b.sim <- Matrix(simulate(e2, nsim = nsim) @ stkmodel @ params @ .Data)

    # matrix of predictions
    Xbeta <- bdiag(Xf, Xny1, Xr) %*% b.sim

    # plusgroup?
    plusgrp <- !is.na(range(e1)["plusgroup"]) && range(e1)["plusgroup"] >= range(e1)["max"]

    # unpack m - good for recycling
    Ms   <- c(m(e1) @ .Data)
 
    # build stock
    Fs <- Ns <- array(exp(Xbeta[1:nrow(Xf),]), dim = c(diff(ages)+1, diff(years)+1, ncol(Xbeta)))
    Ns[] <- NA
    Ns[-1,1,] <- array(exp(Xbeta[nrow(Xf) + 1:nrow(Xny1),]), dim = c(diff(ages), 1, ncol(Xbeta)))
    Ns[1,,] <- array(exp(Xbeta[nrow(Xf) + nrow(Xny1) + 1:nrow(Xr),]), dim = c(1, diff(years)+1, ncol(Xbeta)))
    Zs <- Fs + Ms
    for (a in 2:dim(Ns)[1]) {
      Ns[a,-1,] <- Ns[a-1, 1:diff(years),] * exp( - Zs[a-1, 1:diff(years),] )
    }
    # if plus group
    if (plusgrp) {
      for (y in 1:diff(years)) Ns[a,y+1,] <- Ns[a,y+1,] + Ns[a, y,] * exp( - Zs[a, y,] )
    } 
    # apply centering
    Ns <- Ns * exp(e2 @ stkmodel @ centering)
 
    zfrac <- Fs / Zs * (1 - exp(-Zs))

    dmns <- list(age    = paste(ages[1]:ages[2]), 
                 year   = paste(years[1]:years[2]),
                 unit   = "unique", 
                 season = "all", 
               area   = "unique", 
               iter = paste(1:niters))
               
    dms <- unname(c(dims(e1) $ age, dims(e1) $ year, 1, 1, 1, dims(e1) $ iter))

    stock.n(e1) <- FLQuant(Ns, dim = dms, dimnames = dmns, units = units(catch.n(e1)))
    catch.n(e1) <- FLQuant(zfrac * Ns, dim = dms, dimnames = dmns, units = units(catch.n(e1)))
    harvest(e1) <- FLQuant(Fs, dim = dms, dimnames = dmns, units = "f")
  
    catch(e1) <- computeCatch(e1)
    stock(e1) <- computeStock(e1)
  
    e1
  })

  FLIndices(out)
})

