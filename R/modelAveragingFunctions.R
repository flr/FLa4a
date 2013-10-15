
setMethod("+", c("FLStock", "a4aFit"), function(e1, e2) 
{

  niters <- dims(e1) $ iter
  if (niters > 1) stop("adding a basic a4aFit object only makes sence with 1 iteration")

  years <- range(e1)[c("minyear","maxyear")]
  ages <- range(e1)[c("min","max")]

  dmns <- list(age    = paste(ages[1]:ages[2]), 
               year   = paste(years[1]:years[2]),
               unit   = "unique", 
               season = "all", 
               area   = "unique", 
               iter = paste(1:niters))
               
  dms <- unname(c(dims(e1) $ age, dims(e1) $ year, 1, 1, 1, dims(e1) $ iter))

  stock.n(e1) <- stock.n(e2)
  catch.n(e1) <- catch.n(e2)
  harvest(e1) <- harvest(e2)
  
  catch(e1) <- computeCatch(e1)
  stock(e1) <- computeStock(e1)
  
  e1
})


setMethod("+", c("FLStock", "a4aFitSA"), function(e1, e2) 
{

  niters <- dims(e1) $ iter

  years <- range(e1)[c("minyear","maxyear")]
  ages <- range(e1)[c("min","max")]

  
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
  Xf <- Matrix(getX(e2 @ pars @ stkmodel @ fMod, full.df))

  # initial age structure model matrix
  Xny1 <- getX(e2 @ pars @ stkmodel @ n1Mod, subset(full.df, year == min(year) & age > min(age)))
    
  # now separate the sr model element
  facs <- strsplit(as.character(e2 @ pars @ stkmodel @ srMod)[length(e2 @ pars @ stkmodel @ srMod)], "[+]")[[1]]
  facs <- gsub("(^ )|( $)", "", facs) # remove leading and trailing spaces
  a4as <- grepl(paste("(^",c("bevholt", "ricker","hockey","geomean"),"[(])", collapse = "|", sep = ""), facs)

  # internal r model matrix
  if (sum(a4as) == 0) rmodel <- e2 @ pars @ stkmodel @ srMod else rmodel <- ~ factor(year) 
  Xr <- getX(rmodel, subset(full.df, age == min(age)))

  # reset options
  options(opts)

  #
  # get the f, n1 and r model params and covariance matrix
  #
  which <- grepl("(fMod:)|(n1Mod:)|(rMod:)", dimnames(e2 @ pars @ stkmodel @ params) $ params)
  beta <- c(e2 @ pars @ stkmodel @ params[which])
  Vbeta <- e2 @ pars @ stkmodel @ vcov[which,which,1]
  X <- bdiag(Xf, Xny1, Xr)

  # always simulate from b distribution for SA class.  If you want fitted values do FLStock + a4aFit(a4aFitSA)
  if (niters == 1) {
    b.sim <- Matrix(mvrnorm(1, beta, Vbeta), ncol = 1)
  } else {
    b.sim <- Matrix(t(mvrnorm(niters, beta, Vbeta)))
  }
  
  # matrix of predictions
  Xbeta <- X %*% b.sim

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
  Ns <- Ns * exp(e2 @ pars @ stkmodel @ centering)
 
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


setMethod("mvrnorm", signature("numeric", "a4aFitSA"),
  function(n = 1, mu) 
  {
    #
    # get the f, n1 and r model params and covariance matrix
    #
    which <- grepl("(fMod:)|(n1Mod:)|(rMod:)", dimnames(mu @ pars @ stkmodel @ params) $ params)
    beta <- c(mu @ pars @ stkmodel @ params[which])
    Vbeta <- mu @ pars @ stkmodel @ vcov[which,which,1]


    # always simulate from b distribution for SA class.  If you want fitted values to FLStock + a4aFit(a4aFitSA)
    if (n == 1) {
      b.sim <- Matrix(mvrnorm(1, beta, Vbeta), ncol = 1)
      rownames(b.sim) <- colnames(Vbeta)
    } else {
      b.sim <- Matrix(t(mvrnorm(n, beta, Vbeta)))
    }
    
    b.sim
  }
)


setMethod("mvrnorm", signature("numeric", "a4aFitExt"),
  function(n = 1, mu, tol = 1e-9) 
  {
    # get the mean
    beta <- c(mu @ baseLvlPars)
    # number of params
    p <- length(beta)
    # simlulate iid normals (0,1)
    Z <- matrix(rnorm(p * n), p)
    # work out density using Z - the quick way
    ldens <-  -p/2 * log(2 * pi) + sum(log(diag(mu @ L))) - 0.5 * colSums(Z^2)
    # transform and rescale
    X <- drop(beta) + mu @ L %*% Z
    rownames(X) <- names(beta)
    attr(X, "logdens") <- ldens
    X
  }
)


#' @export
calcLogLik <- function(b, object, stock) # stock is required for M only
{

  if (is.null(dim(b))) b <- matrix(b, ncol = 1, dimnames = list(names(b)))

  #first split params up
  pnames <- names(object @ baseLvlPars)
  bf <- b[grep("fMod:", pnames),]
  bn1 <- b[grep("n1Mod:", pnames),]
  br <- b[grep("rMod:", pnames),]
  bq <- b[grep("qMod:", pnames),]
  bv <- b[grep("vMod:", pnames),]

  # now get the models?  Or do we just use the full design mat...
  pf <- object @ designMatrix $ Xf %*% bf
  pn1 <- object @ designMatrix $ Xn1 %*% bn1
  pr <- object @ designMatrix $ Xr %*% br
  pq <- object @ designMatrix $ Xq %*% bq
  pv <- object @ designMatrix $ Xv %*% bv

  # now convert to stock.n and calculate catch.n and index.n
  plusgrp <- !is.na(range(object)["plusgroup"]) && range(object)["plusgroup"] >= range(object)["max"]

  # unpack m - good for recycling
  Ms   <- c(m(stock) @ .Data)
 
  # build stock
  years <- range(object)[c("minyear","maxyear")]
  ages <- range(object)[c("min","max")]

  Fs <- Ns <- array(exp(pf), dim = c(diff(ages)+1, diff(years)+1, ncol(b)))
  Ns[] <- NA
  Ns[-1,1,] <- array(exp(pn1), dim = c(diff(ages), 1, ncol(b)))
  Ns[1,,] <- array(exp(pr), dim = c(1, diff(years)+1, ncol(b)))
  Zs <- Fs + Ms
  for (a in 2:dim(Ns)[1]) {
    Ns[a,-1,] <- Ns[a-1, 1:diff(years),] * exp( - Zs[a-1, 1:diff(years),] )
  }
  # if plus group
  if (plusgrp) {
    for (y in 1:diff(years)) Ns[a,y+1,] <- Ns[a,y+1,] + Ns[a, y,] * exp( - Zs[a, y,] )
  } 
  Ns <- Ns * exp(object @ pars @ stkmodel @ centering)
 
  # catch.n
  zfrac <- Fs / Zs * (1 - exp(-Zs))
  Cs <- zfrac * Ns

  # index.n
  nindices <- length(object @ pars @ qmodel)
  dim(pq) <- c(diff(ages)+1, diff(years)+1, nindices, ncol(b))
  dimnames(pq) <- list(paste(ages[1]:ages[2]), paste(years[1]:years[2]))
  Is <- exp(pq) * c(Ns[,,rep(1:ncol(b), each = nindices)])

  Is <- lapply(1:nindices, function(i) {
    rng <- range(object @ pars @ qmodel[[i]])
    surveyTime <- mean(rng[c("startf","endf")])
    out <- Is[,,i,,drop=FALSE] * exp(- surveyTime * c(Zs))
    out <- out * exp(object @ pars @ qmodel[[i]] @ centering - object @ pars @ stkmodel @ centering)
    out[paste(rng[1]:rng[2]),paste(rng[4]:rng[5]),,,drop=FALSE]
  })

  dim(pv) <- c(diff(ages)+1, diff(years)+1, nindices + 1, ncol(b))
  dimnames(pv) <- list(paste(ages[1]:ages[2]), paste(years[1]:years[2]))
  Vs <- exp(pv)

  Vs <- lapply(0:nindices, function(i) {
    if (i == 0) {
      Vs[,,1,,drop=FALSE]
    } else {
      rng <- range(object @ pars @ qmodel[[i]])
      Vs[paste(rng[1]:rng[2]),paste(rng[4]:rng[5]),i+1,,drop=FALSE]
    }
  })


  # calculate log likihood
  dat <- data.frame(obs = c(c(catch.n(stock)), unname(unlist(lapply(indices, index)))))
  sapply(1:ncol(b),
  function(i) {
    dat $ pred <- c(c(Cs[,,i]), unlist(lapply(Is, function(x) x[,,,i])))
    dat $ var <- unlist(lapply(Vs, function(x) x[,,,i]))
    llik <- with(dat[complete.cases(dat),], sum(dnorm(log(obs), log(pred), var, log = TRUE)))
    llik
  })
}

#' @export
calcLogPrior <- function(b, object)
{
  sum(dnorm(b, 0, 1e6, log = TRUE))  # flat priors on everthing for now!!  need to change to something better.
}








setMethod("getX", "a4aFitSA",
function(object) 
{

  e2 <- object

  ranges <- c(list(catch = range(e2)), lapply(e2 @ pars @ qmodel, slot, "range"))


  year <- range(e2)[c("minyear","maxyear")]
  age <- range(e2)[c("min","max")]

  #
  # Build design matrix
  #  
  make.df <- function(fleet) {
    expand.grid(age = ranges[[fleet]]["min"]:ranges[[fleet]]["max"], 
                year = ranges[[fleet]]["minyear"]:ranges[[fleet]]["maxyear"])[2:1]
  }
    
  full.df <- do.call(rbind, lapply(1:length(ranges), function(i) cbind(fleet = i, make.df(i))))
  
  # make sure contrasts are set to sumto zero for better performance
  opts <- options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly")) 

  # F model matrix
  Xf <- getX(e2 @ pars @ stkmodel @ fMod, subset(full.df, fleet == 1))
  # F model offsets
  # ...
  
  # Q model matrix  
  fleet.names <- c("catch", names(indices))
  Xqlist <- lapply(seq_along(indices), function(i) getX(e2 @ pars @ qmodel[[i]] @ Mod, subset(full.df, fleet == i + 1)))
  Xq <- as.matrix(do.call(bdiag, Xqlist))  

  
  # var model matrix
  Xvlist <- lapply(1:length(fleet.names), function(i) getX(e2 @ pars @ vmodel[[i]] @ Mod, subset(full.df, fleet == i)))
  Xv <- as.matrix(do.call(bdiag, Xvlist))   

  # initial age structure model matrix
  Xny1 <- getX(e2 @ pars @ stkmodel @ n1Mod, subset(full.df, year == min(year) & age > min(age) & fleet == 1))
  
  # now separate model elements, find SR models, stop if more than one specified
  srmodel <- e2 @ pars @ stkmodel @ srMod
  facs <- strsplit(as.character(srmodel)[length(srmodel)], "[+]")[[1]]
  facs <- gsub("(^ )|( $)", "", facs) # remove leading and trailing spaces
  a4as <- grepl(paste("(^",c("bevholt", "ricker","hockey","geomean"),"[(])", collapse = "|", sep = ""), facs)

  # internal r model matrix
  if (sum(a4as) == 0) rmodel <- srmodel else rmodel <- ~ factor(year) 
  Xr <- getX(rmodel, subset(full.df, age == min(age) & fleet == 1))

 
  # reset contrast options
  options(opts)
  
  X <- bdiag(Xf, Xny1, Xr, Xv, Xq)


  attr(X, "composition") <- list(nf   = nrow(Xf),
                                 nny1 = nrow(Xny1),
                                 nr   = nrow(Xr),
                                 nq   = nrow(Xq),
                                 nv   = nrow(Xv))

  X
})



#' @export
makePrediction <- function(b.sim, X, e1, e2) {
  
  
  years <- range(e2)[c("minyear","maxyear")]
  ages <- range(e2)[c("min","max")]
  comps <- attr(X, "composition")

  
  # matrix of predictions
  Xbeta <- X %*% b.sim
  niters <- ncol(b.sim)

  # plusgroup?
  plusgrp <- !is.na(range(e2)["plusgroup"]) && range(e2)["plusgroup"] >= range(e2)["max"]

  # unpack m - good for recycling
  Ms   <- c(m(e1) @ .Data)
 
  # build stock
  Fs <- Ns <- array(exp(Xbeta[1:comps$nf,]), dim = c(diff(ages)+1, diff(years)+1, ncol(Xbeta)))
  Ns[] <- NA
  Ns[-1,1,] <- array(exp(Xbeta[comps$nf + 1:comps$nny1,]), dim = c(diff(ages), 1, ncol(Xbeta)))
  Ns[1,,] <- array(exp(Xbeta[comps$nf + comps$nny1 + 1:comps$nr,]), dim = c(1, diff(years)+1, ncol(Xbeta)))
  Zs <- Fs + Ms
  for (a in 2:dim(Ns)[1]) {
    Ns[a,-1,] <- Ns[a-1, 1:diff(years),] * exp( - Zs[a-1, 1:diff(years),] )
  }
  # if plus group
  if (plusgrp) Ns[a,-1,] <- Ns[a,-1,] + Ns[a, 1:diff(years),] * exp( - Zs[a, 1:diff(years),] )
 
  # apply centering
  Ns <- Ns * exp(e2 @ pars @ stkmodel @ centering)
 
  zfrac <- Fs / Zs * (1 - exp(-Zs))

  dmns <- list(age    = paste(ages[1]:ages[2]), 
               year   = paste(years[1]:years[2]),
               unit   = "unique", 
               season = "all", 
               area   = "unique", 
               iter = paste(1:niters))
               
  dms <- unname(c(dims(e1) $ age, dims(e1) $ year, 1, 1, 1, dims(e1) $ iter))

  FLQuant(zfrac * Ns, dim = dms, dimnames = dmns, units = units(catch.n(ple4)))
}


