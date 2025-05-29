###############################################################################
# EJ&CM(2012821)
# Methods to add uncertainty to FLQuant objects
# NOTE #1:
# NOTE #2:
###############################################################################

#' Methods to generate FLStock objects
#' @description This method computes the \code{FLStock} slots consistently with the information provided by the \code{FLQuant}. It requires two of the triplet R/C/F to compute the third consistent with Baranov and survival's equations.
#' @param object an FLStock
#' @param R an FLQuant with iterations or missing
#' @param C an FLQuant with iterations or missing
#' @param F an FLQuant with iterations or missing
#' @param ... additional argument list that might not ever be used.
#' @return an FLStock
#' @docType methods
#' @rdname genFLStock-methods
#' @aliases genFLStock genFLStock-methods

setGeneric("genFLStock", function(object, R, C, F, ...) standardGeneric("genFLStock"))

#' @rdname genFLStock-methods
setMethod("genFLStock", c("FLStock", "FLQuant", "FLQuant", "missing"), function(object, R, C, F, ...){
  cat("Not implemented yet\n")
})

#' @rdname genFLStock-methods
setMethod("genFLStock", c("FLStock", "missing", "FLQuant", "FLQuant"), function(object, R, C, F, ...){
  cat("Not implemented yet\n")
})

#' @rdname genFLStock-methods
setMethod("genFLStock", c("FLStock", "FLQuant", "missing", "FLQuant"), function(object, R, C, F, ...){
  args <- list(...)

  # requires checking dimensions
  if(!identical(dim(catch.n(object))[-c(1,6)], dim(R)[-c(1,6)])) stop("Recruitment vector must have consistent dimensions with the stock object")
  if(!identical(dim(catch.n(object))[-6]    , dim(F)[-6])) stop("Harvest matrix must have consistent dimensions with the stock object")
  if(dim(R)[6]!=dim(R)[6]) stop("R and F must have the same number of iterations")

  # get dims and set flq
  dms <- dims(object)
  nages <- dms$age
  nyrs <- dms$year
  niters <- dim(R)[6]
  minyr <- range(object)["minyear"]
  flq <- FLQuant(dimnames=dimnames(F))
  # create id matrix for cohorts
  flcid <- FLCohort(F)
  flcid[!is.na(flcid)] <- 1

  # compute cumulative Z
  Z <- F + m(object)
  Z <- FLCohort(Z)
  Z[is.na(Z)] <- 0
  Z[] <- apply(Z, c(2:6), cumsum)

  # expand variability into [N] by R*[F]
  #ifelse(sum(tolower(names(args))=="ny1")==1, ny1 <- args$ny1, ny1 <- stock.n(object)[,1])
  ny1 <- args$ny1
  Ns <- flq
  Ns[1] <- R
  Ns <- FLCohort(Ns)
  Ns <- Ns[rep(1, nages)]

  for(i in 2:nages){
    Ns[,nages-i+1] <- ny1[i]
  }

  # survivors
  Ns <- Ns*exp(-Z)
  Ns <- as(Ns, "FLQuant")

  # Update object
  stock.n(object) <- flq
  # ny1
  stock.n(object)[,1] <- ny1
  # [R]
  stock.n(object)[1] <- R
  # [N]
  stock.n(object)[-1,-1] <- Ns[-nages,-nyrs]
  # plus group
  stock.n(object)[nages,-1] <- Ns[nages,-nyrs] + stock.n(object)[nages,-1]
  # catch
  stock(object) <- computeStock(object)
  # [F]
  harvest(object) <- F
  # [C]
  Z <- harvest(object) + m(object)
  Cs <- harvest(object)/Z*(1-exp(-Z))*stock.n(object)
  catch.n(object) <- Cs
  catch(object) <- computeCatch(object)
  # [L] & [D] rebuilt from C
  # Ds=D/(D+L)*Cs where Cs is the simulated catch
  D <- discards.n(object)
  L <- landings.n(object)
  discards.n(object) <- D/(D+L)*Cs
  discards(object) <- computeDiscards(object)
  landings.n(object) <- Cs - discards.n(object)
  landings(object) <- computeLandings(object)
  # out
  object
})

#' @name getAcor
#' @rdname getAcor-methods
#' @title compute log-correlation matrix
#' @description Method to compute the log-correlation matrix for the first dimension (\code{quant}) of the \code{FLQuant} object.
#' @template bothargs
#' @return an \code{FLQuant} object with a quant log-correlation matrix
#' @aliases getAcor getAcor-methods
#' @examples
#' data(ple4)
#' getAcor(harvest(ple4))

setGeneric("getAcor", function(object, ...) standardGeneric("getAcor"))
#' @rdname getAcor-methods
setMethod("getAcor", c("FLQuant"), function(object, ...) {
    mu <- log(object)
    Rho <- cor(t(mu[drop = TRUE]))
    return(Rho)
})

#' @title Methods to generate FLQuant objects
#' @description This method uses the quant log-correlation matrix of the \code{FLQuant} object and generates a new \code{FLQuant} using a lognormal multivariate distribution.
#' @param object an FLQuant
#' @param cv the coefficient of variation
#' @param method the method used to compute the correlation matrix; for now only "ac" (autocorrelation) is implemented
#' @param niter the number of iterations to be generated
#' @param ... additional argument list that might not ever be used.
#' @return an FLQuant
#' @docType methods
#' @rdname genFLQuant-methods
#' @aliases genFLQuant genFLQuant-methods
#' @examples
#' data(ple4)
#' sim.F <- genFLQuant(harvest(ple4))
setGeneric("genFLQuant", function(object, ...) standardGeneric("genFLQuant"))

#' @rdname genFLQuant-methods
setMethod("genFLQuant", "FLQuant",
  function(object, cv = 0.2, method = "ac", niter = 250) {
  # use log transform, to be expanded on later versions
  mu <- log(object)
  if(method == "ac") {
    Rho <- cor(t(mu[drop = TRUE]))
    flq <- mvrnorm(niter * dim(mu)[2], rep(0, nrow(Rho)), log(cv^2+1) * Rho)
    mu <- propagate(mu, niter)
    flq <- FLQuant(c(t(flq)), dimnames = dimnames(mu))
    flq <- exp(mu + flq)
  }
  units(flq) <- units(object)
  return(flq)
})

#' @rdname genFLQuant-methods
#' @param type the type of output required. The default is on the scale of the linear predictors (link); 
#'             the alternative "response" is on the scale of the response variable. 
#'             Thus for a model on the log scale the default predictions are of log F (for example) 
#'             and type = "response" gives the predicted F. 
#' @param nsim the number of iterations to simulate, if nsim = 0, then deterministic values are returned
#'             based on the coefficients.  If nsim > 0 then coefficients are simluated using the
#'             covariance slot and distribution slot.
#' @param seed if supplied the random numbers are generate with a fixed seed for repeatablility
#' @param simulate.recruitment if FALSE (default) recruitment is simulated from
#'             the recruitment estimates of recruitment, which may or may not be based on a stock-recruit
#'             model in the origional fit.  If TRUE, then new recruitments are simulated based on the 
#'             stock recruitment model and supplied CV used in the fit, rsulting in a completly different
#'             timeseries of N and Catches.
# if nsim > 0 the simulate nsim times
setMethod("genFLQuant", "submodel",
  function(object, type = c("link", "response"), nsim = 0, seed = NULL) {
      type <- match.arg(type)
      # simulate from submodel?
      if (nsim > 0) {
        object <- simulate(object, nsim = nsim, seed = seed)
      }
      niter <- dim(coef(object))[2]
      # are there iters in centering?
      if (dim(object@centering)[2] == 1) {
        object@centering <- propagate(object@centering, niter)
      } # otherwise rely on propagates error message

      # make empty FLQuant
      flq <- flqFromRange(object)
      df <- as.data.frame(flq)
      # this should have 2 dimensions!
      b <- coef(object)
      # get design matrix
      X <- getX(formula(object), df)
      # predict accross all iters (if dimensions don't match then coefs are the wrong length!)
      pred <- sweep(X %*% as(b, "matrix"), 2, object@centering, "+")
      # add into flq
      flq <- propagate(flq, dims(b)$iter)
      flq[] <- as(pred, "vector")
      # transform if asked
      if (type == "response") {
        object@linkinv(flq)
      } else {
        flq
      }
    }
  )

#' @rdname genFLQuant-methods
# if nsim > 0 the simulate nsim times
setMethod("genFLQuant", "submodels",
  function(object, type = c("link", "response"), nsim = 0, seed = NULL) {
      type <- match.arg(type)
      # simulate from submodels?
      if (nsim > 0) {
        object <- simulate(object, nsim = nsim, seed = seed)
      }
      # convert submodels to FLQuants
      FLQuants(lapply(object, genFLQuant, type = type))
    }
 )



#' @rdname genFLQuant-methods
# if nsim > 0 the simulate nsim times
setMethod("genFLQuant", "a4aStkParams",
  function(object, type = c("link", "response"), nsim = 0, seed = NULL, simulate.recruitment = FALSE) {
    # type is always response for a stock model...
    type <- match.arg(type)
    # simulate from submodels?
    if (nsim > 0) {
      object <- simulate(object, nsim = nsim, seed = seed)
    }
    # predict F, rec, and ny1
    flqs <- predict.stkpars(object)
    # get dims
    dms <- dims(flqs$harvest)

    # Do we want to simulate from the stock recruitment model?
    if (simulate.recruitment) {
      # simulate N based on new recruitments about the estimated SR model
      # this will results in catches and survey indices quite different form that observed

      # get SR model
      srmodel <- geta4aSRmodel(srMod(object))

      # get FLSR definition
      expr_model <- a4aSRmodelDefinitions(srmodel)

      if (is.null(expr_model)) {
        stop("Cannot simulate recruitment from SR relationship as no stock-recruitment model was used")
      }

      # get SR pars
      cnames <- rownames(coef(object))
      parList <- list(a = exp(coef(object)[grep("sraMod", cnames)]), 
                      b = exp(coef(object)[grep("srbMod", cnames)]))
      srrCV <- eval(parse(text = srmodel))$srrCV

      # define the SR prediction
      recPred <- function(ssb) {
        # if ssb is an FLQuant, the return will be an FLQuant
        # stdev = sqrt(cv^2 + 1)
        ssb/ssb * eval(expr_model, c(parList, ssb = ssb)) * 
          exp(rnorm(dms$iter, 0, sqrt(log(srrCV^2 + 1)))) * # random noise
          exp(object@centering)
      }

      # set up N quant
      N <- flqs$harvest
      N[] <- NA
      units(N) <- "1000"

      # initial age structure
      N[1,1] <- flqs$rec[1,1]
      N[-1,1] <- flqs$ny1[-1,]

      # ssb per individual by age and year
      ssbay <- mat(object) * wt(object)
      Z <- flqs$harvest + m(object)

      for (i in 2:dms$year) {
        # predict recruitment
        N[1,i] <- recPred(quantSums(N[,i-1] * ssbay[,i-1]))

        # do some killing
        Nleft <- N[,i-1] * exp(-Z[,i-1])
        N[-1,i] <- Nleft[-dms$age]
        N[dms$age,i] <- N[dms$age,i] + Nleft[dms$age]
        # repeat!
      }

      # a quick debugging check - all looks good
      # plot(quantSums(N * ssbay)[,-dms$year], flqs$rec[,-1], ylim = c(0, max(flqs$rec)))
      # points(quantSums(N * ssbay)[,-dms$year], N[1,-1], col = "red")

    } else {
      # simulate N conditional on the estimated recruitment

      # compute cumulative Z
      cumsumNA <- function(x) {
        x[!is.na(x)] <- cumsum(x[!is.na(x)])
        x
      }
      Z <- flqs$harvest + m(object)

      cumZ <- apply(FLCohort(Z), c(2:6), cumsumNA)

      # expand variability into [N] by R*[F]

      Ns <- FLCohort(flqs$harvest)
      Ns[,-(1:(dms$age-1))] <- flqs$rec[rep(1,dms$age)]
      
      Ns[,1:(dms$age-1)] <- apply(FLCohort(flqs$ny1), 2:6, sum, na.rm = TRUE)[rep(1,dms$age),1:(dms$age-1)]
      
      Ns <- Ns * exp(-cumZ)
      units(Ns) <- object@units
      
      # convert back from cohort shape
      Ns <- as(Ns, "FLQuant")

      # add in recruits and initial age
      N <- Ns
      # [R]
      N[1] <- flqs$rec
      # [N]
      N[-1,-1] <- Ns[-dms$age,-dms$year]
      # [N,1]
      N[-1, 1] <- flqs$ny1[-1,]
      # plus group
      for(y in seq(2, dms$year))
        N[dms$age, y] <- Ns[dms$age - 1, y-1] + N[dms$age, y-1] * exp(-Z[dms$age, y-1])
    }
    # [C]
    Z <- flqs$harvest + m(object)
    C <- flqs$harvest/Z*(1-exp(-Z))*N
    units(C) <- units(N)

    # out
    if (type == "response") {
      FLQuants(harvest = flqs$harvest, stock.n = N, catch.n = C)  
    } else {
      FLQuants(harvest = object@link(flqs$harvest), stock.n = object@link(N), catch.n = object@link(C))  
    }
    
  }
)




#' Methods to generate FLIndex objects
#' @description This method produces an \code{FLIndex} object by using the \code{genFLQuant} method.
#' @param object an \code{FLIndex} object
#' @param cv the coefficient of variation
#' @param niter the number of iterations to be generated
#' @param ... additional argument list that might not ever be used.
#' @return an FLIndex
#' @docType methods
#' @rdname genFLIndex-methods
#' @aliases genFLIndex genFLIndex-methods

setGeneric("genFLIndex", function(object, ...) standardGeneric("genFLIndex"))

#' @rdname genFLIndex-methods
setMethod("genFLIndex", c("FLQuant"), function(object, cv = 0.2, niter = 250) {
      # use log transform, to be expanded on later versions
      mu <- log(object)

#      if(method == "ac") {
        Rho <- cor(t(mu[drop = TRUE]))
        flq <- mvrnorm(niter * dim(mu)[2], rep(0, nrow(Rho)), log(cv^2+1) * Rho)
        mu <- propagate(mu, niter)
        flq <- FLQuant(c(t(flq)), dimnames = dimnames(mu))
        flq <- exp(mu + flq)
#      }
      return(flq)
})

