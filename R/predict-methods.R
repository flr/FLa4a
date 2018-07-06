#====================================================================
#    predict methods
#====================================================================

#' @title Predict methods for SCA
#' @name predict for sca
#' @description Predict methods for a4a stock assessment fits.
#' @aliases predict,a4aFitSA-method
#' @rdname predict-methods
#' @template object
#' @examples
#' data(ple4)
#' data(ple4.index)
#' fmodel <- ~factor(age) + factor(year)
#' qmodel <- list(~factor(age))
#' fit1 <-  sca(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#' flqs <- predict(fit1)
setMethod("predict", signature(object = "a4aFitSA"),
  function(object) {
    obj <- pars(object)
  # need to tag biomass indices
    for(i in seq_along(index(object)))
      attr(obj@vmodel[[i+1]], "FLIndexBiomass") <- attr(obj@qmodel[[i]], "FLIndexBiomass") <- attr(index(object)[[i]], "FLIndexBiomass")
    predict(obj)
  })

#' @rdname predict-methods
setMethod("predict", signature(object = "SCAPars"),
  function(object) {

  sm <- stkmodel(object)
  qm <- qmodel(object)
  vm <- vmodel(object)
  # run predict
    lst <- list(
      stkmodel = predict.stkpars(sm),
      qmodel   = predict.submods(qm, type = "response"),
      vmodel   = predict.submods(vm, type = "response")
    )
  # rec is being estimated from the srmodel so the n1model doesn't
  # have a recruitment estimate, need to update from rec predictions
  lst$stkmodel$ny1[1,1] <- lst$stkmodel$rec[,1]
  lst
})

predict.stkpars <- function(object) {
  ages <- range(object)["min"]:range(object)["max"]
  years <- range(object)["minyear"]:range(object)["maxyear"]
  cnames <- rownames(coef(object))
  df <- expand.grid(age = ages,
                    year = years)
  niter <- dim(coef(object))[2] # reuse this for the others
  if (dim(object@centering)[2] == 1) {
    object@centering <- propagate(object@centering, niter)
  } # otherwise rely on propagates error message

  # predict F
  X <- getX(fMod(object), df)
  b <- coef(object)[grep("fMod", cnames)]
  fit <- object@linkinv(X %*% b)
  harvest <-
    FLQuant(array(fit, dim = c(length(ages), length(years), 1, 1, 1, niter),
                  dimnames = list(age = ages, year = years,
                                  unit = "unique", season = "all", area = "unique",
                                  iter = seq(niter))),
            units = "f")

  # predict N in first year
  X <- getX(n1Mod(object), data.frame(age = ages[-1]))
  b <- coef(object)[grep("n1Mod", cnames)]
  fit <- rbind(NA, object@linkinv(sweep(X %*% b, 2, object @ centering, "+")))
  ny1 <-
    FLQuant(array(fit, dim = c(length(ages), 1, 1, 1, 1, niter),
                  dimnames = list(age = ages, year = years[1],
                                  unit = "unique", season = "all", area = "unique",
                                  iter = seq(niter))))

  # predict Recruitment - no need to worry about the SR model - it is included in the fit :)
  rMod <- if (isPresenta4aSRmodel(srMod(object))) ~ factor(year) else srMod(object)
  X <- getX(rMod, data.frame(year = years))
  b <- coef(object)[grep("rMod", cnames)]
  fit <- object@linkinv(sweep(X %*% b, 2, object @ centering, "+"))
  rec <- FLQuant(array(fit, dim = c(1, length(years), 1, 1, 1, niter),
                       dimnames = list(age = ages[1], year = years,
                                       unit = "unique", season = "all", area = "unique",
                                       iter = seq(niter))))

  FLQuants(harvest = harvest, rec = rec, ny1 = ny1)
}


predict.submods <- function(object, ...) {
  genFLQuant(object, ...)
}

predict.submod <- function(object, ...) {
  genFLQuant(object, ...)
}
