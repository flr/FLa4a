#==================================================================== 
#    predict  methods
#==================================================================== 

#' Predict methods for stock assessment fits
#' @name predict
#' @rdname predict-methods
#' @aliases predict,a4aFitSA-method
setMethod("predict", signature(object = "a4aFitSA"),
  function(object) {
    predict(pars(object))
  })


#' @rdname predict-methods
#' @aliases predict,SCAPars-method
setMethod("predict", signature(object = "SCAPars"),
  function(object) {
    list(
      stkmodel = predict(stkmodel(object)),
      qmodel   = predict(qmodel(object)),
      vmodel   = predict(vmodel(object))
    )
  })

#' @rdname predict-methods
#' @aliases predict,a4aStkParams-method
setMethod("predict", signature(object = "a4aStkParams"),
  function(object) {
      ages <- range(object)["min"]:range(object)["max"]
      years <- range(object)["minyear"]:range(object)["maxyear"]
      cnames <- rownames(coef(object))
      df <- expand.grid(age = ages,
                        year = years)
      X <- getX(object @ fMod, df)
      b <- coef(object)[grep("fMod", cnames)]
      niter <- dim(b)[2] # reuse this for the others
      fit <- exp(c(X %*% b))
      harvest <- 
        FLQuant(array(fit, dim = c(length(ages), length(years), 1, 1, 1, niter), 
                      dimnames = list(age = ages, year = years, 
                                         unit = "unique", season = "all", area = "unique",
                                         iter = seq(niter))),
                units = "f")


      X <- getX(object @ n1Mod, data.frame(age = ages[-1]))
      b <- coef(object)[grep("n1Mod", cnames)]
      fit <- exp(c(rbind(NA, X %*% b)) + object @ centering)
      ny1 <-       
        FLQuant(array(fit, dim = c(length(ages), 1, 1, 1, 1, niter), 
                      dimnames = list(age = ages, year = years[1], 
                                         unit = "unique", season = "all", area = "unique",
                                         iter = seq(niter))))


      X <- getX(object @ srMod, data.frame(year = years))
      b <- coef(object)[grep("rMod", cnames)]
      fit <- c(exp(c(X %*% b) + object @ centering))
      rec <-     
        FLQuant(array(fit, dim = c(1, length(years), 1, 1, 1, niter), 
                      dimnames = list(age = ages[1], year = years, 
                                         unit = "unique", season = "all", area = "unique",
                                         iter = seq(niter))))


      list(harvest = harvest, rec = rec, ny1 = ny1)
})

#' @rdname predict-methods
#' @aliases predict,submodels-method
setMethod("predict", signature(object = "submodels"),
  function(object, ...) {
      lapply(object, predict)
  })


#' @rdname predict-methods
#' @aliases predict,submodel-method
setMethod("predict", signature(object = "submodel"),
  function(object, ...) {
      ages <- range(object)["min"]:range(object)["max"]
      years <- range(object)["minyear"]:range(object)["maxyear"]
      df <- expand.grid(age = ages,
                        year = years)
      X <- getX(object @ Mod, df)
      b <- coef(object)
      niter <- dim(b)[2]
      fit <- exp(c(X %*% b) + object @ centering)
      FLQuant(array(fit, dim = c(length(ages), length(years), 1, 1, 1, niter), 
                         dimnames = list(age = ages, year = years, 
                                         unit = "unique", season = "all", area = "unique",
                                         iter = seq(niter))))
  })

