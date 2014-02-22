






#' Hello 
#' @name logLik
#' @docType methods
#' @rdname logLik-methods
#' @aliases logLik,FLa4aFit-method
setMethod("logLik", signature(object = "a4aFit"),
  function(object, ...) 
  {  
    dim2 <- length(dim(object @ fitSumm))
    if (dim2 == 1) {
      val <- -1 * unname(object @ fitSumm["nlogl"])
      attr(val, "nobs") <- unname(object @ fitSumm["nobs"])
      attr(val, "df") <- unname(object@fitSumm["nopar"])
    } else if (dim2 == 2) {
      val <- -1 * unname(object @ fitSumm["nlogl",])
      attr(val, "nobs") <- unname(object @ fitSumm["nobs",])
      attr(val, "df") <- unname(object@fitSumm["nopar",])
    }
    class(val) <- "logLik"
    val
 })

#' Hello 
#' @name show
#' @docType methods
#' @rdname show-methods
#' @aliases show,FLa4aFit-method
setMethod("plot", signature(x = "a4aFit", y = "FLStock"),
  function (x, y, ratio = 1.5, file = "", onefile = TRUE, 
            #what = c("N","F","Q","Res"),
            what = c("N", "F", "Res"), 
            Ftext = TRUE, ask = TRUE, ...) 
  {

  op <- par(ask = ask, no.readonly = TRUE)

  ages  <- as.numeric(dimnames(stock.n(x)) $ age)
  years <- as.numeric(dimnames(stock.n(x)) $ year)

  cols <- 
     c(rgb(215, 48, 39, max = 255), 
       rgb(252, 141, 89, max = 255),
       rgb(254, 224, 144, max = 255),
       rgb(224, 243, 248, max = 255),
       rgb(145, 191, 219, max = 255),
       rgb(69, 117, 180, max = 255))

  zeros <- function(n) paste(rev(rep(c("0", "0", " 0"), length = n)), collapse = "")
  age.col <- colorRampPalette(cols[c(1, 6)])(length(ages))
      
  if ("N" %in% what) {  
    # N plot
    scale <- ceiling(max(log(stock.n(x)[drop=TRUE], 10))) - 2
    
    matplot(years, t(stock.n(x)[drop=TRUE]) * 10^{-scale},  
            ylab =paste0("N at age ('",zeros(scale),"s)"), xlab = "Year", 
            col = age.col,
            type ='l', lty = 1, lwd = 2, las = 1, main = "N-at-age in the Stock")
    lx = max(years) + diff(range(years))*0.04
    ly = (max(stock.n(x)[drop=TRUE]) + diff(range(stock.n(x)[drop=TRUE]))*0.04 ) * 10^{-scale}
    legend(lx, ly, legend = paste("age", ages), col = age.col, lty = 1, lwd = 2, xjust = 1, yjust = 1)

  }

  if ("F" %in% what) {  

    fest <- harvest(x)[drop=TRUE]

    # F plot
    matplot(years, t(fest),  
            ylab = "F at age", xlab = "Year", 
            col = colorRampPalette(cols[c(1, 6)])(length(ages)),
            type = 'l', lty = 1, lwd = 2, las = 1, main = "F-at-age")
    lx = max(years) + diff(range(years))*0.04
    ly = max(fest) + diff(range(fest))*0.04
    legend(lx, ly, legend = paste("age", ages), col = age.col, lty = 1, lwd = 2, xjust = 1, yjust = 1)


    # another F plot
    p <- matrix.plot(fest, cols = cols, xlab = "Year", ylab = "Age", main = "F-at-age", ymin = min(ages), xmin = min(years), text = Ftext)
    print(p)

    # Yet another F plot
    p <- wireframe(fest, drape = TRUE, colorkey = FALSE,
             screen = list(z = 240, x = -60), 
             col.regions = colorRampPalette(rev(cols))(100),
             panel.aspect = 1/ratio, aspect=c(length(years)/length(ages), 2),
             ylab = "Year", xlab = "Age", zlab = "F", main = "F-at-age",
             #par.settings = list(axis.line = list(col = "transparent")),
             par.box = c(col = "transparent"))
    print(p)

  }

#    sub.dev.new()
    # fbar plot
#    plotError(years, fbar(y)[,,,,,"mean",drop=TRUE], sqrt(fbar(x)[,,,,,"var",drop=TRUE]), 
#              ylab = 'Fbar', xlab = "Year", cols = colorRampPalette(cols[5:6])(3))

#    sub.dev.new()
    # ssb plot
#    plotError(years, ssb(y)[,,,,,"mean",drop=TRUE] * 1e-3, sqrt(ssb(x)[,,,,,"var",drop=TRUE]) * 1e-3, 
#              ylab = "SSB ('000 tonnes)", xlab = "Year", cols = colorRampPalette(cols[5:6])(3))

  if ("Q" %in% what) {  

    # q plots
    for (i in 1:length(x @ logq)) {
      qest <- logq(x)[[i]]
      qages  <- as.numeric(dimnames(qest) $ age)
      qyears <- as.numeric(dimnames(qest) $ year)

      matplot(qages, qest[drop=TRUE],  
              ylab="log catchability", xlab="Age", main = names(x @ logq)[i],
              col = colorRampPalette(cols[c(1, 6)])(length(qyears)),
              type ='l', lty = 1, lwd = 2, las = 1)

      p <- wireframe(qest[drop=TRUE], drape = TRUE, colorkey = FALSE,
             # screen value tuned to ple4.indices[1:2]
             screen = list(z = ifelse(i==1, -1, 1) * 30, x = -60), 
             col.regions = colorRampPalette(rev(cols))(100),
             panel.aspect = 1/ratio, ylab = "Year", xlab = "Age", 
             zlab = "log catchability", main = paste("log catchability-at-age of", names(x @ logq)[i]),
             par.box = c(col = "transparent"))
      print(p) 
    }

    # q plots
    rec.age <- range(y)["min"]
    ssb.x <- ssb(y)[,seq(rec.age) - length(ssb(y)) - 1][drop=TRUE]
    rec.y <- rec(y)[,-seq(rec.age)][drop=TRUE]
    plot(ssb.x, rec.y, pch = 19, ann = FALSE, xlim = c(0, max(ssb.x)), ylim = c(0, max(rec.y)))
    title(main = "Stock and Recruitment", xlab = "SSB", ylab = paste("Recruitment (age ", rec.age, ")", sep = ""))
    # TODO  add in stock recruit fit if a model was used.
  }

  if ("Res" %in% what) {  
    
    par(oma = c(4,0,0,0))
    
    # Residual plot 
    div <- list(c(1,1), c(2,1), c(3,1), c(2,2), c(3,2), c(3,2), c(3,3), c(3,3), c(3,3), c(4,3), c(4,3), c(4,3))
    nind <- length(x @ index)
    op2 <- par(mfrow = div[[nind + 1]], mar = c(4,4,1,2), mgp = c(2,1,0), no.readonly = TRUE)
    ylim <- range(as.numeric(dimnames(stock.n(x)) $ age)) + c(-.5, .5)
    xlim <- range(as.numeric(dimnames(stock.n(x)) $ year)) + c(-.5, .5)

    res <- log(catch.n(y) / catch.n(x))
    ages  <- as.numeric(dimnames(res) $ age)
    years <- as.numeric(dimnames(res) $ year)

    bp(rep(years, each = length(ages)), rep(ages, length(years)), 
       c(res[drop=TRUE]), 
       ylim = ylim, xlim = xlim, 
       xlab = 'year', ylab = 'Age', main = "Catch", 
       scale = 3, las = 1)

    

#    for(i in 1:nind) {
#      res <- index.lres(x)[[i]]
#      ages  <- as.numeric(dimnames(res) $ age)
#      years <- as.numeric(dimnames(res) $ year)

#      bp(rep(years, each = length(ages)), rep(ages, length(years)), 
#         c(res[drop=TRUE]), 
#         ylim = ylim, xlim = xlim, 
#         xlab = 'year', ylab = 'Age', main = names(x @ logq)[i], 
#         scale = 3, las = 1)
#    }

    mtext("Standardised residuals: Positive (red) and negative (blue).  \n50% of dots should be light coloured, 95% should be light and medium, \n5% of dots should be dark.", 
          side = 1, line = 6, adj = 0)
    par(op2)
  }
  
  par(op)
 
})


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
setGeneric("logq", function(object, ...) standardGeneric("logq"))

#' @rdname logq-methods
#' @aliases logq,FLa4aFit-method
setMethod("logq", signature(object = "a4aFit"),
  function(object) {
	object @ logq
  })


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
#' @rdname catch.hat-methods
setGeneric("catch.hat", function(object, ...) standardGeneric("catch.hat"))


#' Hello 
#' @name catch.hat
#' @docType methods
#' @rdname catch.hat-methods
#' @aliases catch.hat,FLa4aFit-method
setMethod("catch.hat", signature(object = "a4aFit"),
  function (object, ...) 
  {
    exp(object @ catch.lhat)
  })


q#' Calculate the median accross iterations
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
#' @rdname index.hat-methods
setGeneric("index.hat", function(object, ...) standardGeneric("index.hat"))


#' Hello 
#' @name index.hat
#' @docType methods
#' @rdname index.hat-methods
#' @aliases index.hat,FLa4aFit-method
setMethod("index.hat", signature(object = "a4aFit"),
  function (object, ...) 
  {
    lapply(object @ index.lhat, exp)
  })


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
#' @rdname catch.lhat-methods
setGeneric("catch.lhat", function(object, ...) standardGeneric("catch.lhat"))


#' Hello 
#' @name catch.lhat
#' @docType methods
#' @rdname catch.lhat-methods
#' @aliases catch.lhat,FLa4aFit-method
setMethod("catch.lhat", signature(object = "a4aFit"),
  function (object, ...) 
  {
    object @ catch.lhat
  })


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
#' @rdname index.lhat-methods
setGeneric("index.lhat", function(object, ...) standardGeneric("index.lhat"))


#' Hello 
#' @name index.lhat
#' @docType methods
#' @rdname index.lhat-methods
#' @aliases index.lhat,FLa4aFit-method
setMethod("index.lhat", signature(object = "a4aFit"),
  function (object, ...) 
  {
    object @ index.lhat
  })





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
#' @rdname catch.lvar-methods
setGeneric("catch.lvar", function(object, ...) standardGeneric("catch.lvar"))


#' Hello 
#' @name catch.lvar
#' @docType methods
#' @rdname catch.lvar-methods
#' @aliases catch.lvar,FLa4aFit-method
setMethod("catch.lvar", signature(object = "a4aFit"),
  function (object, ...) 
  {
    object @ catch.lvar
  })


q#' Calculate the median accross iterations
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
#' @rdname index.lvar-methods
setGeneric("index.lvar", function(object, ...) standardGeneric("index.lvar"))


#' Hello 
#' @name index.lvar
#' @docType methods
#' @rdname index.lvar-methods
#' @aliases index.lvar,FLa4aFit-method
setMethod("index.lvar", signature(object = "a4aFit"),
  function (object, ...) 
  {
    object @ index.lvar
  })


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
#' @rdname index.lres-methods
setGeneric("index.lres", function(object, ...) standardGeneric("index.lres"))


#' Hello 
#' @name index.lres
#' @docType methods
#' @rdname index.lres-methods
#' @aliases index.lres,FLa4aFit-method
setMethod("index.lres", signature(object = "a4aFit"),
  function (object, type = "scaled", ...) 
  {
    type <- match.arg(type, c("scaled", "unscaled"))
    if (type == "scaled") {
      object @ index.lres
    } else if (type == "unscaled") {
      out <- object @ index.lres
      for (i in 1:length(out)) {
        out[[i]] <- out[[i]] * sqrt(object @ index.lvar[[i]])
      }
      out
    }
  })


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
#' @rdname catch.lres-methods
setGeneric("catch.lres", function(object, ...) standardGeneric("catch.lres"))


#' Hello 
#' @name catch.lres
#' @docType methods
#' @rdname catch.lres-methods
#' @aliases catch.lres,FLa4aFit-method
setMethod("catch.lres", signature(object = "a4aFit"),
  function (object, type = "scaled", ...) 
  {
    type <- match.arg(type, c("scaled", "unscaled"))
    if (type == "scaled") {
      object @ catch.lres
    } else if (type == "unscaled") {
      object @ catch.lres * sqrt(object @ catch.lvar)
    }
  })


# -------------------------------------------------------------------
#
#
#    coef  methods
#
#
# -------------------------------------------------------------------


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
#' @rdname coef-methods
#'
#' @examples
#' data(ple4)
#' genFLQuant(harvest(ple4), method = "ac")
setGeneric("coef", function(object, ...) standardGeneric("coef"))

#' @rdname coef-methods
#' @aliases coef,FLa4aFit-method
setMethod("coef", signature(object = "a4aFitSA"),
  function(object) {
	  coef(pars(object))
  })


#' @rdname coef-methods
#' @aliases coef,FLa4aFit-method
setMethod("coef", signature(object = "SCAPars"),
  function(object) {
    list(
      stkmodel = coef(stkmodel(object)),
      qmodel   = coef(qmodel(object)),
      vmodel   = coef(vmodel(object))
    )
  })

#' @rdname coef-methods
#' @aliases coef,FLa4aFit-method
setMethod("coef", signature(object = "a4aStkParams"),
  function(object) {
      object @ params
  })

#' @rdname coef-methods
#' @aliases coef,FLa4aFit-method
setMethod("coef", signature(object = "submodels"),
  function(object) {
      lapply(object, coef)
  })


#' @rdname coef-methods
#' @aliases coef,FLa4aFit-method
setMethod("coef", signature(object = "submodel"),
  function(object) {
      object @ params
  })


# -------------------------------------------------------------------
#
#
#    coef<-  methods
#
#
# -------------------------------------------------------------------

setGeneric("coef<-", function(object, ..., value) standardGeneric("coef<-"))

#' @rdname coef-methods
#' @aliases coef,FLa4aFit-method
setMethod("coef<-", signature(object = "a4aFitSA", value = "numeric"),
  function(object, ..., value) {
    coef(object @ pars) <- value
    object
  })


#' @rdname coef-methods
#' @aliases coef,FLa4aFit-method
setMethod("coef<-", signature(object = "SCAPars", value = "numeric"),
  function(object, ..., value) {
    v <- coef(object)
    old <- unlist(v)
    new <- rep_len(unlist(value), length = length(old))
    
    coef(object @ stkmodel) <- new[grep("stkmodel", names(old))]
    coef(object @ qmodel) <- new[grep("qmodel.", names(old))]
    coef(object @ vmodel) <- new[grep("vmodel.", names(old))]

    object
  })


#' @rdname coef-methods
#' @aliases coef,FLa4aFitSA-method
setMethod("coef<-", signature(object = "a4aStkParams", value = "numeric"),
  function(object, ..., value) {    
    object @ params[] <- value
    object
  })

#' @rdname coef-methods
#' @aliases coef,FLa4aFitSA-method
setMethod("coef<-", signature(object = "submodels", value = "numeric"),
  function(object, ..., value) {
    v <- coef(object)
    old <- unlist(v)
    new <- rep_len(unlist(value), length = length(old))
    
    for (i in seq_along(object)) {
      object[[i]] @ params[] <- new[grep(object[[i]] @ name, names(old))]  
    }
    object
  })


#' @rdname coef-methods
#' @aliases coef,FLa4aFit-method
setMethod("coef<-", signature(object = "submodel", value = "numeric"),
  function(object, ..., value) {
      object @ params[] <- value
      object
  })

# -------------------------------------------------------------------
#
#
#    vcov  methods
#
#
# -------------------------------------------------------------------


#' @rdname vcov-methods
#' @aliases vcov,FLa4aFitSA-method
setMethod("vcov", signature(object = "a4aFitSA"),
  function(object) {
    vcov(pars(object))
  })


#' @rdname vcov-methods
#' @aliases vcov,FLa4aFitSA-method
setMethod("vcov", signature(object = "SCAPars"),
  function(object) {
    list(
      stkmodel = vcov(stkmodel(object)),
      qmodel   = vcov(qmodel(object)),
      vmodel   = vcov(vmodel(object))
    )
  })

#' @rdname vcov-methods
#' @aliases vcov,FLa4aFitSA-method
setMethod("vcov", signature(object = "a4aStkParams"),
  function(object) {
      object @ vcov
  })

#' @rdname vcov-methods
#' @aliases vcov,FLa4aFitSA-method
setMethod("vcov", signature(object = "submodels"),
  function(object) {
      lapply(object, vcov)
  })


#' @rdname coef-methods
#' @aliases coef,FLa4aFit-method
setMethod("vcov", signature(object = "submodel"),
  function(object) {
      object @ vcov
  })


# -------------------------------------------------------------------
#
#
#    vcov<-  methods
#
#
# -------------------------------------------------------------------


#' @rdname vcov-methods
#' @aliases vcov,FLa4aFit-method
setMethod("vcov<-", signature(object = "a4aFitSA", value = "numeric"),
  function(object, ..., value) {
    vcov(object @ pars) <- value
    object
  })


#' @rdname vcov-methods
#' @aliases vcov,FLa4aFit-method
setMethod("vcov<-", signature(object = "SCAPars", value = "numeric"),
  function(object, ..., value) {
    v <- vcov(object)
    old <- unlist(v)
    new <- rep_len(unlist(value), length = length(old))
    
    vcov(object @ stkmodel) <- new[grep("stkmodel", names(old))]
    vcov(object @ qmodel) <- new[grep("qmodel.", names(old))]
    vcov(object @ vmodel) <- new[grep("vmodel.", names(old))]

    object
  })


#' @rdname vcov-methods
#' @aliases vcov,FLa4aFitSA-method
setMethod("vcov<-", signature(object = "a4aStkParams", value = "numeric"),
  function(object, ..., value) {    
    object @ vcov[] <- value
    object
  })

#' @rdname vcov-methods
#' @aliases vcov,FLa4aFitSA-method
setMethod("vcov<-", signature(object = "submodels", value = "numeric"),
  function(object, ..., value) {
    v <- vcov(object)
    old <- unlist(v)
    new <- rep_len(unlist(value), length = length(old))
    
    for (i in seq_along(object)) {
      object[[i]] @ vcov[] <- new[grep(object[[i]] @ name, names(old))]  
    }
    object
  })


#' @rdname coef-methods
#' @aliases coef,FLa4aFit-method
setMethod("vcov<-", signature(object = "submodel", value = "numeric"),
  function(object, ..., value) {
      object @ vcov[] <- value
      object
  })


# -------------------------------------------------------------------
#
#
#    predict  methods
#
#
# -------------------------------------------------------------------


#' @rdname predict-methods
#' @aliases predict,FLa4aFit-method
setMethod("predict", signature(object = "a4aFitSA"),
  function(object) {
    predict(pars(object))
  })


#' @rdname coef-methods
#' @aliases coef,FLa4aFit-method
setMethod("predict", signature(object = "SCAPars"),
  function(object) {
    list(
      stkmodel = predict(stkmodel(object)),
      qmodel   = predict(qmodel(object)),
      vmodel   = predict(vmodel(object))
    )
  })

#' @rdname coef-methods
#' @aliases coef,FLa4aFit-method
setMethod("predict", signature(object = "a4aStkParams"),
  function(object) {
      ages <- range(object)["min"]:range(object)["max"]
      years <- range(object)["minyear"]:range(object)["maxyear"]
      cnames <- rownames(coef(object))
      df <- expand.grid(age = ages,
                        year = years)
      X <- getX(object @ fMod, df)
      b <- coef(object)[grep("fMod", cnames)]
      fit <- exp(c(X %*% b))
      harvest <- FLQuant(array(fit, dim = c(length(ages), length(years)), 
                               dimnames = list(age = ages, year = years)), 
                         units = "f")

      X <- getX(object @ n1Mod, data.frame(age = ages[-1]))
      b <- coef(object)[grep("n1Mod", cnames)]
      fit <- c(NA, exp(c(X %*% b) + object @ centering))
      ny1 <- FLQuant(array(fit, dim = length(ages), 
                               dimnames = list(age = ages)))      

      X <- getX(object @ srMod, data.frame(year = years))
      b <- coef(object)[grep("rMod", cnames)]
      fit <- c(exp(c(X %*% b) + object @ centering))
      rec <- FLQuant(array(fit, dim = c(1,length(years)), 
                               dimnames = list(age = ages[1], year = years)))      


      list(harvest = harvest, rec = rec, ny1 = ny1)
})


#' @rdname coef-methods
#' @aliases coef,FLa4aFit-method
setMethod("predict", signature(object = "submodels"),
  function(object, ...) {
      lapply(object, predict)
  })


#' @rdname coef-methods
#' @aliases coef,FLa4aFit-method
setMethod("predict", signature(object = "submodel"),
  function(object, ...) {
      ages <- range(object)["min"]:range(object)["max"]
      years <- range(object)["minyear"]:range(object)["maxyear"]
      df <- expand.grid(age = ages,
                        year = years)
      X <- getX(object @ Mod, df)
      fit <- exp(c(X %*% coef(object)) + object @ centering)
      FLQuant(array(fit, dim = c(length(ages), length(years)), dimnames = list(age = ages, year = years)))
  })


# -------------------------------------------------------------------
#
#
#    residual  methods
#
#
# -------------------------------------------------------------------



# -------------------------------------------------------------------
#
#
#    simulate  methods
#
#
# -------------------------------------------------------------------

#' @rdname au-methods
#' @aliases au,a4aFitSA,missing,a4aFitSA-method
setMethod("genFLStock", c("a4aFitSA", "missing", "missing", "missing"), 
  function(object, ...){
    simulate(object)
})


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
#' @rdname coef-methods
#'
#' @examples
#' data(ple4)
#' genFLQuant(harvest(ple4), method = "ac")
setGeneric("simulate", useAsDefault = stats::simulate)

#' @rdname coef-methods
#' @aliases coef,FLa4aFit-method
setMethod("simulate", signature(object = "a4aFitSA"),
  function(object, nsim = 1, seed = NULL) {

    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    if (is.null(seed)) {
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    } else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }

    object <- pars(object)

    years <- range(object @ stkmodel)[c("minyear","maxyear")]
    ages <- range(object @ stkmodel)[c("min","max")]
  
    #
    # Build design matrix for catches only
    #
    full.df <- expand.grid(age  = ages[1]:ages[2],
                           year = years[1]:years[2])[2:1]
 
#    if (!is.null(covar)) {
#    # add in covariates to data.frame - it is easiest to provide covariates in one list
#    tmp <- 
#      lapply(seq_along(covar), 
#        function(i) {
#          x <- as.data.frame(covar[[i]])[c(1,2,7)]
#          if (length(unique(x $ age)) == 1) x <- x[names(x) != "age"]
#          if (length(unique(x $ year)) == 1) x <- x[names(x) != "year"]
#          names(x) <- gsub("data", names(covar)[i], names(x))
#          x
#        })
#    covar.df <- tmp[[1]]
#    for (i in seq(length(covar) - 1)) covar.df <- merge(covar.df, tmp[[i + 1]], all = TRUE, sort = FALSE)
#
#    full.df <- merge(full.df, covar.df, all.x = TRUE, all.y = FALSE)
#    } 

    # make sure contrasts are set to sumto zero to match fit
    opts <- options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly")) 
  
    # f model matrix
    Xf <- Matrix(getX(object @ stkmodel @ fMod, full.df))

    # initial age structure model matrix
    Xny1 <- getX(object @ stkmodel @ n1Mod, subset(full.df, year == min(year) & age > min(age)))

    # Q model matrix  
    fleet.names <- c("catch", names(object @ qmodel))
    Xqlist <- lapply(seq_along(object @ qmodel), function(i) getX(object @ qmodel[[i]] @ Mod, subset(full.df, fleet == fleet.names[i+1])))
    Xq <- as.matrix(do.call(bdiag, Xqlist))  
  
    # var model matrix
    Xvlist <- lapply(1:length(fleet.names), function(i) getX(object @ vmodel[[i]] @ Mod, subset(full.df, fleet == fleet.names[i])))
    Xv <- as.matrix(do.call(bdiag, Xvlist))   
    
    # now separate the sr model element
    facs <- strsplit(as.character(object @ stkmodel @ srMod)[length(object @ stkmodel @ srMod)], "[+]")[[1]]
    facs <- gsub("(^ )|( $)", "", facs) # remove leading and trailing spaces
    a4as <- grepl(paste("(^",c("bevholt", "ricker","hockey","geomean"),"[(])", collapse = "|", sep = ""), facs)

    # internal r model matrix
    if (sum(a4as) == 0) rmodel <- object @ stkmodel @ srMod else rmodel <- ~ factor(year) 
    Xr <- getX(rmodel, subset(full.df, age == min(age)))

    # reset options
    options(opts)

    # always simulate from b distribution for SA class.  If you want fitted values do FLStock + a4aFit(a4aFitSA)
    b.sim <- Matrix(simulate(object, nsim = nsim) @ stkmodel @ params @ .Data)

    # matrix of predictions
    Xbeta <- bdiag(Xf, Xny1, Xr) %*% b.sim

    # plusgroup?
    rng <- range(object @ stkmodel)
    plusgrp <- !is.na(rng["plusgroup"]) && rng["plusgroup"] >= rng["max"]

    # unpack m - good for recycling
    Ms   <- c(m(object) @ .Data)
 
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
    Ns <- Ns * exp(object @ stkmodel @ centering)
 
    zfrac <- Fs / Zs * (1 - exp(-Zs))

    dmns <- list(age    = paste(ages[1]:ages[2]), 
                 year   = paste(years[1]:years[2]),
                 unit   = "unique", 
                 season = "all", 
                 area   = "unique", 
                 iter   = paste(seq(dim(b.sim)[2])))
               
    dms <- unname(sapply(dmns, length))

    out <- FLStock(
             stock.n = FLQuant(Ns, dim = dms, dimnames = dmns, units = stkmodel(object) @ units),
             catch.n = FLQuant(zfrac * Ns, dim = dms, dimnames = dmns, units = stkmodel(object) @ units),
             harvest = FLQuant(Fs, dim = dms, dimnames = dmns, units = "f"),
             m       = m(object),
             range   = object @ stkmodel @ range)

    out
  }
)


#' @rdname coef-methods
#' @aliases coef,FLa4aFit-method
setMethod("simulate", signature(object = "SCAPars"),
  function(object, nsim = 1, seed = NULL, iter = NULL) {

    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    if (is.null(seed)) {
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    } else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }

    # sanity checks
    if (is.null(iter)) {
      if (nsim == 1) iter <- seq(dim(object @ stkmodel @ params)[2])
      if (nsim > 1) iter <- 1
    } else
    {
      if (nsim == 1) {
        if (any(iter > dim(object @ stkmodel @ params)[2]) | any(iter < 0)) {
          message("supplied values of iter are not sensible... simulating from all iters")
          iter <- seq_along(dim(object @ stkmodel @ params)[2])
        }
    } else
      {
        if (length(iter) > 1) stop("if nsim > 1 iter must be of length 1")
        if (iter > dim(object @ stkmodel @ params)[2] | iter < 0) {
          message("supplied values of iter are not sensible... simulating from iter = 1")
          iter <- 1
        }
      }
    }
    
    # get parameter estimates
    stkpars <- object @ stkmodel @ params
    qpars   <- lapply(object @ qmodel, function(x) x @ params)
    vpars   <- lapply(object @ vmodel, function(x) x @ params)
    
    # get parameter variance matrices
    stkvars <- object @ stkmodel @ vcov
    qvars   <- lapply(object @ qmodel, function(x) x @ vcov)
    vvars   <- lapply(object @ vmodel, function(x) x @ vcov)

    # simulate some new params from the first iteration only!
    if (dim(object @ stkmodel @ vcov)[3] == 1) {
      itervar <- rep(1, length(iter))
    } else {
      itervar = iter
    }
    stkparsim <-
      sapply(seq_along(iter),
        function(i) 
          t(mvrnorm(nsim, c(stkpars[,iter[i]]), stkvars[,,itervar[i]])))

    qparsim <- 
      lapply(seq_along(qpars), 
        function(j) 
          sapply(seq_along(iter), 
            function(i) 
              t(mvrnorm(nsim, c(qpars[[j]][,iter[i]]), qvars[[j]][,,itervar[i]]))))

    vparsim <- 
      lapply(seq_along(vpars), 
        function(j) 
          sapply(seq_along(iter), 
            function(i) 
              t(mvrnorm(nsim, c(vpars[[j]][,iter[i]]), vvars[[j]][,,itervar[i]]))))
    
    # load simpars into a SCAPars object and return
    out <- object
    
    if (nsim == 1) { 
      out @ stkmodel @ params <- object @ stkmodel @ params
    } else {
      out @ stkmodel @ params <- propagate(object @ stkmodel @ params[,iter], nsim)
    }
    out @ stkmodel @ params[] <- c(stkparsim)

    for (i in seq_along(out @ qmodel)) {
      if (nsim == 1) { 
        out @ qmodel[[i]] @ params <- object @ qmodel[[i]] @ params
      } else {
        out @ qmodel[[i]] @ params <- propagate(object @ qmodel[[i]] @ params[,1], nsim)
      }
      out @ qmodel[[i]] @ params[] <- c(qparsim[[i]])
    }

    for (i in seq_along(out @ vmodel)) {
      if (nsim == 1) { 
        out @ vmodel[[i]] @ params <- object @ vmodel[[i]] @ params
      } else {
        out @ vmodel[[i]] @ params <- propagate(object @ vmodel[[i]] @ params[,1], nsim)
      }
      out @ vmodel[[i]] @ params[] <- c(vparsim[[i]])
    }

    ####
    # note we set the variance matrices to zero
    # since having a variance no longer makes sense...
    ####
    vcov(out) <- 0

    return(out)
  })







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
  e1 + pars(e2)
})



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

setMethod("+", c("FLIndices", "a4aFitSA"), function(e1, e2) 
{
  e1 + pars(e2)
})


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

