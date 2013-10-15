






#' Hello 
#' @name logLik
#' @docType methods
#' @rdname logLik-methods
#' @aliases logLik,FLa4aFit-method
setMethod("logLik", signature(object = "a4aFit"),
  function(object, ...) 
  {  
    val <- -1 * unname(object @ fitSumm["nlogl",])
    attr(val, "nobs") <- unname(object @ fitSumm["nobs",])
    attr(val, "df") <- unname(object@fitSumm["nopar",])
    class(val) <- "logLik"
    val
 })

#' Hello 
#' @name show
#' @docType methods
#' @rdname show-methods
#' @aliases show,FLa4aFit-method
setMethod("plot", signature(x = "a4aFit", y = "FLStock"),
  function (x, y, ratio = 1.5, file = "", onefile = TRUE, what = c("N","F","Q","Res"), Ftext = FALSE, ...) 
  {

  # some checks
  
 
  # add fit to stock in case it hasnt been done
  y <- merge(y, x)


    if (file == "") {
      sub.dev.new <- function() dev.new(width = 7 * ratio, height = 7)
    }
    else { # will need different dev.new if file == pdf
      sub.dev.new <- function() NULL
      file <- if (onefile) paste0(file, ".pdf") else paste0(file, "%03d.pdf")
      pdf(onefile = onefile, file = file, width = 6 * ratio, height = 6)
    }

    ages  <- as.numeric(dimnames(stock.n(x)) $ age)
    years <- as.numeric(dimnames(stock.n(x)) $ year)

    cols <- 
     c(rgb(215, 48, 39, max = 255), 
       rgb(252, 141, 89, max = 255),
       rgb(254, 224, 144, max = 255),
       rgb(224, 243, 248, max = 255),
       rgb(145, 191, 219, max = 255),
       rgb(69, 117, 180, max = 255))

    # N plot
    zeros <- function(n) paste(rev(rep(c("0", "0", " 0"), length = n)), collapse = "")
    age.col <- colorRampPalette(cols[c(1, 6)])(length(ages))
      
  if ("N" %in% what) {  
    
    scale <- ceiling(max(log(stock.n(x)[drop=TRUE], 10))) - 2
    
    sub.dev.new()
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

    sub.dev.new()
    # F plot
    matplot(years, t(fest),  
            ylab = "F at age", xlab = "Year", 
            col = colorRampPalette(cols[c(1, 6)])(length(ages)),
            type = 'l', lty = 1, lwd = 2, las = 1, main = "F-at-age")
    lx = max(years) + diff(range(years))*0.04
    ly = max(fest) + diff(range(fest))*0.04
    legend(lx, ly, legend = paste("age", ages), col = age.col, lty = 1, lwd = 2, xjust = 1, yjust = 1)


    sub.dev.new()
    # another F plot
    p <- matrix.plot(fest, cols = cols, xlab = "Year", ylab = "Age", main = "F-at-age", ymin = min(ages), xmin = min(years), text = Ftext)
    print(p)

    sub.dev.new()
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

      sub.dev.new()
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
      sub.dev.new()
      print(p) 
    }

    # q plots
    sub.dev.new()
    rec.age <- range(y)["min"]
    ssb.x <- ssb(y)[,seq(rec.age) - length(ssb(y)) - 1][drop=TRUE]
    rec.y <- rec(y)[,-seq(rec.age)][drop=TRUE]
    plot(ssb.x, rec.y, pch = 19, ann = FALSE, xlim = c(0, max(ssb.x)), ylim = c(0, max(rec.y)))
    title(main = "Stock and Recruitment", xlab = "SSB", ylab = paste("Recruitment (age ", rec.age, ")", sep = ""))
    # TODO  add in stock recruit fit if a model was used.
  }

  if ("Res" %in% what) {  

    sub.dev.new()
    
    par(oma = c(4,0,0,0))
    
    # Residual plot 
    div <- list(c(1,1), c(2,1), c(3,1), c(2,2), c(3,2), c(3,2), c(3,3), c(3,3), c(3,3), c(4,3), c(4,3), c(4,3))
    nind <- length(x @ logq)
    op <- par(mfrow = div[[nind + 1]], mar = c(4,4,1,2), mgp = c(2,1,0))
    ylim <- range(as.numeric(dimnames(stock.n(x)) $ age)) + c(-.5, .5)
    xlim <- range(as.numeric(dimnames(stock.n(x)) $ year)) + c(-.5, .5)

    res <- catch.lres(x)
    ages  <- as.numeric(dimnames(res) $ age)
    years <- as.numeric(dimnames(res) $ year)

    bp(rep(years, each = length(ages)), rep(ages, length(years)), 
       c(res[drop=TRUE]), 
       ylim = ylim, xlim = xlim, 
       xlab = 'year', ylab = 'Age', main = "Catch", 
       scale = 3, las = 1)

    

    for(i in 1:nind) {
      res <- index.lres(x)[[i]]
      ages  <- as.numeric(dimnames(res) $ age)
      years <- as.numeric(dimnames(res) $ year)

      bp(rep(years, each = length(ages)), rep(ages, length(years)), 
         c(res[drop=TRUE]), 
         ylim = ylim, xlim = xlim, 
         xlab = 'year', ylab = 'Age', main = names(x @ logq)[i], 
         scale = 3, las = 1)
    }

    mtext("Standardised residuals: Positive (red) and negative (blue).  \n50% of dots should be light coloured, 95% should be light and medium, \n5% of dots should be dark.", 
          side = 1, line = 6, adj = 0)

  }

    if (file != "") dev.off()
  }
)


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
setMethod("coef", signature(object = "a4aFit"),
  function(object) {
	object @ coefficients[1:object@fit.sum["nopar"]]
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
#' @rdname vcov-methods
#'
#' @examples
#' data(ple4)
#' genFLQuant(harvest(ple4), method = "ac")
setGeneric("vcov", function(object, ...) standardGeneric("vcov"))

#' @rdname vcov-methods
#' @aliases vcov,FLa4aFit-method
setMethod("vcov", signature(object = "a4aFit"),
  function(object) {
	  object @ covariance[1:object@fit.sum["nopar"], 1:object@fit.sum["nopar"]]
  })


