####################################################################
# Experimental code
# A confusing set of experimental code 
####################################################################

#==================================================================== 
# plot
#==================================================================== 
#setMethod("plot", signature(x = "a4aFit", y = "FLStock"),
#  function (x, y, ratio = 1.5, file = "", onefile = TRUE, 
#            #what = c("N","F","Q","Res"),
#            what = c("N", "F", "Res"), 
#            Ftext = TRUE, ask = TRUE, ...) 
#  {

#  op <- par(ask = ask, no.readonly = TRUE)

#  ages  <- as.numeric(dimnames(stock.n(x)) $ age)
#  years <- as.numeric(dimnames(stock.n(x)) $ year)

#  cols <- 
#     c(rgb(215, 48, 39, max = 255), 
#       rgb(252, 141, 89, max = 255),
#       rgb(254, 224, 144, max = 255),
#       rgb(224, 243, 248, max = 255),
#       rgb(145, 191, 219, max = 255),
#       rgb(69, 117, 180, max = 255))

#  zeros <- function(n) paste(rev(rep(c("0", "0", " 0"), length = n)), collapse = "")
#  age.col <- colorRampPalette(cols[c(1, 6)])(length(ages))
#      
#  if ("N" %in% what) {  
#    # N plot
#    scale <- ceiling(max(log(stock.n(x)[drop=TRUE], 10))) - 2
#    
#    matplot(years, t(stock.n(x)[drop=TRUE]) * 10^{-scale},  
#            ylab =paste0("N at age ('",zeros(scale),"s)"), xlab = "Year", 
#            col = age.col,
#            type ='l', lty = 1, lwd = 2, las = 1, main = "N-at-age in the Stock")
#    lx = max(years) + diff(range(years))*0.04
#    ly = (max(stock.n(x)[drop=TRUE]) + diff(range(stock.n(x)[drop=TRUE]))*0.04 ) * 10^{-scale}
#    legend(lx, ly, legend = paste("age", ages), col = age.col, lty = 1, lwd = 2, xjust = 1, yjust = 1)

#  }

#  if ("F" %in% what) {  

#    fest <- harvest(x)[drop=TRUE]

#    # F plot
#    matplot(years, t(fest),  
#            ylab = "F at age", xlab = "Year", 
#            col = colorRampPalette(cols[c(1, 6)])(length(ages)),
#            type = 'l', lty = 1, lwd = 2, las = 1, main = "F-at-age")
#    lx = max(years) + diff(range(years))*0.04
#    ly = max(fest) + diff(range(fest))*0.04
#    legend(lx, ly, legend = paste("age", ages), col = age.col, lty = 1, lwd = 2, xjust = 1, yjust = 1)


#    # another F plot
#    p <- matrix.plot(fest, cols = cols, xlab = "Year", ylab = "Age", main = "F-at-age", ymin = min(ages), xmin = min(years), text = Ftext)
#    print(p)

#    # Yet another F plot
#    p <- wireframe(fest, drape = TRUE, colorkey = FALSE,
#             screen = list(z = 240, x = -60), 
#             col.regions = colorRampPalette(rev(cols))(100),
#             panel.aspect = 1/ratio, aspect=c(length(years)/length(ages), 2),
#             ylab = "Year", xlab = "Age", zlab = "F", main = "F-at-age",
#             #par.settings = list(axis.line = list(col = "transparent")),
#             par.box = c(col = "transparent"))
#    print(p)

#  }

##    sub.dev.new()
#    # fbar plot
##    plotError(years, fbar(y)[,,,,,"mean",drop=TRUE], sqrt(fbar(x)[,,,,,"var",drop=TRUE]), 
##              ylab = 'Fbar', xlab = "Year", cols = colorRampPalette(cols[5:6])(3))

##    sub.dev.new()
#    # ssb plot
##    plotError(years, ssb(y)[,,,,,"mean",drop=TRUE] * 1e-3, sqrt(ssb(x)[,,,,,"var",drop=TRUE]) * 1e-3, 
##              ylab = "SSB ('000 tonnes)", xlab = "Year", cols = colorRampPalette(cols[5:6])(3))

#  if ("Q" %in% what) {  

#    # q plots
#    for (i in 1:length(x @ logq)) {
#      qest <- logq(x)[[i]]
#      qages  <- as.numeric(dimnames(qest) $ age)
#      qyears <- as.numeric(dimnames(qest) $ year)

#      matplot(qages, qest[drop=TRUE],  
#              ylab="log catchability", xlab="Age", main = names(x @ logq)[i],
#              col = colorRampPalette(cols[c(1, 6)])(length(qyears)),
#              type ='l', lty = 1, lwd = 2, las = 1)

#      p <- wireframe(qest[drop=TRUE], drape = TRUE, colorkey = FALSE,
#             # screen value tuned to ple4.indices[1:2]
#             screen = list(z = ifelse(i==1, -1, 1) * 30, x = -60), 
#             col.regions = colorRampPalette(rev(cols))(100),
#             panel.aspect = 1/ratio, ylab = "Year", xlab = "Age", 
#             zlab = "log catchability", main = paste("log catchability-at-age of", names(x @ logq)[i]),
#             par.box = c(col = "transparent"))
#      print(p) 
#    }

#    # q plots
#    rec.age <- range(y)["min"]
#    ssb.x <- ssb(y)[,seq(rec.age) - length(ssb(y)) - 1][drop=TRUE]
#    rec.y <- rec(y)[,-seq(rec.age)][drop=TRUE]
#    plot(ssb.x, rec.y, pch = 19, ann = FALSE, xlim = c(0, max(ssb.x)), ylim = c(0, max(rec.y)))
#    title(main = "Stock and Recruitment", xlab = "SSB", ylab = paste("Recruitment (age ", rec.age, ")", sep = ""))
#    # TODO  add in stock recruit fit if a model was used.
#  }

#  if ("Res" %in% what) {  
#    
#    par(oma = c(4,0,0,0))
#    
#    # Residual plot 
#    div <- list(c(1,1), c(2,1), c(3,1), c(2,2), c(3,2), c(3,2), c(3,3), c(3,3), c(3,3), c(4,3), c(4,3), c(4,3))
#    nind <- length(x @ index)
#    op2 <- par(mfrow = div[[nind + 1]], mar = c(4,4,1,2), mgp = c(2,1,0), no.readonly = TRUE)
#    ylim <- range(as.numeric(dimnames(stock.n(x)) $ age)) + c(-.5, .5)
#    xlim <- range(as.numeric(dimnames(stock.n(x)) $ year)) + c(-.5, .5)

#    res <- log(catch.n(y) / catch.n(x))
#    ages  <- as.numeric(dimnames(res) $ age)
#    years <- as.numeric(dimnames(res) $ year)

#    bp(rep(years, each = length(ages)), rep(ages, length(years)), 
#       c(res[drop=TRUE]), 
#       ylim = ylim, xlim = xlim, 
#       xlab = 'year', ylab = 'Age', main = "Catch", 
#       scale = 3, las = 1)

#    

##    for(i in 1:nind) {
##      res <- index.lres(x)[[i]]
##      ages  <- as.numeric(dimnames(res) $ age)
##      years <- as.numeric(dimnames(res) $ year)

##      bp(rep(years, each = length(ages)), rep(ages, length(years)), 
##         c(res[drop=TRUE]), 
##         ylim = ylim, xlim = xlim, 
##         xlab = 'year', ylab = 'Age', main = names(x @ logq)[i], 
##         scale = 3, las = 1)
##    }

#    mtext("Standardised residuals: Positive (red) and negative (blue).  \n50% of dots should be light coloured, 95% should be light and medium, \n5% of dots should be dark.", 
#          side = 1, line = 6, adj = 0)
#    par(op2)
#  }
#  
#  par(op)
# 
#})




#==================================================================== 
#    simulate  methods
#==================================================================== 

# setMethod("simulate", signature(object = "a4aFitSA"),
#   function(object, nsim = 1, seed = NULL) {


#     object <- pars(object)

#     years <- range(object @ stkmodel)[c("minyear","maxyear")]
#     ages <- range(object @ stkmodel)[c("min","max")]
  
#     #
#     # Build design matrix for catches only
#     #
#     full.df <- expand.grid(age  = ages[1]:ages[2],
#                            year = years[1]:years[2])[2:1]
 
# #    if (!is.null(covar)) {
# #    # add in covariates to data.frame - it is easiest to provide covariates in one list
# #    tmp <- 
# #      lapply(seq_along(covar), 
# #        function(i) {
# #          x <- as.data.frame(covar[[i]])[c(1,2,7)]
# #          if (length(unique(x $ age)) == 1) x <- x[names(x) != "age"]
# #          if (length(unique(x $ year)) == 1) x <- x[names(x) != "year"]
# #          names(x) <- gsub("data", names(covar)[i], names(x))
# #          x
# #        })
# #    covar.df <- tmp[[1]]
# #    for (i in seq(length(covar) - 1)) covar.df <- merge(covar.df, tmp[[i + 1]], all = TRUE, sort = FALSE)
# #
# #    full.df <- merge(full.df, covar.df, all.x = TRUE, all.y = FALSE)
# #    } 

#     # make sure contrasts are set to sumto zero to match fit
#     opts <- options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly")) 
  
#     # f model matrix
#     Xf <- Matrix(getX(object @ stkmodel @ fMod, full.df))

#     # initial age structure model matrix
#     Xny1 <- getX(object @ stkmodel @ n1Mod, subset(full.df, year == min(year) & age > min(age)))

#     # Q model matrix  
#     fleet.names <- c("catch", names(object @ qmodel))
#     Xqlist <- lapply(seq_along(object @ qmodel), function(i) getX(object @ qmodel[[i]] @ Mod, subset(full.df, fleet == fleet.names[i+1])))
#     Xq <- as.matrix(do.call(bdiag, Xqlist))  
  
#     # var model matrix
#     Xvlist <- lapply(1:length(fleet.names), function(i) getX(object @ vmodel[[i]] @ Mod, subset(full.df, fleet == fleet.names[i])))
#     Xv <- as.matrix(do.call(bdiag, Xvlist))   
    
#     # now separate the sr model element
#     facs <- strsplit(as.character(object @ stkmodel @ srMod)[length(object @ stkmodel @ srMod)], "[+]")[[1]]
#     facs <- gsub("(^ )|( $)", "", facs) # remove leading and trailing spaces
#     a4as <- grepl(paste("(^",c("bevholt", "ricker","hockey","geomean"),"[(])", collapse = "|", sep = ""), facs)

#     # internal r model matrix
#     if (sum(a4as) == 0) rmodel <- object @ stkmodel @ srMod else rmodel <- ~ factor(year) 
#     Xr <- getX(rmodel, subset(full.df, age == min(age)))

#     # reset options
#     options(opts)

#     # always simulate from b distribution for SA class.  If you want fitted values do FLStock + a4aFit(a4aFitSA)
#     b.sim <- Matrix(simulate(object, nsim = nsim) @ stkmodel @ params @ .Data)

#     # matrix of predictions
#     Xbeta <- bdiag(Xf, Xny1, Xr) %*% b.sim

#     # plusgroup?
#     rng <- range(object @ stkmodel)
#     plusgrp <- !is.na(rng["plusgroup"]) && rng["plusgroup"] >= rng["max"]

#     # unpack m - good for recycling
#     Ms   <- c(m(object) @ .Data)
 
#     # build stock
#     Fs <- Ns <- array(exp(Xbeta[1:nrow(Xf),]), dim = c(diff(ages)+1, diff(years)+1, ncol(Xbeta)))
#     Ns[] <- NA
#     Ns[-1,1,] <- array(exp(Xbeta[nrow(Xf) + 1:nrow(Xny1),]), dim = c(diff(ages), 1, ncol(Xbeta)))
#     Ns[1,,] <- array(exp(Xbeta[nrow(Xf) + nrow(Xny1) + 1:nrow(Xr),]), dim = c(1, diff(years)+1, ncol(Xbeta)))
#     Zs <- Fs + Ms
#     for (a in 2:dim(Ns)[1]) {
#       Ns[a,-1,] <- Ns[a-1, 1:diff(years),] * exp( - Zs[a-1, 1:diff(years),] )
#     }
#     # if plus group
#     if (plusgrp) {
#       for (y in 1:diff(years)) Ns[a,y+1,] <- Ns[a,y+1,] + Ns[a, y,] * exp( - Zs[a, y,] )
#     } 
#     # apply centering
#     Ns <- Ns * exp(object @ stkmodel @ centering)
 
#     zfrac <- Fs / Zs * (1 - exp(-Zs))

#     dmns <- list(age    = paste(ages[1]:ages[2]), 
#                  year   = paste(years[1]:years[2]),
#                  unit   = "unique", 
#                  season = "all", 
#                  area   = "unique", 
#                  iter   = paste(seq(dim(b.sim)[2])))
               
#     dms <- unname(sapply(dmns, length))

#     out <- FLStock(
#              stock.n = FLQuant(Ns, dim = dms, dimnames = dmns, units = stkmodel(object) @ units),
#              catch.n = FLQuant(zfrac * Ns, dim = dms, dimnames = dmns, units = stkmodel(object) @ units),
#              harvest = FLQuant(Fs, dim = dms, dimnames = dmns, units = "f"),
#              m       = m(object),
#              range   = object @ stkmodel @ range)

#     out
#   }
# )


