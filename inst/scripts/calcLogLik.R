#' Extracts the ADMB output from the standard ADMB output files given a directory
#'
#'
#' @param stock an FLStock object containing catch and stock information
#' @param indices an FLStock object containing catch and stock information
#' @return a list containing the summary of the ADMB optimisation
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
calcLogLik <- function(coefs, object, stock, indices, srrCV = -1,
                   fmodel.extra = NULL, qmodel.extra = NULL, rmodel.extra = NULL,
                   wkdir = NULL, verbose = FALSE, eval.at = NULL,  vform = list( ~ 1, ~ 1, ~ 1))
{

  # set up observation data frame
  basedata <- c(list(catch=catch.n(stock)[drop=TRUE]),
                     lapply(indices, function(x) index(x)[drop=TRUE] ))
  center.log <- sapply(basedata, function(x) mean(log(x), na.rm = TRUE))

  doone <- function(f)
  {
    x <- basedata[[f]]
    obs <- as.vector(x) * exp(- center.log[f]) # TODO convert a4a executable to take log obs as input
    y <- as.numeric(colnames(x)[col(x)])
    a <- as.numeric(rownames(x)[row(x)])
    ret <- data.frame(fleet = as.integer(f), year = as.integer(y), age = as.integer(a), obs = obs)
    ret <- ret[complete.cases(ret),]
  }

  all <- do.call(rbind,lapply(1:length(basedata), doone))

  all.df <- all
  all.df $ x <- 1

  # set up auxilliary data
  stockInfo <- sapply(c("m","mat","stock.wt","catch.wt"), function(x) get(x)(stock)[drop=TRUE], simplify = FALSE)
  stockInfo $ fbar <-  unname(range(stock)[c("minfbar","maxfbar")])
  stockInfo $ plusgroup  <- as.integer( !is.na(range(stock)["plusgroup"]), range(stock)["plusgroup"] >= range(stock)["max"] )

  # other arguments
  stockInfo $ surveytime <- unname(sapply(indices, function(x) mean(c(dims(x) $ startf, dims(x) $ endf))))
  
  f.df <- expand.grid(age = as.numeric(rownames(basedata[[1]])),
                         year = as.numeric(colnames(basedata[[1]])), 
                          x = 1)[c(2,1,3)]

  opts <- options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly")) 

  Xf <- getX(object @ models $ fmodel, f.df)

  
  Xq1 <- as.matrix(do.call(bdiag, lapply(1:length(indices), 
                                        function(i) getX( object @ models $ qmodel[[i]], 
                                                      subset(all.df, fleet == (i+1)) ) 
                                        )))

  Xq <- matrix(0, nrow(f.df), ncol(Xq1))
  Xq[apply(f.df[1:2], 1, paste, collapse = ".") %in% apply(subset(all.df, fleet > 1)[2:3], 1, paste, collapse = "."),] <- Xq1
     

  Xr <- getX(object @ models $ rmodel, 
             expand.grid(year = as.numeric(colnames(basedata[[1]])), x = 1))

  cv.df <- expand.grid(age = as.numeric(rownames(basedata[[1]])),
                       year = as.numeric(colnames(basedata[[1]])), 
                         fleet = 1 + 0:length(indices), x = 1)[c(3,2,1,4)]

  Xcv1 <- as.matrix(do.call(bdiag, lapply(1 + 0:length(indices), function(f) getX(vform[[f]], subset(all.df, fleet == f)))))

  Xcv <- matrix(0, nrow(cv.df), ncol(Xcv1))
  Xcv[apply(cv.df[1:3], 1, paste, collapse = ".") %in% apply(all[1:3], 1, paste, collapse = "."),] <- Xcv1
  
  options(opts) # reset options
                         
  X <- list(Xf = Xf, Xq = Xq, Xr = Xr, Xcv = Xcv)

  if (is.matrix(coefs)) {
    par <- list(f = coefs[,grep("fpar", colnames(coefs))],
                q = coefs[,grep("qpar", colnames(coefs))],
                r = coefs[,grep("rpar", colnames(coefs))],
                ry1 = coefs[,grep("ry1", colnames(coefs))],
                cv = coefs[,grep("logSd", colnames(coefs))])
  } else {
    par <- list(f = coefs[grep("fpar", names(coefs))],
                q = coefs[grep("qpar", names(coefs))],
                r = coefs[grep("rpar", names(coefs))],
                ry1 = coefs[grep("ry1", names(coefs))],
                cv = coefs[grep("logSd", names(coefs))])
    par <- lapply(par, matrix, nrow = 1)
  }


  out <- .Call("calcLogLik", all, stockInfo, X, par, PACKAGE = "FLa4a")

  out $ ll
}

