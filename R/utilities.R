#' @title Assorted methods needed by FLa4a
#' @docType methods
#' @name methods
#' @rdname assorted-methods
#' @aliases is.empty
#' @section is.empty:
#' Method \code{is.empty} checks if object is empty. Takes any object returns a logical, \code{TRUE} if object is of length 0.
#' @examples
#' is.empty(list())
#' is.empty(list(a=2))

is.empty <- function(object) {
	length(object) == 0
}

#' @rdname assorted-methods
#' @aliases pars2dim,FLPar-method
#' @section pars2dim:
#' Checks that the name of the second dimension in params is "iter". For internal use, not very interesting for users. It takes a \code{FLPar} object and returns a \code{logical}
#' @examples
#' pars2dim(FLPar())
#' pars2dim(FLPar(array(dim=c(1,1,1))))
setMethod("pars2dim", "FLPar", function(object) {

	dnm <- dimnames(object)
	names(dnm)[2]=="iter"

})

#' @rdname assorted-methods
#' @aliases getYidx getYidx-methods getYidx,FLQuant-method
#' @section getYidx:
#' Gets the FLQuant's numeric id for a vector of "years". For internal use, not very interesting for users. It takes a \code{FLQuant} and a \code{vector} of years and returns a \code{numeric vector} that can be used to subset the \code{FLQuant}.
#' @examples
#' data(ple4)
#' flq <- catch(ple4)
#' getYidx(flq, 2000:2004)
#' flq[, getYidx(flq, 2000:2004)]

setGeneric("getYidx", function(object, ...) standardGeneric("getYidx"))
setMethod("getYidx", "FLQuant", function(object, year) {
	yrs <- dimnames(object)[[2]]
	if(sum(year>0)>0){
		idx <- match(as.character(year), yrs, nomatch=0)
		if(sum(idx)>0){
			idx
		} else {
			year
		}
	} else {
		length(yrs)+year+1
	} 

})

#' @rdname assorted-methods
#' @aliases niters niters-methods niters,FLModelSim-method
#' @section niters:
#' Compute number of iterations. Takes an object of any \code{FLR} class and returns a \code{numeric}.
setGeneric("niters", function(object, ...) standardGeneric("niters"))
setMethod("niters", "FLModelSim", function(object){
	dim(params(object))[2]
})





#' Adds the location of the a4a executable to the search path.  This function is called when FLa4a is attached
#'
#'
#' @param libname required
#' @param pkgname required
#' @return a pointer to the environment in which summaries of the data reside
#' @note \code{.onAttach} is intended to be used internally
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
#' @examples
#' # set up grapical options
#'  lattice.options(default.theme = geta4aLatticeOptions())
geta4aLatticeOptions <- function(...) 
{
  opts <- standard.theme(color = FALSE)
  opts <- modifyList(opts, list(strip.background = list(col = grey(0.95))))
  opts <- modifyList(opts, list(plot.symbol = list(pch = 16, col = "black", cex = 0.7)))
  opts <- modifyList(opts, list(axis.line = list(col = grey(0.8))))
  opts <- modifyList(opts, list(strip.border = list(col = grey(0.8))))

  opts
}

#' Adds the location of the a4a executable to the search path.  This function is called when FLa4a is attached
#'
#'
#' @param libname required
#' @param pkgname required
#' @return a pointer to the environment in which summaries of the data reside
#' @note \code{.onAttach} is intended to be used internally
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
plotIters <- function(object, nsamples = 5, ylab = "", by = "quant", col = 1, ...) {

  by <- match.arg(by, c("year","quant"))
  quant(object) <- "quant"

  x <- setdiff(c("year","quant"), by)

  doOne <- function(p, object) cbind(as.data.frame(apply(object, 1:5, quantile, p)), p = p)

  dat <- do.call(rbind, lapply(c(0.025, 0.5, 0.975), doOne, object = object))

  dat $ by <- dat[[by]]
  dat $ x <- dat[[x]]


  p1 <- xyplot(data ~ x | factor(by), group = p, data = dat,
               type = c("l","g"), scales = list(y = list(relation = "free")), 
               lwd = c(1,2,1), lty = c(2,1,2), col = col,
                ylab = ylab, xlab = x, as.table = TRUE, ...)
  if (nsamples > 0) {
  
    nsamples <- min(nsamples, min(9, dim(object)[6]))
    dat2 <- do.call(rbind, lapply(sample(dims(object) $ iter, nsamples), function(i) cbind(as.data.frame(object[,,,,,i]), p = i)))
    dat2 $ by <- dat2[[by]]
    dat2 $ x <- dat2[[x]]          

    p2 <- xyplot(data ~ x | factor(by), group = p, data = dat2,
                 type = "l", col =  brewer.pal(max(3, nsamples), "Set1"), lty = 1)
               
    p1 + p2
  } else {
    p1
  }
}


#' Reads the hessian file from any ADMB fit.  Used here with the a4a model.
#'
#' @param wkdir the location of the admb output
#' @return a list with the following elements
#' @note \code{getADMBHessian} is intended to be used internally
#' @author not Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
#' @examples
#' # load some data
#' data(ple4)
#' data(ple4.indices)
#' # choose a working directory
#' wkdir <- tempfile()
#' # do an 'assessment' fit wth default settings (not recomended!) and keep results in wkdir
#' fit <- a4a(stock = ple4, indices = ple4.indices, wkdir = wkdir, fit = "assessment")
#' 
#' hessInfo <- getADMBHessian(wkdir)
#' str(hessInfo)
#'
#' # calculate covariance matrix
#' Sigma <- solve(hessInfo $ hes)
#'
#' # plot correlation matrix of parameters
#' Cor <- cov2cor(Sigma)
#' diag(Cor) <- 0
#' colors <- colorRampPalette(c("blue","grey99","red"))(100)
#' image(Matrix(Cor), main = "correlations between parameter estimates", lwd = 0, col.regions = colors)
#' # end of example
getADMBHessian <- function(wkdir) {
## This function reads in all of the information contained in the
## admodel.hes file. Some of this is needed for relaxing the covariance
## matrix, and others just need to be recorded and rewritten to file so ADMB
## "sees" what it’s expecting.
  filename <- file(paste0(wkdir,"/admodel.hes"), "rb")
  on.exit(close(filename))
  num.pars <- readBin(filename, "integer", 1)
  hes.vec <- readBin(filename, "numeric", num.pars^2)
  hes <- matrix(hes.vec, ncol=num.pars, nrow=num.pars)
  hybrid_bounded_flag <- readBin(filename, "integer", 1)
  scale <- readBin(filename, "numeric", num.pars)

  list(num.pars = num.pars, hes = hes, hybrid_bounded_flag = hybrid_bounded_flag, scale = scale)
}

#' Adds the location of the a4a executable to the search path.  This function is called when FLa4a is attached
#'
#'
#' @param libname required
#' @param pkgname required
#' @return a pointer to the environment in which summaries of the data reside
#' @note \code{.onAttach} is intended to be used internally
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
.onAttach <- function(libname, pkgname)
{
  ## TODO find out sep char for environment vars on macs
  sep <- if (os.type("linux")) ":" else if (os.type("windows")) ";" else ","
  path <- paste0(a4a.dir(), sep, Sys.getenv("PATH"))
  Sys.setenv(PATH=path)

  ## message with version number
  tbl <- library(help = FLa4a)$info[[1]]
  version <- gsub(" |[a-zA-z]|:", "", tbl[grep("Version:",tbl)])
  msg <- paste0("This is FLa4a ", version,". For overview type \'help(\"FLa4a-package\")\'\n")
  packageStartupMessage(msg)

  #
  check.executable()
}



#' returns the location on the file system of the ADMB executable
#'
#'
#' @param stock an FLStock object containing catch and stock information
#' @param index an FLIndex object containing survey indices 
#' @return a pointer to the environment in which summaries of the data reside
#' @note \code{a4a.dir} is intended to be used internally and based on the same function in INLA
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
a4a.dir <- function () 
{
  if (os.type("mac")) {
    fnm <- system.file(paste("bin/mac/", os.32or64bit(), "bit", sep = ""), package = "FLa4a")
  }
  else if (os.type("linux")) {
    fnm <- system.file("bin/linux", package = "FLa4a")
  }
  else if (os.type("windows")) {
    fnm <- system.file("bin/windows", package = "FLa4a")
  }
  else {
    stop("Unknown OS")
  }
  if (file.exists(fnm)) {
    return(fnm)
  }
  else {
    stop(paste("FLa4a installation error; no such file", fnm))
  }
}


#' returns TRUE if correct operating system is passed as an argument
#'
#'
#' @param type character string of operating system type
#' @return logical
#' @note \code{os.type} is intended to be used internally and based on the same function in INLA
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
os.type <- function (type = c("linux", "mac", "windows", "else")) 
{
  type = match.arg(type)
  if (type == "windows") {
    return(.Platform$OS.type == "windows")
  }
  else if (type == "mac") {
   result = (file.info("/Library")$isdir && file.info("/Applications")$isdir)
    if (is.na(result)) {
      result = FALSE
    }
    return(result)
  }
  else if (type == "linux") {
    return((.Platform$OS.type == "unix") && !os.type("mac"))
  }
  else if (type == "else") {
    return(TRUE)
  }
  else {
    stop("This shouldn't happen.")
  }
}


#' finds the size of the operating system addresses
#'
#' @return a character giving the size of the operating system address bus
#' @note \code{extractData} is intended to be used internally and based on the same function in INLA
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
os.32or64bit <- function () 
{
  return(ifelse(.Machine$sizeof.pointer == 4, "32", "64"))
}


#' Checks that the executable can be run by the user
#'
#' @return TRUE or FALSE.. primary function is to message the user.
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
#' @export
check.executable <- function() {
 if (os.type("linux")) {
   system(paste0("ls -l ", a4a.dir(), "/a4a > syslog.txt"))
   syslog <- readLines("syslog.txt")
   unlink("syslog.txt")
   
   is.x <- grepl("x", substring(syslog, 1,10))
 
   if (!is.x) {
     message(paste0(
       "Something has gone wrong!\n",
       "the a4a executable has the wrong permissions:\n\t",
          substring(syslog, 1,10), 
     "\nPlease change permissions (in a terminal) to a+x using\n",
       "\tchmod a+x ", a4a.dir(), "/a4a\n",
       "if you installed under sudo you will have to run:\n",
       "\tsudo chmod a+x ", a4a.dir(), "/a4a"))
     }
   
   return(is.x)
 } else { # windows
    return(TRUE)
 }   
}



#' coerces FLPar into list for call
#'
#' @return a list with each parametr as one element
#' @author EJ \email{ernesto.jardim@@jrc.ec.europa.eu}
#setAs("FLPar", "list", function(from){
#	lst <- split(from@.Data, 1:nrow(from))
#	names(lst) <- dimnames(from)[[1]]	
#	lst
#})

#setMethod("iterMedians", "FLPar", function(x){
#	apply(x, -match("iter", names(dimnames(x))), median)
#})

###################################################################################
# internal functions
###################################################################################

# utility to convert to a 2d array

quant2mat <- function(x) {
	out <- x[drop=TRUE]
	dim(out) <- dim(x)[1:2]
	dimnames(out) <- dimnames(x)[1:2]
	if (nrow(out) == 1 && dimnames(out)[[1]] == "all") dimnames(out)[[1]] <- NA_character_  # "all" denotes a biomass survey
	out 
}

# convert to dataframe
list2df <- function(fleet, list.obs, list.var, center.log) {
	x <- list.obs[[fleet]]
	v <- as.vector(list.var[[fleet]])
	year <- as.numeric(colnames(x)[col(x)])
	age <- as.numeric(rownames(x)[row(x)])
	obs <- log(as.vector(x)) - center.log[fleet]
	if (all(is.na(v))) {
		wts <- 1
	} else { # use inverse variance weighting
		wts <-  1 / v # inverse variance weigting
	}
	ret <- data.frame(fleet = fleet, year = year, age = age, obs = obs, weights = wts)
	ret <- ret[!is.na(ret $ obs), ]
	if (any(is.na(ret[,5])) || any(ret[,5] <= 0)) {
		ret[,5] <- 1
		warning("*** NA and/or non-positive variances found in: ", names(list.obs)[fleet], " - all variances set to 1", call. = FALSE)
	}
	ret
}

# build a full data frame first (we will use this for the variance model so it is not a waste)
make.df <- function(fleet, stock, indices) {
	thing <- if (fleet == 1) stock else indices[[fleet - 1]]
	expand.grid(age = if (is.na(range(thing)["min"])) NA else range(thing)["min"]:range(thing)["max"], 
				year = range(thing)["minyear"]:range(thing)["maxyear"])[2:1]
}

# local utility
write.t <- function(x, file, ...) write.table(x, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t', file = file, append = TRUE)

write.t.sparse <- function(x, file, ...) {
	x <- as(x, "dsCMatrix")
	cat("\n# i\n", x @ i, "\n# p\n", x @ p, "\n# x\n", x @ x, file = file, append = TRUE)  
}  

  







