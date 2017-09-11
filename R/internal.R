###################################################################################
# internal functions
###################################################################################

# startup message
.onAttach <- function(libname, pkgname)
{
  ## TODO find out sep char for environment vars on macs
  sep <- if (os.type("linux")) ":" else if (os.type("windows")) ";" else ","
  path <- paste0(a4a.dir(), sep, Sys.getenv("PATH"))
  Sys.setenv(PATH=path)

  ## message with version number
  tbl <- library(help = FLa4a)$info[[1]]
  version <- gsub(" |[a-zA-z]|:", "", tbl[grep("Version:",tbl)])
  msg <- paste0("This is FLa4a ", version,". For overview type \'help(package=\"FLa4a\")\'\n")
  packageStartupMessage(msg)

  # check 64 bit platform in windows
  if(os.type("windows") && !grep("x86", sessionInfo()$running))
    stop("a4a executable in this package has been compiled for a 64 bit OS,
      please get the i386 version on the FLa4a release page at
      https://github.com/flr/FLa4a/releases")
  
  #
  check.executable()
}

# returns the location on the file system of the ADMB executable
a4a.dir <- function () 
{
#  if (os.type("mac")) {
#    fnm <- system.file(paste("bin/mac/", os.32or64bit(), "bit", sep = ""), package = "FLa4a")
#  }
#  else if (os.type("linux")) {
  if (os.type("linux")) {
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


# returns TRUE if correct operating system is passed as an argument
#os.type <- function (type = c("linux", "mac", "windows", "else")) 
os.type <- function (type = c("linux", "windows", "else")) 
{
  type = match.arg(type)
  if (type == "windows") {
    return(.Platform$OS.type == "windows")
  }
#  else if (type == "mac") {
#   result = (file.info("/Library")$isdir && file.info("/Applications")$isdir)
#    if (is.na(result)) {
#      result = FALSE
#    }
#    return(result)
#  }
  else if (type == "linux") {
#    return((.Platform$OS.type == "unix") && !os.type("mac"))
    return(.Platform$OS.type == "unix")
  }
  else if (type == "else") {
    return(TRUE)
  }
  else {
    stop("This shouldn't happen.")
  }
}

# finds the size of the operating system addresses
os.32or64bit <- function () 
{
  return(ifelse(.Machine$sizeof.pointer == 4, "32", "64"))
}


# Checks that the executable can be run by the user
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
	expand.grid(age = if (is(thing, 'FLIndexBiomass') ) NA else range(thing)["min"]:range(thing)["max"], 
				year = range(thing)["minyear"]:range(thing)["maxyear"])[2:1]
}

# local utility
write.t <- function(x, file, ...) write.table(x, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t', file = file, append = TRUE)

write.t.sparse <- function(x, file, ...) {
	x <- as(x, "dsCMatrix")
	cat("\n# i\n", x @ i, "\n# p\n", x @ p, "\n# x\n", x @ x, file = file, append = TRUE)  
}  

  
# simulate mvnorm with empirical T (fixes bug in mvrnorm)

mvrEmpT <- function(n, mu, Sigma, tol = 1e-6, empirical=TRUE){
	if(empirical){
		if(n>length(mu)){
			mm <- mvrnorm(n, mu, Sigma, tol=tol, empirical=T)
		} else {
			mm <- mvrnorm(length(mu)+1, mu, Sigma, tol=tol, empirical=T)
			mm <- mm[1:n,]	
		}
	} else {
			mm <- mvrnorm(n, mu, Sigma, tol=tol, empirical=FALSE)
	}
	
	# output with right dims for FLPar
	if(is(mm, "matrix")) t(mm) else (t(t(mm)))

}


# if FLPar param is of dim 1 coerce to matrix and name "intercept"
par2mat <- function(object){
	p0 <- object@params
	dims <- dim(p0)
	if(dims[1]==1){
		m0 <- t(t(p0[drop=TRUE]))
		dimnames(m0)[[2]] <- dimnames(p0)[[1]]
	} else {
		m0 <- t(p0[drop=T])
	} 
	m0
}



