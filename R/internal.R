################################################################################
# internal functions
################################################################################

# startup message
.onAttach <- function(libname, pkgname) {
  ## TODO find out sep char for environment vars on macs
  sep <-
    if (os.type("linux") | os.type("osx")) {
      ":"
    } else if (os.type("windows")){
      ";"
    } else {
      ","
    }

  path <- paste0(a4a.dir(), sep, Sys.getenv("PATH"))
  Sys.setenv(PATH = path)

  ## message with version number
  tbl <- library(help = FLa4a)$info[[1]]
  version <- gsub(" |[a-zA-z]|:", "", tbl[grep("Version:", tbl)])
  msg <- paste0("This is FLa4a ", version,
                ". For overview type \'help(package=\"FLa4a\")\'\n")
  packageStartupMessage(msg)

  # check 64 bit platform in windows
  if (os.type("windows") && grepl("x86", sessionInfo()$running))
    stop("a4a executable in this package has been compiled for a 64 bit OS,\n",
         "       please get the i386 version on the FLa4a release page at\n",
         "       https://github.com/flr/FLa4a/releases")

  check.executable()
}

# returns the location on the file system of the ADMB executable
a4a.dir <- function () {
  fnm <-
    if (os.type("linux")) {
      system.file("bin/linux", package = "FLa4a")
    } else if (os.type("osx")) {
      system.file("bin/osx", package = "FLa4a")
    } else if (os.type("windows")) {
      system.file("bin/windows", package = "FLa4a")
    } else {
      stop("Unknown OS")
    }

  if (!file.exists(fnm))
    stop(paste("FLa4a installation error; no such file", fnm))

  fnm
}


# returns TRUE if correct operating system is passed as an argument
os.type <- function (type = c("linux", "windows", "osx", "else")) {
  type <- match.arg(type)
  if (type == "windows") {
    .Platform$OS.type == "windows"
  }
#  else if (type == "mac") {
#   result = (file.info("/Library")$isdir && file.info("/Applications")$isdir)
#    if (is.na(result)) {
#      result = FALSE
#    }
#    return(result)
#  }
  else if (type == "osx") {
    return(Sys.info()["sysname"] == "Darwin")
  }
  else if (type == "linux") {
#    return((.Platform$OS.type == "unix") && !os.type("mac"))
    return(.Platform$OS.type == "unix" && Sys.info()["sysname"] != "Darwin")
  }
  else if (type == "else") {
    return(TRUE)
  }
  else {
    stop("This shouldn't happen.")
  }
}

# finds the size of the operating system addresses
os.32or64bit <- function () {
  ifelse(.Machine$sizeof.pointer == 4, "32", "64")
}


# Checks that the executable can be run by the user
check.executable <- function() {
 if (os.type("linux")) {
   is.x <- file.access(paste0(a4a.dir(), "/a4a"), 1)[[1]]

   if (is.x != 0) {
     message(paste0(
       "Something has gone wrong!\n",
       "the a4a executable has the wrong permissions:\n\t",
     "\nPlease change permissions (in a terminal) to a+x using\n",
       "\tchmod a+x ", a4a.dir(), "/a4a\n",
       "if you installed under sudo you will have to run:\n",
       "\tsudo chmod a+x ", file.path(a4a.dir(), "a4a")))
     }

   return(is.x)
  }
  # else windows or osx
  TRUE
}

# utility to convert to a 2d array

quant2mat <- function(x) {
  out <- x[drop = TRUE]
  dim(out) <- dim(x)[1:2]
  dimnames(out) <- dimnames(x)[1:2]
  # "all" denotes a biomass survey
  if (nrow(out) == 1 && dimnames(out)[[1]] == "all")
    dimnames(out)[[1]] <- NA_character_
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
  } else {
    # use inverse variance weighting
    wts <-  1 / v # inverse variance weigting
  }
  ret <- data.frame(fleet = fleet, year = year, age = age,
                    obs = obs, weights = wts)
  ret <- ret[!is.na(ret $ obs), ]
  if (any(is.na(ret[, 5])) || any(ret[, 5] <= 0)) {
    ret[, 5] <- 1
    warning("*** NA and/or non-positive variances found in: ",
            names(list.obs)[fleet],
            " - all variances set to 1", call. = FALSE)
  }
  ret
}

# build a full data frame
make.df <- function(fleet, stock, indices) {
  if (fleet <= 2)
    thing <- stock
  else
    thing <- indices[[fleet - 2]]

  if (is(thing, "FLIndexBiomass") || fleet == 2) {
    age <- NA
  } else {
    age <- range(thing)["min"]:range(thing)["max"]
  }
  expand.grid(
    age =  age,
    year = range(thing)["minyear"]:range(thing)["maxyear"]
  )[2:1]
}

# local utility
write.t <- function(x, file, ...) {
  write.table(x, file = file,
              row.names = FALSE, col.names = FALSE, quote = FALSE,
              sep = "\t", append = TRUE)
}

write.t.sparse <- function(x, file, ...) {
  x <- as(x, "dsCMatrix")
  cat("# i\n", x@i,
    "\n# p\n", x@p,
    "\n# x\n", x@x,
    file = file, append = TRUE)
}


# simulate mvnorm with empirical T (fixes bug in mvrnorm)

mvrEmpT <- function(n, mu, Sigma, tol = 1e-6, empirical = TRUE) {
  if (empirical) {
    if (n > length(mu)) {
      mm <- mvrnorm(n, mu, Sigma, tol = tol, empirical = TRUE)
    } else {
      mm <- mvrnorm(length(mu) + 1, mu, Sigma,
                          tol = tol, empirical = TRUE)
      mm <- mm[1:n, ]
    }
  } else {
      mm <- mvrnorm(n, mu, Sigma, tol = tol, empirical = FALSE)
  }

  # output with right dims for FLPar
  if (is(mm, "matrix")) t(mm) else (t(t(mm)))

}


# if FLPar param is of dim 1 coerce to matrix and name "intercept"
par2mat <- function(object){
  p0 <- object@coefficients
  dims <- dim(p0)
  if (dims[1] == 1){
    m0 <- t(t(p0[drop = TRUE]))
    dimnames(m0)[[2]] <- dimnames(p0)[[1]]
  } else {
    m0 <- t(p0[drop = TRUE])
  }
  m0
}

flq_from_range <- function(object) {
  range <- range(object)
  if (all(is.na(range[c("min", "max")])) |
      isTRUE(attr(object, "FLIndexBiomass"))) {
    # fix for biomass indices or any quant that has "all" for the first dim
    FLQuant(
      matrix(NA,
        nrow = 1,
        ncol = range["maxyear"] - range["minyear"] + 1),
        dimnames = list(age = "all",
                        year = range["minyear"]:range["maxyear"]
      )
    )
  } else {
    # the normal case
    FLQuant(
      matrix(NA,
        nrow = range["max"] - range["min"] + 1,
        ncol = range["maxyear"] - range["minyear"] + 1),
        dimnames = list(age = range["min"]:range["max"],
                        year = range["minyear"]:range["maxyear"]
      )
    )
  }
}

formula_has_covariates <- function(x, ok_vars) {
  # asssign default value  
  if (missing(ok_vars)) {
    ok_vars <- c("age", "year", "unit", "season", "area", "iter")
  }
  covars <- setdiff(all.vars(x), ok_vars)
  length(covars) > 0
}


formula_has_covariate_factors <- function(x, ok_vars, warnings = TRUE) {
  # asssign default value  
  if (missing(ok_vars)) {
    ok_vars <- c("age", "year", "unit", "season", "area", "iter")
  }

  funcs <- setdiff(all.names(x), all.vars(x))
  if ("factor" %in% funcs) {
    trms <- terms(x, specials = "factor")
    facs <- attr(trms, "term.labels")[attr(trms, "specials")$factor]
    # only warn if thre are factors for covars
    fac_vars <- sapply(facs, function(form) all.vars(as.formula(paste("~", form))))
    cov_facs <- setdiff(fac_vars, ok_vars)

    if (length(cov_facs)) {
      if (warnings) warning("factor() found! - need to deal with this!")
      return(TRUE)
    }
  }
  return(FALSE)
}



dropMatrixIter <- function(object, iter = 1) {
  dims <- dim(object)
  if (!inherits(object, "array") || length(dims) != 3)
    stop("object must be an array")
  out <- object[,, iter]
  dim(out) <- dim(object)[1:2]
  dimnames(out) <- dimnames(object)[1:2]
  out
}


# ----------------------------------------------
#
#  bevholt, ricker, geomean share definitions with FLSR
#  the others need to be defined by the FLa4a package
#
#  It is the responsibility of FLa4a to maintain the defninition of
#  these SR functoions in the tpl file to match the definition in FLCore
#
# ----------------------------------------------

check_cv <- function(CV) {
  if (!is.numeric(CV))
    stop("CV must be a numeric value not: ", capture.output(print(CV)))
  if (CV <= 0)
    stop ("CV in stock recruit relationship cannot be less than zero")
}

bevholt <- function(CV = 0.5, a = ~ 1, b = ~ 1) {
  check_cv(CV)
  list(srr = "bevholt", a = a, b = b, SPR0 = 1, srrCV = CV, ID = 1)
}

bevholtSV <- function(CV = 0.5, SPR0 = 1, a = ~ 1, b = ~ 1) {
  check_cv(CV)
  list(srr = "bevholtSV", a = a, b = b, SPR0 = SPR0, srrCV = CV, ID = 5)
}

ricker <- function(CV = 0.5, a = ~ 1, b = ~ 1) {
  check_cv(CV)
  list(srr = "ricker", a = a, b = b, SPR0 = 1, srrCV = CV, ID = 2)
}

hockey <- function(CV = 0.5, a = ~ 1, b = ~ 1) {
  check_cv(CV)
  list(srr = "hockey", a = a, b = b, SPR0 = 1, srrCV = CV, ID = 3)
}

geomean <- function(CV = 0.5, a = ~ 1, ...) {
  check_cv(CV)
  list(srr = "geomean", a = a, b = ~ 1, SPR0 = 1, srrCV = CV, ID = 4)
}

none <- function(...) {
  list(srr = "geomean", a = ~ 1, b = ~ 1, SPR0 = 1, srrCV = -1, ID = 4)
}

a4a_srmodel_list <- c("bevholt", "bevholtSV", "ricker", "hockey", "geomean")
flc_srmodel_list <- c("bevholt", "ricker", "geomean")


a4aSRmodelDefinitions <- function(srmodel) {
  srmodelname <- gsub("\\([^()]*\\)", "", srmodel)
  if (srmodelname %in% flc_srmodel_list) {
    # get FLSR definition
    eval(parse(text = paste0("FLCore::", srmodelname, "()")))$model[[3]]
  } else {
    # use a4a definition
    switch(srmodelname,
      # for bevholtSV, spr0 is provided by user,
      bevholtSV = (~ 6 * h * b * ssb /
                     spr0 /
                     ( (b + 5 * ssb) * (a / (1 + a) * 0.8 + 0.2) + b - ssb)
                  )[[2]],
      hockey = (~ a *
                 (ssb + sqrt(b ^ 2 + 0.0025) - sqrt( (ssb - b) ^ 2 + 0.0025))
               )[[2]],
      stop("unknown SR model")
    )
  }
}


isPresenta4aSRmodel <- function(srmodel) {
  facs <- strsplit(as.character(srmodel)[length(srmodel)], "[+]")[[1]]
  facs <- gsub("(^ )|( $)", "", facs) # remove leading and trailing spaces
  grepl(paste("(^", a4a_srmodel_list, "[(])", collapse = "|", sep = ""), facs)
}


geta4aSRmodel <- function(srmodel) {
  facs <- strsplit(as.character(srmodel)[length(srmodel)], "[+]")[[1]]
  # remove leading and trailing spaces
  facs <- gsub("(^ )|( $)", "", facs)
  a4as <-
    grepl(paste("(^", a4a_srmodel_list, "[(])", collapse = "|", sep = ""),
          facs)
  if (sum(a4as) > 1)
    stop("you can only specify one type of stock recruit relationship.")
  if (sum(a4as) == 0)
    "none()"
  else
    facs[a4as]
}
