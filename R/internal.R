###################################################################################
# internal functions
###################################################################################

# startup message
.onAttach <- function(libname, pkgname)
{
  ## message with version number
  tbl <- library(help = FLa4a)$info[[1]]
  version <- gsub(" |[a-zA-z]|:", "", tbl[grep("Version:",tbl)])
  msg <- paste0("This is FLa4a ", version,". For overview type \'help(package=\"FLa4a\")\'\n")
  packageStartupMessage(msg)
}


# utility to convert to a 2d array
quant2mat <- function(x) {
  out <- x[drop = TRUE]
  dim(out) <- dim(x)[1:2]
  dimnames(out) <- dimnames(x)[1:2]
  if (nrow(out) == 1 && dimnames(out)[[1]] == "all") dimnames(out)[[1]] <- NA_character_ # "all" denotes a biomass survey
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
    wts <- 1 / v # inverse variance weigting
  }
  ret <- data.frame(fleet = fleet, year = year, age = age, obs = obs, weights = wts)
  ret <- ret[!is.na(ret$ obs), ]
  if (any(is.na(ret[, 5])) || any(ret[, 5] <= 0)) {
    ret[, 5] <- 1
    warning("*** NA and/or non-positive variances found in: ", names(list.obs)[fleet], " - all variances set to 1", call. = FALSE)
  }
  ret
}

# build a full data frame first (we will use this for the variance model so it is not a waste)
make.df <- function(fleet, stock, indices) {
  thing <- if (fleet == 1) stock else indices[[fleet - 1]]
  expand.grid(
    age = if (is(thing, "FLIndexBiomass")) NA else range(thing)["min"]:range(thing)["max"],
    year = range(thing)["minyear"]:range(thing)["maxyear"]
  )[2:1]
}

# local utility
write.t <- function(x, file, ...) write.table(x, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t', file = file, append = TRUE)

write.t.sparse <- function(x, file, ...) {
  x <- as(x, "symmetricMatrix")
  cat("\n# i\n", x @ i, "\n# p\n", x @ p, "\n# x\n", x @ x, file = file, append = TRUE)
}


# simulate mvnorm with empirical T (fixes bug in mvrnorm)

mvrEmpT <- function(n, mu, Sigma, tol = 1e-6, empirical = TRUE) {
  if (empirical) {
    if (n > length(mu)) {
      mm <- mvrnorm(n, mu, Sigma, tol = tol, empirical = T)
    } else {
      mm <- mvrnorm(length(mu) + 1, mu, Sigma, tol = tol, empirical = T)
      mm <- mm[1:n, ]
    }
  } else {
    mm <- mvrnorm(n, mu, Sigma, tol = tol, empirical = FALSE)
  }

  # output with right dims for FLPar
  if (is(mm, "matrix")) t(mm) else (t(t(mm)))
}


# if FLPar param is of dim 1 coerce to matrix and name "intercept"
par2mat <- function(object) {
  p0 <- object@coefficients
  dims <- dim(p0)
  if (dims[1] == 1) {
    m0 <- t(t(p0[drop = TRUE]))
    dimnames(m0)[[2]] <- dimnames(p0)[[1]]
  } else {
    m0 <- t(p0[drop = T])
  }
  m0
}

flqFromRange <- function(object) {
  range <- range(object)
#  if (all(is.na(range[c("min", "max")]))) {
  if (all(is.na(range[c("min", "max")])) | isTRUE(attr(object, "FLIndexBiomass"))) {
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


dropMatrixIter <- function(object, iter = 1) {
  dims <- dim(object)
  if (!inherits(object, "array") || length(dims) != 3) stop("object must be an array")
  out <- object[,, iter]
  dim(out) <- dim(object)[1:2]
  dimnames(out) <- dimnames(object)[1:2]
  out
}


# these are left over from when you could set linear models for SR params
  #bevholt <- function(a = ~ 1, b = ~ 1, CV = 0.5) {
  #  if (CV <= 0) stop ("CV in stock recruit relationship cannot be less than zero")
  #  list(srr = "bevholt", a = a, b = b, SPR0 = 1, srrCV = CV, ID = 1)
  #}
  #bevholtSV <- function(h = ~ 1, v = ~ 1, SPR0 = 1, CV = 0.5) {
  #  if (CV <= 0) stop ("CV in stock recruit relationship cannot be less than zero")
  #  list(srr = "bevholtSV", a = h, b = v, SPR0 = SPR0, srrCV = CV, ID = 5)
  #}
  #ricker <- function(a = ~ 1, b = ~ 1, CV = 0.5) {
  #  if (CV <= 0) stop ("CV in stock recruit relationship cannot be less than zero")
  #  list(srr = "ricker", a = a, b = b, SPR0 = 1, srrCV = CV, ID = 2)
  #}
  #hockey <- function(a = ~ 1, b = ~ 1, CV = 0.5) {
  #  if (CV <= 0) stop ("CV in stock recruit relationship cannot be less than zero")
  #  list(srr = "hockey", a = a, b = b, SPR0 = 1, srrCV = CV, ID = 3)
  #}
  #geomean <- function(a = ~ 1, CV = 0.5) {
  #  if (CV <= 0) stop ("CV in stock recruit relationship cannot be less than zero")
  #  list(srr = "geomean", a = a, b = ~ 1, SPR0 = 1, srrCV = CV, ID = 4)
  #}
  #none <- function() list(srr = "geomean", a = ~ 1, b = ~ 1, srrCV = -1, ID = 4)

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
  if (CV != -1 && CV <= 0) {
    stop("CV in stock recruit relationship cannot be less than zero")
  }
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
  list(srr = "geomean", a = ~ 1, b = ~ 1, SPR0 = 1, srrCV = 1, ID = 0)
}

a4aSRmodelList <- c("bevholt", "bevholtSV", "ricker", "hockey", "geomean")
flcSRmodelList <- c("bevholt", "ricker", "geomean")


a4aSRmodelDefinitions <- function(srmodel) {
  srmodelName <- gsub("\\([^()]*\\)", "", srmodel)
  if (srmodelName %in% flcSRmodelList) {
    # get FLSR definition
    eval(parse(text=paste0("FLCore::", srmodelName, "()")))$model[[3]]
  } else {
    # use a4a definition
    switch(srmodelName,
      bevholtSV = (~ (6*h*b*ssb) / (spr0* (((a/(1+a)*0.8 + 0.2) + 1)*b + (5*(a/(1+a)*0.8 + 0.2) - 1)*ssb) ))[[2]], # spr0 is provided by user,
      hockey = (~ a * (ssb + sqrt(b^2 + 0.0025) - sqrt((ssb - b)^2 + 0.0025)))[[2]],
      none = NULL,
      stop("unknown SR model")
    )
  }
}


isPresenta4aSRmodel <- function(srMod) {
  #facs <- strsplit(as.character(srMod)[length(srMod)], "[+]")[[1]]
  facs <- as.character(srMod)[2]
  facs <- gsub("(^ )|( $)", "", facs) # remove leading and trailing spaces
  grepl(paste("(^",a4aSRmodelList,"[(])", collapse = "|", sep = ""), facs)
}


geta4aSRmodel <- function(srMod) {
  #facs <- strsplit(as.character(srMod)[length(srMod)], "[+]")[[1]]
  facs <- as.character(srMod)[2]
  facs <- gsub("(^ )|( $)", "", facs) # remove leading and trailing spaces
  a4as <- grepl(paste("(^",a4aSRmodelList,"[(])", collapse = "|", sep = ""), facs)
  if (sum(a4as) > 1) stop("you can only specify one type of stock recruit relationship.")
  if (sum(a4as) == 0) "none()" else facs[a4as]
}
