#' @title Get model matrix
#' @name getX
#' @rdname getX-methods
#' @description Uses the user-specified formula to build a model matrix.
#' @template bothargs
#' @param df the data.frame to build the model matrix against.
#' @return a matrix.
#' @note \code{getX} is intended to be used internally
#' @aliases getX getX-methods
setGeneric("getX", function(object, ...) standardGeneric("getX"))

#' @rdname getX-methods
setMethod("getX", "submodel",
  function(object) {
    # make empty FLQuant and make a data.frame from it
    flq <- flq_from_range(object)
    df <- as.data.frame(flq)

    # add in covariates
    if (length(object@covariates) > 0) {
      covar <- object@covariates
      # add in covariates to data.frame
      tmp <- lapply(seq_along(covar), function(i) {
        x <- as.data.frame(covar[[i]])[c("age", "year", "data")]
        if (length(unique(x$age)) == 1) x <- x[names(x) != "age"]
        if (length(unique(x$year)) == 1) x <- x[names(x) != "year"]
        names(x) <- gsub("data", names(covar)[i], names(x))
        x
      })
      covar.df <- tmp[[1]]
      for (i in seq_along(tmp[-1]))
        covar.df <- merge(covar.df, i, all = TRUE, sort = FALSE)

      df <- merge(df, covar.df, all.x = TRUE, all.y = FALSE)
    }

    # get design matrix
    getX(formula(object), df)
  }
)

#' @rdname getX-methods
setMethod("getX", "formula",
  function(object, df) {
    opts <-
      options(contrasts = c(unordered = "contr.sum", ordered = "contr.poly"))

    model <- object

    # drop unused factor levels
    facs <- which(sapply(df, is.factor))
    df[facs] <- lapply(df[facs], function(x) x[drop = TRUE])

    # quick fix for problems predicting with smooths...
    # this will fail in some instances when covariates are included,
    olddf <- df
    df <- unique(df)

    model.type <- deparse(substitute(model))

    # step 1 - separate out elements
    facs <- strsplit(as.character(model)[length(model)], "[+]")[[1]]
    facs <- gsub("(^ )|( $)", "", facs) # remove leading and trailing spaces

    # some 'model builder' functions.  Here for now, but maybe move them
    # somewhere else later.
    # note they use df through their scope - so moving might be troublesome
    # they just need to return a character version of the formula element
    # they code for
    trawl <- function(plateau, selectivity = "fixed", ...) {
      selectivity <- match.arg(selectivity, c("fixed", "variable"))

      # implement plateau and calculate appropriate degrees of freedom for age
      if (missing(plateau)) {
        ka <- ceiling(0.5 * length(unique(df$age)))
        var <- "age"
      } else {
        ka <-
          ceiling(0.5 *
                  length(unique(replace(df$age, df$age > plateau, plateau))))
        var <- paste("replace(age, age >", plateau, ",", plateau, ")")
      }

      # apply a fixed or evolving selectivity pattern and get some
      # appropriate degreed of freedom
      if (selectivity == "fixed") {
        ky <- ceiling(0.5 * length(unique(df$year)))
        out <- paste("s(", var, ", k =", ka, ") + s(year, k =", ky, ")")
      } else {
        ky <- ceiling(0.35 * length(unique(df$year)))
        out <- paste("te(", var, ",year, k = c(", ka, ",", ky, "))")
      }

      out
    }

    # some other a4a functions that operate on inputs directly
    breakpts <- function(var, breaks, ...) {
      if (min(var, na.rm = TRUE) < min(breaks))
        breaks <- c(min(var, na.rm = TRUE) - 1, breaks)
      if (max(var, na.rm = TRUE) > max(breaks))
        breaks <- c(breaks, max(var, na.rm = TRUE))
      label <- paste0("(", breaks[-length(breaks)], ",", breaks[-1], "]")
      cut(var, breaks = breaks, label = label)
    }

    # look for builder functions like trawl and substitute with their definition
    a4as <- grepl("(^trawl[(])", facs)
    if (any(a4as)) {
      facs[a4as] <- sapply(facs[a4as], function(x) eval(parse(text = x)))
    }
    # and split by "+" again
    model <- eval(parse(text = paste("~", paste(facs, collapse = "+"))))
    facs <- strsplit(as.character(model)[length(model)], "[+]")[[1]]
    facs <- gsub("(^ )|( $)", "", facs) # remove leading and trailing spaces


    # evaluate by argument in gam elements
    gams <- grepl("(^s[(])|(^te[(])|(^ti[(])", facs)
    if (any(gams)) {
      # function to return by argument data
      tmp.sfunc <- function(..., by = NULL) {
        out <- eval(substitute(by), df)
        if (is.logical(out)) {
          factor(out)
        } else {
          out
        }
      }
      dummy.gams <- gsub("(^s[(])|(^te[(])|(^ti[(])", "tmp.sfunc(", facs[gams])
      gmf <- lapply(dummy.gams, function(x) eval(parse(text = x)))
      bygams <- !sapply(gmf, is.null)
      if (any(bygams)) {
        # if there are by arguments then add them to df otherwise do nothing
        gmf <- gmf[bygams]
        gmf <- as.data.frame(gmf)
        names(gmf) <- paste("by", 1:ncol(gmf), sep = "")
        if (!all(complete.cases(gmf)))
          stop("Some 'by' arguments in smoothers evaluate to NA:",
               "cannot proceed.")
        df <- cbind(df, gmf)
        # now replace by = ... to by = by1. etc and rebuild the formula :)
        # function to get by argument NAME
        tmp.sfunc <- function(..., by = NULL) deparse(substitute(by))
        replace.by <-
          lapply(dummy.gams[bygams],
            function(x) {
              gsub("[)]", "[)]", gsub("[(]", "[(]", eval(parse(text = x))))
            })
        facs[gams][bygams] <-
          sapply(seq(ncol(gmf)),
            function(i) {
              gsub(replace.by[[i]], names(gmf)[i], facs[gams][bygams][i])
            })
        model <- eval(parse(text = paste("~", paste(facs, collapse = "+"))))
      }
    }

    # keep final model after processing?
    # Or provide formula processing as a seprate function?

    # check non gam covariates for NAs
    if (any(!gams)) {
      model_sansgam <-
        eval(parse(text = paste("~", paste(facs[!gams], collapse = "+"))))
        modelframe_sansgam <- try(
          model.frame(model_sansgam, df)
        )
      if (inherits(modelframe_sansgam, "try-error")) {
        stop("Please check your model formula: ",
             print(model_sansgam, showEnv = FALSE), ", ",
             "there seems to be an error.",
             call. = FALSE)
      }
      if (nrow(modelframe_sansgam) != nrow(df))
        stop("something went wrong - check for NAs in covariates",
             call. = FALSE)
    }

    # now run formula through model.matrix
    if (!any(gams)) {
      X <- model.matrix(model, df)
    } else {
      deparsed_model <- deparse(model[[length(model)]], width.cutoff = 500L)
      model <- formula(paste("fake.obs ~", deparsed_model))
      G <- try(
          gam(model, data = cbind(fake.obs = 1, df), fit = FALSE)
        )
      if (inherits(G, "try-error")) {
        stop("Please check your (gam) model formula: ",
             print(model[-2], showEnv = FALSE), ", ",
             "there seems to be an error.",
             call. = FALSE)
      }
      X <- G$X
      colnames(X) <- G$term.names
    }
    # a double check
    if (nrow(X) != nrow(df))
      stop("something went wrong - check for NAs in covariates")

    # check model for redundant parameters
    qr.X <- qr(X)
    rank.deficient <- qr.X$pivot[abs(diag(qr.X$qr)) < 1e-7]
    if (length(rank.deficient)) {
      droppar <- paste(colnames(X)[rank.deficient], collapse = "\n\t")
      warning("*** ", model.type, " has ",  length(rank.deficient),
              " too many parameter(s)!!\n",
              "    i will remove the redundant ones:\n\t",
              droppar, call. = FALSE)
      X <- X[, -rank.deficient]
    }

    # reset options
    options(opts)

    # now remap from df to olddf
    df$id <- 1:nrow(df)
    olddf$oldid <- 1:nrow(olddf)
    olddf <- merge(olddf, df, all = TRUE)
    olddf <- olddf[order(olddf$oldid), ]

    X[olddf$id,, drop = FALSE]
  }
)


#' @title Get covariance matrix
#' @name getCov
#' @rdname getCov-methods
#' @description Returns the covariance matrix of the specified Gaussian markov random field model.
#' @param n integer giving the size of the random feild
#' @param model chatacter giving the name of the GMRF
#' @param tau numeric giving the multiplier of the structure matrix for the model
#' @return a covariance matrix
#' @aliases getCov
getCov <- function(n, model, tau) {
  # in the future thinka about rw1, rw2 etc.
  model <- match.arg(model, "iid")

  as(tau * Matrix::Diagonal(n), "CsparseMatrix")
}
