#' @export
smooth.construct.cgmrf.smooth.spec <- function (object, data, knots) {

  # the order of the gmrf
  m <- object$p.order[1]

  # get unique x locations
  if (is.null(object$xt$xrange)) {
    x <- sort(unique(data[[object$term]]))
  } else {
    x <- object$xt$xrange[1]:object$xt$xrange[2]
  }

  # size of full GMRF
  nx <- length(x)

  # contruct an m_th order distance matrix
  D <- Drw(nx, order = m)

  # set up the step variances (nx - m steps)
  weights <- rep(1, nx - m)
  # make the 'gap' years stiffer / less bendy
  if (is.null(object$xt$weight)) {
    # set a sensible default?
    ## object$xt$weight <- 0.5
    # or go for a hard penalty of contant (m=1) or
    # linear (m=2) or quadratic (m=3) interpolation
    object$xt$weight <- 0.01
  }
  # apply weight restriction so that we link missing years to
  # present years.  This means that m = 2 will link back
  # two years to allow the slope to be inferered,
  # and for m = 3 three years to get at the quadratic
  g <- object$xt$gaps
  for (i in seq_len(m)) {
    g <- c(g, unique(object$xt$gaps) - i)
  }
  g <- sort(unique(g))
  weights[which(x %in% g)] <- object$xt$weight
  weights <- 1/weights^2

  Q <- (t(D) %*% diag(weights)) %*% D
  eQ <- eigen(Q)
  X <- eQ$vectors

  # add row names for Predict method
  row.names(X) <- paste(x)
  # scale X !!
  #X <- sweep(X, 2, rev(eQ$values)[-1])

  # store X in xt for the Predict method
  object$xt$X <- X

  # set up X in object so it matches the data
  object$X <- Predict.matrix.cgmrf.smooth(object, data)

  # set penalty matrix
  object$S <- list(diag(object$bs.dim))

  # set smoother properties
  object$rank <- object$bs.dim
  object$null.space.dim <- 0
  object$C <- matrix(0, 0, ncol(object$X))
  object$side.constrain <- FALSE
  object$plot.me <- FALSE
  object$te.ok <- 2
  object$noterp <- TRUE
  object$random <- TRUE
  class(object) <- "cgmrf.smooth"
  object
}

#' @export
Predict.matrix.cgmrf.smooth <- function (object, data)
{
  X <- object$xt$X[,ncol(object$xt$X) - 1:(object$bs.dim-1)]
  X[paste(data[[object$term]]),]
}


# taken from colins gmrf package - once on CRAN will link directly
Drw <- function (n, order = 2, cyclic = FALSE)
{
    stopifnot(order > 0)
    out <- diag(n)
    for (i in 1:order) {
        out <- out %*% Drw1(n, cyclic = TRUE)
    }
    if (cyclic)
        out
    else out[1:(n - order), ]
}

Drw1 <- function (n, cyclic = FALSE)
{
    out <- diag(n)
    diag(out[, -1]) <- -1
    out[n, 1] <- -1
    if (cyclic)
        out
    else out[-n, ]
}
