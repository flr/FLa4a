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
  D <- gmrf::Drw(nx, order = m)

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



smooth.construct.bs_knots.smooth.spec <- function(object, data, knots) {

  if (!is.null(object$xt$gaps)) {

    # spread knots accross data range avoiding the gaps
    x <- unique(data[[object$term]])
    xg <- x[!x %in% object$xt$gaps]

    k <- getknots(x, object$p.order[1], object$bs.dim, xg)

    knots <- list(k)
    names(knots) <- object$term
  } else {
    knots <- object$xt$knots
  }

  mgcv::smooth.construct.bs.smooth.spec(object, data, knots)
}


getknots <- function(x, m, bs.dim, available_x = x) {
  # number of required knots
  nk <- bs.dim - m + 1

  # what is the average knot distance
  xl <- min(x)
  xu <- max(x)
  xr <- xu - xl
  xl <- xl - xr * 0.001
  xu <- xu + xr * 0.001
  dx <- (xu - xl) / (nk - 1)

  # outer knots are easy:
  k_outer_lower <- xl - (m:1) * dx
  k_outer_upper <- xu + (1:m) * dx

  ## get knots based in distribution of available data
  xdens <- density(available_x, bw = dx / 2)
  # restrict to data range
  xdens_which <- which(xdens$x >=xl & xdens$x <= xu)
  xcumdens <- cumsum(xdens$y[xdens_which]) / sum(xdens$y[xdens_which])

  k_inner <-
    approx(x = c(0, xcumdens),
           y = c(xl, xdens$x[xdens_which]),
           xout = seq(0, 1, length = nk))$y
  k_inner[1] <- xl
  k_inner[nk] <- xu

  # construct all knots
  c(k_outer_lower, k_inner, k_outer_upper)
}
