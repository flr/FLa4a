#' creates a matrix operator for the nth difference
#'
#' @param n the length of the vector for difference are to be calculated
#' @return a matrix which functions as the difference operator
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
Dfunc <- function(n) {
  D <- diag(n)
  diag(D[-1,]) <- -1
  D[-1,]
}

#' Makes a precision matrix for a specific GMRF model
#'
#'
#' @param n the size of the GMRF
#' @param type a character giving the name of a GMRF (only rw1 and rw2 are currently available)
#' @param weights numeric to supply weights if GMRF is to have a different variance through the feild
#' @return the covariance matrix of a GMRF
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
Qfunc <- function(n, type = "rw1", weights = 1) {
  type <- match.arg(type, c("rw1", "rw2"))
  
  switch(type,
      rw1 = {
        if (!(length(weights) %in% c(1, n-1))) 
          stop("supplied weights should be of length n - 1 forRW1 model")
        D <- Dfunc(n)
        Q <- t(D) %*% (diag(weights, nrow = n-1) %*% D)  
      },
      rw2 = {
        if (!(length(weights) %in% c(1, n-2)))
          stop("supplied weights should be of length n - 2 for RW2 model")
        D <- Dfunc(n-1) %*% Dfunc(n)
        Q <- t(D) %*% (diag(weights, nrow = n-2) %*% D) 
      })
  Q
}

#' STILL EXPERIMENTAL provides the covariance of a composite GMRF model designed to produce 
#' families of selectivity curves such as flat top
#'
#'
#' @param n the size of the GMRF to prodce, usually the number of ages
#' @param type the model
#' @param edf NOT USED will be the expected degrees of freedom of the model
#' @return the covariance matrix of a GMRF
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
selQfunc <- function(n, type = "smooth", edf) {
  type <- match.arg(type, c("smooth", "flattop"))
  
  switch(type,
      smooth = {
        Q <- Qfunc(n, type = "rw2")  
      },
      flattop = {
        weights <- rep(0, n - 1)
        filt <- floor(n/2):(n-1)
        weights[filt] <-  exp( 5 * (seq_along(filt) - 1) /(length(filt - 1)) )
        Q <- Qfunc(n, type = "rw2", weights = weights[-1] + 1) + Qfunc(n, type = "rw1", weights = weights) 
      })
  Q * n * n
}

#' samples from a multivariate normal given a mean and precision
#'
#'
#' @param n the number of samples to take
#' @param mu the mean vector
#' @param Q the precision matrix
#' @param tol the tolerance for what is positive definate
#' @param return.density whether to return the log density of the simulation.  Only makes sence when n is 1.
#' @return a pointer to the environment in which summaries of the data reside
#' @author Colin Millar \email{colin.millar@@jrc.ec.europa.eu}
mvrnormQ <- function (n, mu, Q, tol = 1e-06, return.density = FALSE) {
  # 
  #  Sample from a GMRF with mean mu and precision Q
  #  based on mvrnorm from library MASS
  #
  p <- length(mu)
  if (!all(dim(Q) == c(p, p))) stop("incompatible arguments")
  eQ <- eigen(Q, symmetric = TRUE, EISPACK = TRUE)
  ev <- eQ $ values
  if (!all(ev >= -tol * abs(ev[1L])))
    stop("'Q' is not positive definite")
  X <- matrix(rnorm(p * n), nrow = n)
  evinv <- ifelse(ev < tol, 0, 1/ev)
  
  X <- drop(mu) + (eQ $ vectors %*% diag(sqrt(evinv), p)) %*% t(X)

  if (return.density & n == 1) {
    ldet <- sum(log(ev[ev>0]))
    x <- X - mu
    X <- list(x = X, log.dens = 0.5 * ldet - 0.5 * t(x) %*% (Q%*%x))
  }
  
  X
}

