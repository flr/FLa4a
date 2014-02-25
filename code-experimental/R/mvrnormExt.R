setMethod("mvrnorm", signature("numeric", "a4aFitExt"),
  function(n = 1, mu, tol = 1e-9) 
  {
    # get the mean
    beta <- c(mu @ baseLvlPars)
    # number of params
    p <- length(beta)
    # simlulate iid normals (0,1)
    Z <- matrix(rnorm(p * n), p)
    # work out density using Z - the quick way
    ldens <-  -p/2 * log(2 * pi) + sum(log(diag(mu @ L))) - 0.5 * colSums(Z^2)
    # transform and rescale
    X <- drop(beta) + mu @ L %*% Z
    rownames(X) <- names(beta)
    attr(X, "logdens") <- ldens
    X
  }
)

