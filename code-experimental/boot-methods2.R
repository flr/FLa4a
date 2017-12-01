library(FLa4a)
data(ple4)
data(ple4.indices)
# number of resamples
it <- 250
# trimming indices to simplify example
ple4.indices <- ple4.indices[1]
fit <- sca(ple4, ple4.indices)
res <- residuals(fit, ple4, ple4.indices)
# standard deviation of residuals
# we'll use a parametric bootstrap instead of resampling the residuals
# because you don't have lots of datapoints
sdres <- lapply(res, sd, na.rm=TRUE)
# control de seed for replicability
set.seed(1234)
# simulate a lognormal distribution with median equal to estimates
catch.n(ple4) <- rlnorm(it, log(catch.n(fit)), sdres[[1]])
# control de seed for replicability
set.seed(1234)
# simulate a lognormal distribution with median equal to estimates
index(ple4.indices[[1]]) <- rlnorm(it, log(index(fit)[[1]]), sdres[[2]])
# refit
fitboot <- sca(ple4, ple4.indices)
# magic :)
plot(ple4+fitboot)
# medians should be very similar
iterMedians(catch.n(fitboot))/catch.n(fit)
iterMedians(harvest(fitboot))/harvest(fit)
iterMedians(stock.n(fitboot))/stock.n(fit)



