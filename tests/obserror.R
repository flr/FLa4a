# tests for observation error

library(FLa4a)
data(ple4)
data(ple4.indices)
nsim <- 50
fit <- sca(ple4, ple4.indices)
res <- residuals(fit, ple4, ple4.indices, type="deviances")
stk00 <- ple4 + fit
# simulate
stk01 <- ple4 + simulate(fit, nsim, seed=1234)
# add oe to catches and indices using vmodel estimates
fits <- simulate(fit, nsim, seed=1234, obserror=TRUE)
stk02 <- ple4 + fits
idx02 <- ple4.indices + fits
fit02 <- sca(stk02, idx02)
stk02b <- stk02 + fit02

# check
v0 <- ssb(stk00)
ul <- apply(ssb(stk02b), 1:5, quantile, p=0.975)
ll <- apply(ssb(stk02b), 1:5, quantile, p=0.025)
v0 <- v0 < ul & v0 > ll
sum(v0) == length(v0)

v0 <- rec(stk00)
ul <- apply(rec(stk02b), 1:5, quantile, p=0.975)
ll <- apply(rec(stk02b), 1:5, quantile, p=0.025)
v0 <- v0 < ul & v0 > ll
sum(v0) == length(v0)

v0 <- fbar(stk00)
ul <- apply(fbar(stk02b), 1:5, quantile, p=0.975)
ll <- apply(fbar(stk02b), 1:5, quantile, p=0.025)
v0 <- v0 < ul & v0 > ll
sum(v0) == length(v0)

v0 <- catch(stk00)
ul <- apply(catch(stk02b), 1:5, quantile, p=0.975)
ll <- apply(catch(stk02b), 1:5, quantile, p=0.025)
v0 <- v0 < ul & v0 > ll
sum(v0) == length(v0)
