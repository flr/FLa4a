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

# check
v0 <- ssb(stk00)


