# tests for observation error

library(FLa4a)
library(ggplotFL)
data(ple4)
data(ple4.indices)
nsim <- 250
fit <- sca(ple4, ple4.indices)
res <- residuals(fit, ple4, ple4.indices, type="deviances")
# simulate
stk01 <- ple4 + simulate(fit, nsim, seed=1234)
# add oe to catches and indices using vmodel estimates
fits <- simulate(fit, nsim, seed=1234, obserror=TRUE)
stk02 <- ple4 + fits
idx02 <- ple4.indices + fits
# add oe to catches using residual variance
stk03 <- propagate(ple4, nsim) + res
# add oe to catches using residual variance
stk04 <- stk01 + res

# plot
stks <- FLStocks(s01 = stk01, s02 = stk02, s03 = stk03, s04 = stk04)
cths <- lapply(stks, "slot", "catch.n")

system.time(fit02 <- sca(stk02, idx02))
