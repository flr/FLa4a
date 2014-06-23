library(FLa4a)
data(ple4)
data(ple4.indices)
err <- 0.05
nits <- 10

fit <-  a4aSCA(ple4, ple4.indices[1], qmodel=list(~s(age, k=4)))
stk <- ple4 + fit
fits <- simulate(fit, nits, 1234)
stks <- propagate(ple4, nits) + fits 
idxs <- ple4.indices[1]
index(idxs[[1]]) <- index(fits)[[1]]

fits <-  a4aSCA(stks, idxs, qmodel=list(~s(age, k=4)))
stks2 <- stks + fits

# are the simulated values unbiased
stk.fit <- stock.n(fit)
stk.sim <- stock.n(stks2)
stk.rat <- iterMedians(stk.sim)/stk.fit
f.fit <- harvest(fit)
f.sim <- harvest(stks2)
f.rat <- iterMedians(f.sim)/f.fit
cth.fit <- catch.n(fit)
cth.sim <- catch.n(stks2)
cth.rat <- iterMedians(cth.sim)/cth.fit

max(stk.rat)-min(stk.rat) < err
max(f.rat)-min(f.rat) < err
max(cth.rat)-min(cth.rat) < err

#plot(FLStocks(orig=stk, sim=stks, fitsim=stks+fits))

