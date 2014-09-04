library(FLa4a)
data(ple4)
data(ple4.indices)
err <- 0.05
nits <- 250

fit <-  a4aSCA(ple4, ple4.indices[1], qmodel=list(~s(age, k=4)))

fit.sim <- simulate(fit, nits, 1234)

idx.sim <- ple4.indices[1]
index(idx.sim[[1]]) <- index(fit.sim)[[1]]
stk.sim <- propagate(ple4, nits) + fit.sim 

fit2 <-  a4aSCA(stk.sim, idx.sim, qmodel=list(~s(age, k=4)), fit="MP")

# are the simulated values unbiased
stk.fit <- stock.n(fit)
stk.sim <- stock.n(fit2)
stk.rat <- iterMedians(stk.sim)/stk.fit
f.fit <- harvest(fit)
f.sim <- harvest(fit2)
f.rat <- iterMedians(f.sim)/f.fit
cth.fit <- catch.n(fit)
cth.sim <- catch.n(fit2)
cth.rat <- iterMedians(cth.sim)/cth.fit

max(stk.rat)-min(stk.rat) < err
max(f.rat)-min(f.rat) < err
max(cth.rat)-min(cth.rat) < err



