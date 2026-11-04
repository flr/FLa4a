library(FLa4a)
data(ple4)
data(ple4.indices)
ple4.indices <- ple4.indices[c(3,4,5)]
err <- 0.1
nits <- 250

fit <-  sca(ple4, ple4.indices[1], qmodel=list(~s(age, k=4)))
fit.sim <- simulate(fit, nits, 1234)

idx.sim <- ple4.indices[1]
index(idx.sim[[1]]) <- index(fit.sim)[[1]]
stk.sim <- ple4 + fit.sim

fit2 <-  sca(stk.sim, idx.sim, qmodel=list(~s(age, k=4)), fit="MP")

# are the simulated values unbiased
stk.fit <- stock.n(fit)
stk.q <- apply(stock.n(fit2), c(1:5), quantile, probs=c(err/2, 1-err/2))
f.fit <- harvest(fit)
f.q <- apply(harvest(fit2), c(1:5), quantile, probs=c(err/2, 1-err/2))
cth.fit <- catch.n(fit)
cth.q <- apply(catch.n(fit2), c(1:5), quantile, probs=c(err/2, 1-err/2))

mean(iter(stk.q, 1) < stk.fit & iter(stk.q,2) > stk.fit)==1
mean(iter(f.q, 1) < f.fit & iter(f.q,2) > f.fit)==1
mean(iter(cth.q, 1) < cth.fit & iter(cth.q,2) > cth.fit)==1
