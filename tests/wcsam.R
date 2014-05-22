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

lst <- lapply(split(1:nits, 1:nits), function(x){
	out <- try(a4aSCA(iter(stks, x), FLIndices(iter(idxs[[1]], x)), fit="MP")) 
	if(is(out, "try-error")) NULL else out
})

stks2 <- stks
for(i in 1:nits){
	iter(catch.n(stks2), i) <- catch.n(lst[[i]])
	iter(stock.n(stks2), i) <- stock.n(lst[[i]])
	iter(harvest(stks2), i) <- harvest(lst[[i]])
} 
catch(stks2) <- computeCatch(stks2) 
stock(stks2) <- computeStock(stks2) 

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

plot(FLStocks(orig=stk, sim=stks, fitsim=stks2))

