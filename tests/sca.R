#====================================================================
# tests for sca
#====================================================================
library(FLa4a)
data(ple4)
data(ple4.indices)
data(ple4.index)
nits <- 2

#====================================================================
# run sca
#====================================================================
fit0 <-  sca(ple4, FLIndices(ple4.index))

#--------------------------------------------------------------------
# iters
#--------------------------------------------------------------------

idx2 <- propagate(ple4.index, nits)
stk2 <- propagate(ple4, nits)

# Nx1
fit <- sca(stk2, FLIndices(ple4.index))
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

# 1xN
fit <- sca(ple4, FLIndices(idx2))
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

# NxN
fit <- sca(stk2, FLIndices(idx2))
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

#====================================================================
# run a4aSCA
#====================================================================
fit0 <-  a4aSCA(ple4, FLIndices(ple4.index), qmodel=list(~s(age, k=4)))

#--------------------------------------------------------------------
# iters
#--------------------------------------------------------------------

idx2 <- propagate(ple4.index, nits)
stk2 <- propagate(ple4, nits)

# Nx1
fit <- a4aSCA(stk2, FLIndices(ple4.index), qmodel=list(~s(age, k=4)))
dims(fit)$iter==nits
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

dim(fit@pars@stkmodel@vcov)[3]==nits
dim(fit@pars@stkmodel@params)[2]==nits
dim(fit@pars@qmodel[[1]]@vcov)[3]==nits
dim(fit@pars@qmodel[[1]]@params)[2]==nits
dim(fit@pars@vmodel[[1]]@vcov)[3]==nits
dim(fit@pars@vmodel[[1]]@params)[2]==nits
dim(fit@pars@vmodel[[2]]@vcov)[3]==nits
dim(fit@pars@vmodel[[2]]@params)[2]==nits

# 1xN
fit <- a4aSCA(ple4, FLIndices(idx2), qmodel=list(~s(age, k=4)))
dims(fit)$iter==nits
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

dim(fit@pars@stkmodel@vcov)[3]==nits
dim(fit@pars@stkmodel@params)[2]==nits
dim(fit@pars@qmodel[[1]]@vcov)[3]==nits
dim(fit@pars@qmodel[[1]]@params)[2]==nits
dim(fit@pars@vmodel[[1]]@vcov)[3]==nits
dim(fit@pars@vmodel[[1]]@params)[2]==nits
dim(fit@pars@vmodel[[2]]@vcov)[3]==nits
dim(fit@pars@vmodel[[2]]@params)[2]==nits

# NxN
fit <- a4aSCA(stk2, FLIndices(idx2), qmodel=list(~s(age, k=4)))
dims(fit)$iter==nits
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

dim(fit@pars@stkmodel@vcov)[3]==nits
dim(fit@pars@stkmodel@params)[2]==nits
dim(fit@pars@qmodel[[1]]@vcov)[3]==nits
dim(fit@pars@qmodel[[1]]@params)[2]==nits
dim(fit@pars@vmodel[[1]]@vcov)[3]==nits
dim(fit@pars@vmodel[[1]]@params)[2]==nits
dim(fit@pars@vmodel[[2]]@vcov)[3]==nits
dim(fit@pars@vmodel[[2]]@params)[2]==nits

#====================================================================
# run a4aSCA with simulate
#====================================================================
fit0 <-  a4aSCA(ple4, FLIndices(ple4.index), qmodel=list(~s(age, k=4)))

stk2 <- ple4
idx2 <- ple4.index
catch.n(stk2) <- genFLQuant(catch.n(stk2), 0.1, niter=nits)
index(idx2) <- genFLQuant(index(fit0)[[1]], 0.1, niter=nits)

# check iters match 
fit <- a4aSCA(stk2, FLIndices(idx2), qmodel=list(~s(age, k=4)))
fit0a <- a4aSCA(iter(stk2,1), FLIndices(iter(idx2,1)), qmodel=list(~s(age, k=4)))
fit0b <- a4aSCA(iter(stk2,2), FLIndices(iter(idx2,2)), qmodel=list(~s(age, k=4)))

identical(catch.n(fit)[,,,,,1], catch.n(fit0a))
identical(stock.n(fit)[,,,,,1], stock.n(fit0a))
identical(harvest(fit)[,,,,,1], harvest(fit0a))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0b)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0b)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0b)[drop=TRUE])

#====================================================================
# retrospective analysis
#====================================================================

# sca with defaults
fit0 <-  sca(ple4, FLIndices(ple4.index))
ret <- ra(ple4, FLIndices(ple4.index), 3)
stk <- window(ple4, end=range(ple4)["maxyear"]-2)
idx <- window(ple4.index, end=range(ple4)["maxyear"]-2)
fit3 <- sca(stk, FLIndices(idx), fit="assessment")
identical(stk + fit3, ret[[3]])

# a4aSCA with setted up models
fmod <- ~s(age, k=4)+s(year, k=20)
qmod <- list(~s(age, k=4))
n1mod <- ~s(age, k=4)
vmod <- list(~1, ~1)
srmod <- ~factor(year)

fit0 <-  a4aSCA(ple4, FLIndices(ple4.index), fmodel=fmod, qmodel=qmod, n1model=n1mod, vmodel=vmod, srmodel=srmod)
ret <- ra(ple4, FLIndices(ple4.index), 3, fmodel=fmod, qmodel=qmod, n1model=n1mod, vmodel=vmod, srmodel=srmod)
stk <- window(ple4, end=range(ple4)["maxyear"]-2)
idx <- window(ple4.index, end=range(ple4)["maxyear"]-2)
fit3 <- a4aSCA(stk, FLIndices(idx), fmodel=fmod, qmodel=qmod, n1model=n1mod, vmodel=vmod, srmodel=srmod)
identical(stk + fit3, ret[[3]])







































