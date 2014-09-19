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

# check that indices have attr biomass
"FLIndexBiomass" %in%  names(attributes(index(fit0)[[1]]))
"range" %in%  names(attributes(index(fit0)[[1]]))

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
# check that indices have attr biomass
"FLIndexBiomass" %in%  names(attributes(index(fit0)[[1]]))
"range" %in%  names(attributes(index(fit0)[[1]]))

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
# bug in name matching between idxs names and pars names
#====================================================================
fit0 <-  a4aSCA(ple4, FLIndices(a=ple4.index), qmodel=list(~s(age, k=4)))
# this bug concatenated catch pars with qpars
length(params(vmodel(pars(fit0))[[2]]))==1

#====================================================================
# hessian non-positive definite
#====================================================================

fit0 <- FLa4a:::a4aInternal(ple4, FLIndices(ple4.index), fmodel=~factor(age)+ factor(year), qmodel=list(~factor(age)+ s(year, k=20)))
sum(stock.n(fit0), na.rm=T)==0
sum(catch.n(fit0), na.rm=T)==0
sum(index(fit0)[[1]], na.rm=T)==0
sum(pars(fit0)@stkmodel@params, na.rm=T)==0
sum(pars(fit0)@qmodel[[1]]@params, na.rm=T)==0
sum(pars(fit0)@vmodel[[1]]@params, na.rm=T)==0
sum(pars(fit0)@vmodel[[2]]@params, na.rm=T)==0
sum(pars(fit0)@stkmodel@vcov, na.rm=T)==0
sum(pars(fit0)@qmodel[[1]]@vcov, na.rm=T)==0
sum(pars(fit0)@vmodel[[1]]@vcov, na.rm=T)==0
sum(pars(fit0)@vmodel[[2]]@vcov, na.rm=T)==0

fit <- a4aSCA(ple4, FLIndices(ple4.index), fmodel=~factor(age)+ factor(year), qmodel=list(~factor(age)+ s(year, k=20)))
sum(stock.n(fit), na.rm=T)==0
sum(catch.n(fit), na.rm=T)==0
sum(index(fit)[[1]], na.rm=T)==0
sum(pars(fit)@stkmodel@params, na.rm=T)==0
sum(pars(fit)@qmodel[[1]]@params, na.rm=T)==0
sum(pars(fit)@vmodel[[1]]@params, na.rm=T)==0
sum(pars(fit)@vmodel[[2]]@params, na.rm=T)==0
sum(pars(fit)@stkmodel@vcov, na.rm=T)==0
sum(pars(fit)@qmodel[[1]]@vcov, na.rm=T)==0
sum(pars(fit)@vmodel[[1]]@vcov, na.rm=T)==0
sum(pars(fit)@vmodel[[2]]@vcov, na.rm=T)==0

fit1 <- sca(ple4, FLIndices(ple4.index), fmodel=~factor(age)+ factor(year), qmodel=list(~factor(age)+ s(year, k=20)), fit="assessment")
sum(stock.n(fit1), na.rm=T)==0
sum(catch.n(fit1), na.rm=T)==0
sum(index(fit1)[[1]], na.rm=T)==0
sum(pars(fit1)@stkmodel@params, na.rm=T)==0
sum(pars(fit1)@qmodel[[1]]@params, na.rm=T)==0
sum(pars(fit1)@vmodel[[1]]@params, na.rm=T)==0
sum(pars(fit1)@vmodel[[2]]@params, na.rm=T)==0
sum(pars(fit1)@stkmodel@vcov, na.rm=T)==0
sum(pars(fit1)@qmodel[[1]]@vcov, na.rm=T)==0
sum(pars(fit1)@vmodel[[1]]@vcov, na.rm=T)==0
sum(pars(fit1)@vmodel[[2]]@vcov, na.rm=T)==0

#====================================================================
# run sca with recruitment index
#====================================================================
fit0 <-  sca(ple4, FLIndices(ple4.index[1]), qmodel=list(~1))

#--------------------------------------------------------------------
# iters
#--------------------------------------------------------------------

idx2 <- propagate(ple4.index[1], nits)
stk2 <- propagate(ple4, nits)

# Nx1
fit <- sca(stk2, FLIndices(ple4.index[1]), qmodel=list(~1))
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

# 1xN
fit <- sca(ple4, FLIndices(idx2), qmodel=list(~1))
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

# NxN
fit <- sca(stk2, FLIndices(idx2), qmodel=list(~1))
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

#====================================================================
# run sca with FLQuantDistr
#====================================================================
data(ple4)
data(ple4.index)
catch.n(ple4) <- FLQuantDistr(catch.n(ple4), (0.2/catch.n(ple4))^2)
index.var(ple4.index) <- (0.2/index(ple4.index))^2

fit0 <-  sca(ple4, FLIndices(ple4.index), qmodel=list(~s(age, k=4)))

#--------------------------------------------------------------------
# iters
#--------------------------------------------------------------------

idx2 <- propagate(ple4.index, nits)
stk2 <- propagate(ple4, nits)
# is propagate working 
is(catch.n(stk2), "FLQuantDistr")
identical(c(catch.n(stk2)[,,,,,1]), c(catch.n(stk2)[,,,,,2]))
identical(c(var(catch.n(stk2)[,,,,,1])), c(var(catch.n(stk2)[,,,,,2])))
identical(c(catch.n(stk2)[,,,,,1]), c(catch.n(ple4)))
identical(c(var(catch.n(stk2)[,,,,,1])), c(var(catch.n(ple4))))

# Nx1
fit <- sca(stk2, FLIndices(ple4.index), qmodel=list(~s(age, k=4)))
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

# 1xN
fit <- sca(ple4, FLIndices(idx2), qmodel=list(~s(age, k=4)))
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

# NxN
fit <- sca(stk2, FLIndices(idx2), qmodel=list(~s(age, k=4)))
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

#--------------------------------------------------------------------
# zeros in var (must fail with error message)
#--------------------------------------------------------------------

var(catch.n(ple4))[4,10] <- 0
fit0 <-  try(sca(ple4, FLIndices(ple4.index), qmodel=list(~s(age, k=4))))
is(fit0, "try-error")
var(catch.n(ple4))[4,10] <- 0.1

index.var(ple4.index)[4,10] <- 0
fit0 <-  try(sca(ple4, FLIndices(ple4.index), qmodel=list(~s(age, k=4))))
is(fit0, "try-error")

#--------------------------------------------------------------------
# reset data
#--------------------------------------------------------------------
data(ple4)
data(ple4.index)

#====================================================================
# run a4aSCA with biomass index
#====================================================================
data(ple4)
bioidx <- FLIndexBiomass(FLQuant(NA, dimnames=list(age="all", year=range(ple4)["minyear"]:range(ple4)["maxyear"])), name="bioidx")
index(bioidx) <- stock(ple4)*0.001
index(bioidx) <- index(bioidx)*exp(rnorm(index(bioidx), sd=0.1))
range(bioidx)[c("startf","endf")] <- c(0,0)

# fitting the model
fit0 <- a4aSCA(ple4, FLIndices(bioidx), qmodel=list(~1))

#--------------------------------------------------------------------
# iters
#--------------------------------------------------------------------

idx2 <- propagate(bioidx, nits)
stk2 <- propagate(ple4, nits)

# Nx1
fit <- a4aSCA(stk2, FLIndices(bioidx), qmodel=list(~1))
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

# 1xN
fit <- a4aSCA(ple4, FLIndices(idx2), qmodel=list(~1))
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

# NxN
fit <- a4aSCA(stk2, FLIndices(idx2), qmodel=list(~1))
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

#====================================================================
# run a4aSCA with biomass index in specific ages
#====================================================================
bioidx <- FLIndexBiomass(FLQuant(NA, dimnames=list(age="all", year=range(ple4)["minyear"]:range(ple4)["maxyear"])), name="bioidx")
index(bioidx) <- stock(ple4)*0.001
index(bioidx) <- index(bioidx)*exp(rnorm(index(bioidx), sd=0.1))
range(bioidx)[c("min","max","startf","endf")] <- c(2,5,0,0)

# fitting the model
fit1 <- a4aSCA(ple4, FLIndices(bioidx), qmodel=list(~1))

# should be different from fit0
identical(fit1,fit0)

# set fit1 as fit0 to avoid changing all the code
fit0 <- fit1

#--------------------------------------------------------------------
# iters
#--------------------------------------------------------------------

idx2 <- propagate(bioidx, nits)
stk2 <- propagate(ple4, nits)

# Nx1
fit <- a4aSCA(stk2, FLIndices(bioidx), qmodel=list(~1))
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

# 1xN
fit <- a4aSCA(ple4, FLIndices(idx2), qmodel=list(~1))
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

# NxN
fit <- a4aSCA(stk2, FLIndices(idx2), qmodel=list(~1))
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

#====================================================================
# bug in indices naming 
#====================================================================
data(ple4)
data(ple4.index)
biofull <- 0.001*stock(ple4)
biofull <- FLIndexBiomass(index=biofull)
range(biofull)[c("startf","endf")] <- c(0,0)

fit0 <- sca(ple4, FLIndices(ple4.index, biofull), qmodel=list(~s(age, k=4), ~1))
fit1 <- sca(ple4, FLIndices(biofull, ple4.index), qmodel=list(~1, ~s(age, k=4)))

identical(stock.n(fit0), stock.n(fit1))
identical(harvest(fit0), harvest(fit1))


