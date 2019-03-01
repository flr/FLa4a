#====================================================================
# tests for sca
#====================================================================
library(FLa4a)
data(ple4)
data(ple4.indices)
ple4.indices <- ple4.indices[c(3,4,5)]
data(ple4.index)
nits <- 2

#====================================================================
# run sca MP
#====================================================================
fit0 <-  sca(ple4, FLIndices(ple4.index))
is(fit0, "a4aFit")

# check that indices have attr biomass, set to FALSE
"FLIndexBiomass" %in%  names(attributes(index(fit0)[[1]]))
!attr(index(fit0)[[1]], "FLIndexBiomass")

# check that indices have attr range
"range" %in%  names(attributes(index(fit0)[[1]]))

# check convergence info
fitSumm(fit0)["convergence",]==0

#--------------------------------------------------------------------
# iters
#--------------------------------------------------------------------

idx2 <- propagate(ple4.index, nits)
stk2 <- propagate(ple4, nits)

# Nx1
fit <- sca(stk2, FLIndices(ple4.index), fit = "MP")
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

# 1xN
fit <- sca(ple4, FLIndices(idx2), fit = "MP")
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

# NxN
fit <- sca(stk2, FLIndices(idx2), fit = "MP")
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

#--------------------------------------------------------------------
# check qmodel defaults
#--------------------------------------------------------------------
fit <- sca(ple4, FLIndices(ple4.index[1]))
all.equal(formula(qmodel(pars(fit))[[1]]), formula("~1"))
fit <- sca(ple4, FLIndices(ple4.index[1:2]))
all.equal(formula(qmodel(pars(fit))[[1]]), formula("~age"))
fit <- sca(ple4, FLIndices(ple4.index[1:3]))
all.equal(formula(qmodel(pars(fit))[[1]]), formula("~age"))
fit <- sca(ple4, FLIndices(ple4.index[1:4]))
all.equal(formula(qmodel(pars(fit))[[1]]), formula("~s(age, k=3)"))
fit <- sca(ple4, FLIndices(ple4.index))
all.equal(formula(qmodel(pars(fit))[[1]]), formula("~s(age, k=6)"))

#====================================================================
# run sca assessment
#====================================================================

# run
fit0 <-  sca(ple4, FLIndices(ple4.index), qmodel=list(~s(age, k=4)))
is(fit0, "a4aFitSA")

# check that indices have attr biomass, set to FALSE
"FLIndexBiomass" %in%  names(attributes(index(fit0)[[1]]))
!attr(index(fit0)[[1]], "FLIndexBiomass")

# check that indices have attr range
"range" %in%  names(attributes(index(fit0)[[1]]))

# check convergence info
fitSumm(fit0)["convergence",]==0

#--------------------------------------------------------------------
# iters
#--------------------------------------------------------------------

idx2 <- propagate(ple4.index, nits)
stk2 <- propagate(ple4, nits)

# Nx1
fit <- sca(stk2, FLIndices(ple4.index), qmodel=list(~s(age, k=4)), fit = "assessment")
dims(fit)$iter==nits
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

dim(vcov(stkmodel(pars(fit))))[3]==nits
dim(coef(stkmodel(pars(fit))))[2]==nits
dim(vcov(qmodel(pars(fit))[[1]]))[3]==nits
dim(coef(qmodel(pars(fit))[[1]]))[2]==nits
dim(vcov(vmodel(pars(fit))[[1]]))[3]==nits
dim(coef(vmodel(pars(fit))[[1]]))[2]==nits
dim(vcov(vmodel(pars(fit))[[2]]))[3]==nits
dim(coef(vmodel(pars(fit))[[2]]))[2]==nits


# 1xN
fit <- sca(ple4, FLIndices(idx2), qmodel=list(~s(age, k=4)), fit = "assessment")
dims(fit)$iter==nits
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

dim(vcov(stkmodel(pars(fit))))[3]==nits
dim(coef(stkmodel(pars(fit))))[2]==nits
dim(vcov(qmodel(pars(fit))[[1]]))[3]==nits
dim(coef(qmodel(pars(fit))[[1]]))[2]==nits
dim(vcov(vmodel(pars(fit))[[1]]))[3]==nits
dim(coef(vmodel(pars(fit))[[1]]))[2]==nits
dim(vcov(vmodel(pars(fit))[[2]]))[3]==nits
dim(coef(vmodel(pars(fit))[[2]]))[2]==nits


# NxN
fit <- sca(stk2, FLIndices(idx2), qmodel=list(~s(age, k=4)), fit = "assessment")
dims(fit)$iter==nits
dim(fitSumm(fit))[2]==nits
identical(catch.n(fit)[,,,,,1], catch.n(fit0))
identical(stock.n(fit)[,,,,,1], stock.n(fit0))
identical(harvest(fit)[,,,,,1], harvest(fit0))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0)[drop=TRUE])

dim(vcov(stkmodel(pars(fit))))[3]==nits
dim(coef(stkmodel(pars(fit))))[2]==nits
dim(vcov(qmodel(pars(fit))[[1]]))[3]==nits
dim(coef(qmodel(pars(fit))[[1]]))[2]==nits
dim(vcov(vmodel(pars(fit))[[1]]))[3]==nits
dim(coef(vmodel(pars(fit))[[1]]))[2]==nits
dim(vcov(vmodel(pars(fit))[[2]]))[3]==nits
dim(coef(vmodel(pars(fit))[[2]]))[2]==nits

#--------------------------------------------------------------------
# equal submodels check
#--------------------------------------------------------------------

# sca defaults
fit0 <-  sca(ple4, ple4.indices)
# default fit is "assessment" class should be "a4aFitSA"
is(fit0, "a4aFitSA")

# when fit="MP" class is "a4aFit"
fit0 <-  sca(ple4, ple4.indices, fit="MP")
is(fit0, "a4aFit")

# when fit="assessment" class is "a4aFitSA"
fit1 <-  sca(ple4, ple4.indices, fit="assessment")
is(fit1, "a4aFitSA")

# both must have the same results
all.equal(fitSumm(fit0)[-6,], fitSumm(fit1)[-6,])
all.equal(harvest(fit0), harvest(fit1))
all.equal(stock.n(fit0), stock.n(fit1))
all.equal(catch.n(fit0), catch.n(fit1))

#====================================================================
# run sca with simulate
#====================================================================
fit0 <-  sca(ple4, FLIndices(ple4.index), qmodel=list(~s(age, k=4)))

stk2 <- ple4
idx2 <- ple4.index
catch.n(stk2) <- genFLQuant(catch.n(stk2), 0.1, niter=nits)
index(idx2) <- genFLQuant(index(fit0)[[1]], 0.1, niter=nits)

# check iters match
fit <- sca(stk2, FLIndices(idx2), qmodel=list(~s(age, k=4)))
fit0a <- sca(iter(stk2,1), FLIndices(iter(idx2,1)), qmodel=list(~s(age, k=4)))
fit0b <- sca(iter(stk2,2), FLIndices(iter(idx2,2)), qmodel=list(~s(age, k=4)))

identical(catch.n(fit)[,,,,,1], catch.n(fit0a))
identical(stock.n(fit)[,,,,,1], stock.n(fit0a))
identical(harvest(fit)[,,,,,1], harvest(fit0a))
identical(catch.n(fit)[,,,,,2, drop=TRUE], catch.n(fit0b)[drop=TRUE])
identical(stock.n(fit)[,,,,,2, drop=TRUE], stock.n(fit0b)[drop=TRUE])
identical(harvest(fit)[,,,,,2, drop=TRUE], harvest(fit0b)[drop=TRUE])

#====================================================================
# bug in name matching between idxs names and pars names
#====================================================================
fit0 <-  sca(ple4, FLIndices(a=ple4.index), qmodel=list(~s(age, k=4)))
# this bug concatenated catch pars with qpars
length(params(vmodel(pars(fit0))[[2]]))==1

#====================================================================
# bug in gcv with NA
#====================================================================
stk <- ple4
catch.n(stk)[,"2000"] <- NA
fit <- sca(stk, ple4.indices)
!is.na(fitSumm(fit)["gcv",])

#====================================================================
# hessian non-positive definite
#====================================================================

fit0 <- FLa4a:::a4aInternal(ple4, FLIndices(ple4.index), fmodel=~factor(age)+ factor(year), qmodel=list(~factor(age)+ s(year, k=20)))
sum(stock.n(fit0), na.rm=T)==0
sum(catch.n(fit0), na.rm=T)==0
sum(index(fit0)[[1]], na.rm=T)==0
sum(params(stkmodel(pars(fit0))), na.rm=T)==0
sum(params(qmodel(pars(fit0))[[1]]), na.rm=T)==0
sum(params(vmodel(pars(fit0))[[1]]), na.rm=T)==0
sum(params(vmodel(pars(fit0))[[2]]), na.rm=T)==0
sum(vcov(stkmodel(pars(fit0))), na.rm=T)==0
sum(vcov(qmodel(pars(fit0))[[1]]), na.rm=T)==0
sum(vcov(vmodel(pars(fit0))[[1]]), na.rm=T)==0
sum(vcov(vmodel(pars(fit0))[[2]]), na.rm=T)==0


fit <- sca(ple4, FLIndices(ple4.index), fmodel=~factor(age)+ factor(year), qmodel=list(~factor(age)+ s(year, k=20)))
# check convergence info
fitSumm(fit)["convergence",]==1
sum(stock.n(fit), na.rm=T)==0
sum(catch.n(fit), na.rm=T)==0
sum(index(fit)[[1]], na.rm=T)==0
sum(params(stkmodel(pars(fit))), na.rm=T)==0
sum(params(qmodel(pars(fit))[[1]]), na.rm=T)==0
sum(params(vmodel(pars(fit))[[1]]), na.rm=T)==0
sum(params(vmodel(pars(fit))[[2]]), na.rm=T)==0
sum(vcov(stkmodel(pars(fit))), na.rm=T)==0
sum(vcov(qmodel(pars(fit))[[1]]), na.rm=T)==0
sum(vcov(vmodel(pars(fit))[[1]]), na.rm=T)==0
sum(vcov(vmodel(pars(fit))[[2]]), na.rm=T)==0

fit1 <- sca(ple4, FLIndices(ple4.index), fmodel=~factor(age)+ factor(year), qmodel=list(~factor(age)+ s(year, k=20)), fit="assessment")
sum(stock.n(fit1), na.rm=T)==0
sum(catch.n(fit1), na.rm=T)==0
sum(index(fit1)[[1]], na.rm=T)==0
sum(params(stkmodel(pars(fit1))), na.rm=T)==0
sum(params(qmodel(pars(fit1))[[1]]), na.rm=T)==0
sum(params(vmodel(pars(fit1))[[1]]), na.rm=T)==0
sum(params(vmodel(pars(fit1))[[2]]), na.rm=T)==0
sum(vcov(stkmodel(pars(fit1))), na.rm=T)==0
sum(vcov(qmodel(pars(fit1))[[1]]), na.rm=T)==0
sum(vcov(vmodel(pars(fit1))[[1]]), na.rm=T)==0
sum(vcov(vmodel(pars(fit1))[[2]]), na.rm=T)==0

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
catch.n(ple4) <- FLQuantDistr(catch.n(ple4), catch.n(ple4))
var(catch.n(ple4))[] <- 0.2
index.var(ple4.index)[] <- 0.2

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
# run sca with biomass index
#====================================================================
data(ple4)
bioidx <- FLIndexBiomass(FLQuant(NA, dimnames=list(age="all", year=range(ple4)["minyear"]:range(ple4)["maxyear"])), name="bioidx")
index(bioidx) <- stock(ple4)*0.001
index(bioidx) <- index(bioidx)*exp(rnorm(index(bioidx), sd=0.1))
range(bioidx)[c("startf","endf")] <- c(0,0)

# fitting the model
fit0 <- sca(ple4, FLIndices(bioidx), qmodel=list(~1))

#--------------------------------------------------------------------
# iters
#--------------------------------------------------------------------

idx2 <- propagate(bioidx, nits)
stk2 <- propagate(ple4, nits)

# Nx1
fit <- sca(stk2, FLIndices(bioidx), qmodel=list(~1))
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
# run sca with biomass index in specific ages
#====================================================================
bioidx <- FLIndexBiomass(FLQuant(NA, dimnames=list(age="all", year=range(ple4)["minyear"]:range(ple4)["maxyear"])), name="bioidx")
index(bioidx) <- stock(ple4)*0.001
index(bioidx) <- index(bioidx)*exp(rnorm(index(bioidx), sd=0.1))
range(bioidx)[c("min","max","startf","endf")] <- c(2,5,0,0)

# fitting the model
fit1 <- sca(ple4, FLIndices(bioidx), qmodel=list(~1))

# fit1 should be different from fit0
!identical(fit1,fit0)

# set fit0 as fit1 to avoid changing all the code
fit0 <- fit1

#--------------------------------------------------------------------
# iters
#--------------------------------------------------------------------

idx2 <- propagate(bioidx, nits)
stk2 <- propagate(ple4, nits)

# Nx1
fit <- sca(stk2, FLIndices(bioidx), qmodel=list(~1))
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
# bug in indices naming
#====================================================================
data(ple4)
data(ple4.index)
biofull <- 0.001*stock(ple4)
biofull <- FLIndexBiomass(index=biofull)
range(biofull)[c("startf","endf")] <- c(0,0)

fit0 <- sca(ple4, FLIndices(ple4.index, biofull), qmodel=list(~s(age, k=3), ~1))
fit1 <- sca(ple4, FLIndices(biofull, ple4.index), qmodel=list(~1, ~s(age, k=3)))

identical(stock.n(fit0), stock.n(fit1))
all.equal(c(harvest(fit0)), c(harvest(fit1)), tolerance = 1e-5)

#====================================================================
# center argument
#====================================================================
data(ple4)
data(ple4.indices)
fit0 <- sca(ple4, ple4.indices)
fit1 <- sca(ple4, ple4.indices, center=FALSE)
fit2 <- sca(ple4, ple4.indices, center=1)

# maxgrad is different
!identical(fitSumm(fit0)["maxgrad",], fitSumm(fit1)["maxgrad",])
!identical(fitSumm(fit0)["maxgrad",], fitSumm(fit2)["maxgrad",])
!identical(fitSumm(fit1)["maxgrad",], fitSumm(fit2)["maxgrad",])

# likelihood is equal
identical(fitSumm(fit0)["nlogl",], fitSumm(fit1)["nlogl",])
identical(fitSumm(fit0)["nlogl",], fitSumm(fit2)["nlogl",])
identical(fitSumm(fit1)["nlogl",], fitSumm(fit2)["nlogl",])

