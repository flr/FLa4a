#====================================================================
# tests for simulate
#====================================================================
library(FLa4a)
data(ple4)
data(ple4.index)
data(ple4.indices)
nits <- 1000
err <- 0.05

#====================================================================
# check stkmodel
#====================================================================
fit <-  a4aSCA(ple4, FLIndices(ple4.index), qmodel=list(~s(age, k=4)), vmodel=list(~s(age, k=4), ~s(age, k=4)))

# replicate with set.seed
stkm <- stkmodel(pars(fit))
obj <- simulate(stkm, seed=1234)
identical(dim(params(obj)), dim(params(stkm)))
identical(dim(vcov(obj)), dim(vcov(stkm)))

# empirical is different
obj1 <- simulate(stkm, empirical=FALSE, seed=1234)
!identical(obj1, obj)

# but seed sould replicate
obj1 <- simulate(stkm, seed=1234)
identical(obj1, obj)

# vcov should be the same ...
obj <- simulate(stkm, nits, seed=1234)
dim(params(obj))["iter"] == nits
identical(dim(vcov(obj)), dim(vcov(stkm)))
rat <- cov(t(params(obj)))/vcov(stkm)[,,1,drop=T]
max(rat) - min(rat) < err/10

# while here it shouldn't ...
obj <- simulate(stkm, nits, empirical = FALSE)
rat <- cov(t(params(obj)))/vcov(stkm)[,,1,drop=T]
max(rat) - min(rat) > err*4

#====================================================================
# abundance indices
#====================================================================

#--------------------------------------------------------------------
# 1
#--------------------------------------------------------------------
fit <-  a4aSCA(ple4, FLIndices(ple4.index), qmodel=list(~s(age, k=4)), vmodel=list(~s(age, k=4), ~s(age, k=4)))

# check it runs
obj <- simulate(fit)
is(obj, "a4aFitSA")
validObject(obj)

# check vcov is not full of zeros
sum(vcov(stkmodel(pars(obj)))==0)!=length(vcov(stkmodel(pars(obj))))
sum(vcov(qmodel(pars(obj))[[1]])==0)!=length(vcov(qmodel(pars(obj))[[1]]))
sum(vcov(vmodel(pars(obj))[[1]])==0)!=length(vcov(vmodel(pars(obj))[[1]]))
sum(vcov(vmodel(pars(obj))[[2]])==0)!=length(vcov(vmodel(pars(obj))[[2]]))

# simulate from simulated object
obj <- simulate(obj)
is(obj, "a4aFitSA")
validObject(obj)

# check it runs with nits
obj <- simulate(fit, nits)
is(obj, "a4aFitSA")
validObject(obj)
dim(catch.n(obj))[6] == nits

# check vcov is not full of zeros
sum(vcov(stkmodel(pars(obj)))==0)!=length(vcov(stkmodel(pars(obj))))
sum(vcov(qmodel(pars(obj))[[1]])==0)!=length(vcov(qmodel(pars(obj))[[1]]))
sum(vcov(vmodel(pars(obj))[[1]])==0)!=length(vcov(vmodel(pars(obj))[[1]]))
sum(vcov(vmodel(pars(obj))[[2]])==0)!=length(vcov(vmodel(pars(obj))[[2]]))

# check vcov has only one iter
dim(vcov(stkmodel(pars(obj))))[3]==1
dim(vcov(qmodel(pars(obj))[[1]]))[3]==1
dim(vcov(vmodel(pars(obj))[[1]]))[3]==1
dim(vcov(vmodel(pars(obj))[[2]]))[3]==1

# simulate from simulated object
obj <- simulate(obj, nits)
is(obj, "a4aFitSA")
validObject(obj)
dim(catch.n(obj))[6] == nits

# can the seed be controled ?
obj0 <- simulate(fit, nits, 1234)
obj1 <- simulate(fit, nits, 1234)
all.equal(obj0, obj1)

# are the simulated values unbiased
obj <- simulate(fit, nits, seed=1234)
stk.fit <- stock.n(fit)
stk.sim <- stock.n(obj)
stk.rat <- iterMedians(stk.sim)/stk.fit
f.fit <- harvest(fit)
f.sim <- harvest(obj)
f.rat <- iterMedians(f.sim)/f.fit
idx.fit <- index(fit)[[1]]
idx.sim <- index(obj)[[1]]
idx.rat <- iterMedians(idx.sim)/idx.fit

#max(stk.rat)-min(stk.rat) < err
#max(f.rat)-min(f.rat) < err
#max(idx.rat)-min(idx.rat) < err

mean(stk.rat) > 1 - err & mean(stk.rat) < 1 + err
mean(f.rat) > 1 - err & mean(f.rat) < 1 + err
mean(idx.rat) > 1 - err & mean(idx.rat) < 1 + err

# is the vcov matrix ok
vce <- c(cov(t(params(qmodel(pars(obj))[[1]]))))
vc <- vcov(qmodel(pars(obj))[[1]])[,,1,drop=TRUE]
vrat <- vce/c(vc)
max(vrat)-min(vrat) < err/10

vce <- c(cov(t(params(vmodel(pars(obj))[[1]]))))
vc <- vcov(vmodel(pars(obj))[[1]])[,,1,drop=TRUE]
vrat <- vce/c(vc)
max(vrat)-min(vrat) < err/10

vce <- c(cov(t(params(vmodel(pars(obj))[[2]]))))
vc <- vcov(vmodel(pars(obj))[[2]])[,,1,drop=TRUE]
vrat <- vce/c(vc)
max(vrat)-min(vrat) < err/10

vce <- c(cov(t(params(stkmodel(pars(obj))))))
vc <- vcov(stkmodel(pars(obj)))[,,1,drop=TRUE]
vrat <- vce/c(vc)
max(vrat)-min(vrat) < err/10

#--------------------------------------------------------------------
# several
#--------------------------------------------------------------------
fit <-  a4aSCA(ple4, ple4.indices, qmodel=list(~s(age, k=4), ~s(age, k=4), ~s(age, k=3)), vmodel=list(~s(age, k=4), ~s(age, k=4), ~s(age, k=4), ~s(age, k=3)))

# check
obj <- simulate(fit)
is(obj, "a4aFitSA")
validObject(obj)

obj <- simulate(fit, nits)
is(obj, "a4aFitSA")
validObject(obj)
dim(catch.n(obj))[6] == nits

# can the seed be controled ?
obj0 <- simulate(fit, nits, 1234)
obj1 <- simulate(fit, nits, 1234)
all.equal(obj0, obj1)

# are the simulated values unbiased
obj <- simulate(fit, 1000)
stk.fit <- stock.n(fit)
stk.sim <- stock.n(obj)
stk.rat <- iterMedians(stk.sim)/stk.fit
f.fit <- harvest(fit)
f.sim <- harvest(obj)
f.rat <- iterMedians(f.sim)/f.fit
idx.fit <- index(fit)[[1]]
idx.sim <- index(obj)[[1]]
idx.rat <- iterMedians(idx.sim)/idx.fit
idx2.fit <- index(fit)[[2]]
idx2.sim <- index(obj)[[2]]
idx2.rat <- iterMedians(idx2.sim)/idx2.fit
idx3.fit <- index(fit)[[3]]
idx3.sim <- index(obj)[[3]]
idx3.rat <- iterMedians(idx3.sim)/idx3.fit

mean(stk.rat) > 1 - err & mean(stk.rat) < 1 + err
mean(f.rat) > 1 - err & mean(f.rat) < 1 + err
mean(idx.rat) > 1 - err & mean(idx.rat) < 1 + err
mean(idx2.rat) > 1 - err & mean(idx2.rat) < 1 + err
mean(idx3.rat) > 1 - err & mean(idx3.rat) < 1 + err

#====================================================================
# biomass index
#====================================================================
bioidx <- FLIndexBiomass(FLQuant(NA, dimnames=list(age="all", year=range(ple4)["minyear"]:range(ple4)["maxyear"])), name="bioidx")
index(bioidx) <- stock(ple4)*0.001
index(bioidx) <- index(bioidx)*exp(rnorm(index(bioidx), sd=0.1))
range(bioidx)[c("startf","endf")] <- c(0,0)

# fitting the model
fit <- a4aSCA(ple4, FLIndices(bioidx), qmodel=list(~1))

#--------------------------------------------------------------------
# check
#--------------------------------------------------------------------

# now the classes
obj <- simulate(fit)
is(obj, "a4aFitSA")
validObject(obj)

obj <- simulate(fit, nits)
is(obj, "a4aFitSA")
validObject(obj)
dim(catch.n(obj))[6] == nits

# can the seed be controled ?
set.seed(1234)
obj0 <- simulate(fit, nits)
set.seed(1234)
obj1 <- simulate(fit, nits)
all.equal(obj0, obj1)

# are the simulated values unbiased
obj <- simulate(fit, 1000)
stk.fit <- stock.n(fit)
stk.sim <- stock.n(obj)
stk.rat <- iterMedians(stk.sim)/stk.fit
f.fit <- harvest(fit)
f.sim <- harvest(obj)
f.rat <- iterMedians(f.sim)/f.fit
idx.fit <- index(fit)[[1]]
idx.sim <- index(obj)[[1]]
idx.rat <- iterMedians(idx.sim)/idx.fit

mean(stk.rat) > 1 - err & mean(stk.rat) < 1 + err
mean(f.rat) > 1 - err & mean(f.rat) < 1 + err
mean(idx.rat) > 1 - err & mean(idx.rat) < 1 + err

#====================================================================
# biomass and abundance indices
#====================================================================

# fitting the model
fit <- sca(ple4, FLIndices(bioidx, ple4.index), qmodel=list(~1, ~s(age, k=4)), fit="assessment")

#--------------------------------------------------------------------
# check
#--------------------------------------------------------------------

# now the classes
obj <- simulate(fit)
is(obj, "a4aFitSA")
validObject(obj)

obj <- simulate(fit, nits)
is(obj, "a4aFitSA")
validObject(obj)
dim(catch.n(obj))[6] == nits

# can the seed be controled ?
obj0 <- simulate(fit, nits, 1234)
obj1 <- simulate(fit, nits, 1234)
all.equal(obj0, obj1)

# are the simulated values unbiased
obj <- simulate(fit, 1000)
stk.fit <- stock.n(fit)
stk.sim <- stock.n(obj)
stk.rat <- iterMedians(stk.sim)/stk.fit
f.fit <- harvest(fit)
f.sim <- harvest(obj)
f.rat <- iterMedians(f.sim)/f.fit
idx.fit <- index(fit)[[1]]
idx.sim <- index(obj)[[1]]
idx.rat <- iterMedians(idx.sim)/idx.fit
idx2.fit <- index(fit)[[2]]
idx2.sim <- index(obj)[[2]]
idx2.rat <- iterMedians(idx2.sim)/idx2.fit

mean(stk.rat) > 1 - err & mean(stk.rat) < 1 + err
mean(f.rat) > 1 - err & mean(f.rat) < 1 + err
mean(idx.rat) > 1 - err & mean(idx.rat) < 1 + err
mean(idx2.rat) > 1 - err & mean(idx2.rat) < 1 + err

# is the vcov matrix ok
vce <- c(cov(t(params(qmodel(pars(obj))[[2]]))))
vc <- vcov(qmodel(pars(obj))[[2]])[,,1,drop=TRUE]
vrat <- vce/c(vc)
max(vrat)-min(vrat) < err

vce <- c(cov(t(params(vmodel(pars(obj))[[1]]))))
vc <- vcov(vmodel(pars(obj))[[1]])[,,1,drop=TRUE]
vrat <- vce/c(vc)
max(vrat)-min(vrat) < err

vce <- c(cov(t(params(stkmodel(pars(obj))))))
vc <- vcov(stkmodel(pars(obj)))[,,1,drop=TRUE]
vrat <- vce/c(vc)
max(vrat)-min(vrat) < err

#====================================================================
# more than one iter in vcov, coming from assessments with iters
#====================================================================

fit0 <-  a4aSCA(ple4, FLIndices(ple4.index), qmodel=list(~s(age, k=4)))
nits <- 3

stk2 <- ple4
idx2 <- ple4.index
catch.n(stk2) <- genFLQuant(catch.n(stk2), 0.1, niter=nits)
index(idx2) <- genFLQuant(index(fit0)[[1]], 0.1, niter=nits)

fit <- a4aSCA(stk2, FLIndices(idx2), qmodel=list(~s(age, k=4)))

# check it runs
obj <- simulate(fit)
is(obj, "a4aFitSA")
validObject(obj)

# check it fails with nsim != nits
is(try(simulate(fit, nits+1)), "try-error")

# check vcov has nits iters
dim(vcov(stkmodel(pars(obj))))[3]==nits
dim(vcov(qmodel(pars(obj))[[1]]))[3]==nits
dim(vcov(vmodel(pars(obj))[[1]]))[3]==nits
dim(vcov(vmodel(pars(obj))[[2]]))[3]==nits

# check vcov is not full of zeros
sum(vcov(stkmodel(pars(obj)))==0)!=length(vcov(stkmodel(pars(obj))))
sum(vcov(qmodel(pars(obj))[[1]])==0)!=length(vcov(qmodel(pars(obj))[[1]]))
sum(vcov(vmodel(pars(obj))[[1]])==0)!=length(vcov(vmodel(pars(obj))[[1]]))
sum(vcov(vmodel(pars(obj))[[2]])==0)!=length(vcov(vmodel(pars(obj))[[2]]))

# simulate from simulated object
obj <- simulate(obj)
is(obj, "a4aFitSA")
validObject(obj)

# check vcov has nits iters
dim(vcov(stkmodel(pars(obj))))[3]==nits
dim(vcov(qmodel(pars(obj))[[1]]))[3]==nits
dim(vcov(vmodel(pars(obj))[[1]]))[3]==nits
dim(vcov(vmodel(pars(obj))[[2]]))[3]==nits

# check it runs with nits
obj <- simulate(fit, nits)
is(obj, "a4aFitSA")
validObject(obj)
dim(catch.n(obj))[6] == nits

# check vcov has nits iters
dim(vcov(stkmodel(pars(obj))))[3]==nits
dim(vcov(qmodel(pars(obj))[[1]]))[3]==nits
dim(vcov(vmodel(pars(obj))[[1]]))[3]==nits
dim(vcov(vmodel(pars(obj))[[2]]))[3]==nits

# check vcov is not full of zeros
sum(vcov(stkmodel(pars(obj)))==0)!=length(vcov(stkmodel(pars(obj))))
sum(vcov(qmodel(pars(obj))[[1]])==0)!=length(vcov(qmodel(pars(obj))[[1]]))
sum(vcov(vmodel(pars(obj))[[1]])==0)!=length(vcov(vmodel(pars(obj))[[1]]))
sum(vcov(vmodel(pars(obj))[[2]])==0)!=length(vcov(vmodel(pars(obj))[[2]]))

# simulate from simulated object
obj <- simulate(obj, nits)
is(obj, "a4aFitSA")
validObject(obj)
dim(catch.n(obj))[6] == nits

# can the seed be controled ?
obj0 <- simulate(fit, nits, 1234)
obj1 <- simulate(fit, nits, 1234)
all.equal(obj0, obj1)

# are the simulated values unbiased
obj <- simulate(fit, nits)
stk.fit <- stock.n(fit)
stk.sim <- stock.n(obj)
stk.rat <- iterMedians(stk.sim)/iterMedians(stk.fit)
f.fit <- harvest(fit)
f.sim <- harvest(obj)
f.rat <- iterMedians(f.sim)/iterMedians(f.fit)
idx.fit <- index(fit)[[1]]
idx.sim <- index(obj)[[1]]
idx.rat <- iterMedians(idx.sim)/iterMedians(idx.fit)

mean(stk.rat) > 1 - err & mean(stk.rat) < 1 + err
mean(f.rat) > 1 - err & mean(f.rat) < 1 + err
mean(idx.rat) > 1 - err & mean(idx.rat) < 1 + err

