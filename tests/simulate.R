#====================================================================
# tests for simulate
#====================================================================
library(FLa4a)
data(ple4)
data(ple4.index)
data(ple4.indices)
nits <- 10
err <- 0.05
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

max(stk.rat)-min(stk.rat) < err
max(f.rat)-min(f.rat) < err
max(idx.rat)-min(idx.rat) < err

# is the vcov matrix ok
vce <- c(cov(t(params(qmodel(pars(obj))[[1]]))))
vc <- vcov(qmodel(pars(obj))[[1]])[,,1,drop=TRUE]
# check only those that have correlations stronger than 50%
i <- c(abs(cov2cor(vc))>0.5)
vrat <- vce[i]/c(vc)[i]
max(vrat)-min(vrat) < 2*err

vce <- c(cov(t(params(vmodel(pars(obj))[[1]]))))
vc <- vcov(vmodel(pars(obj))[[1]])[,,1,drop=TRUE]
# check only those that have correlations stronger than 50%
i <- c(abs(cov2cor(vc))>0.5)
vrat <- vce[i]/c(vc)[i]
max(vrat)-min(vrat) < 3*err

vce <- c(cov(t(params(vmodel(pars(obj))[[2]]))))
vc <- vcov(vmodel(pars(obj))[[2]])[,,1,drop=TRUE]
# check only those that have correlations stronger than 50%
i <- c(abs(cov2cor(vc))>0.5)
vrat <- vce[i]/c(vc)[i]
max(vrat)-min(vrat) < 2*err

vce <- c(cov(t(params(stkmodel(pars(obj))))))
vc <- vcov(stkmodel(pars(obj)))[,,1,drop=TRUE]
# check only those that have correlations stronger than 50%
i <- c(abs(cov2cor(vc))>0.5)
vrat <- vce[i]/c(vc)[i]
max(vrat)-min(vrat) < 4*err

# what about variance 
vce <- c(diag(cov(t(params(qmodel(pars(obj))[[1]])))))
vc <- diag(vcov(qmodel(pars(obj))[[1]])[,,1,drop=TRUE])
vrat <- vce/c(vc)
max(vrat)-min(vrat) < 2*err

vce <- c(diag(cov(t(params(vmodel(pars(obj))[[1]])))))
vc <- diag(vcov(vmodel(pars(obj))[[1]])[,,1,drop=TRUE])
vrat <- vce/c(vc)
max(vrat)-min(vrat) < 2*err

vce <- c(diag(cov(t(params(vmodel(pars(obj))[[2]])))))
vc <- diag(vcov(vmodel(pars(obj))[[2]])[,,1,drop=TRUE])
vrat <- vce/c(vc)
max(vrat)-min(vrat) < 2*err

vce <- c(diag(cov(t(params(stkmodel(pars(obj)))))))
vc <- diag(vcov(stkmodel(pars(obj)))[,,1,drop=TRUE])
vrat <- vce/c(vc)
max(vrat)-min(vrat) < 4*err

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

max(stk.rat)-min(stk.rat) < err
max(f.rat)-min(f.rat) < err
max(idx.rat)-min(idx.rat) < err
max(idx2.rat)-min(idx2.rat) < err
max(idx3.rat)-min(idx3.rat) < err

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

max(stk.rat)-min(stk.rat) < err
max(f.rat)-min(f.rat) < err
max(idx.rat)-min(idx.rat) < err

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

max(stk.rat)-min(stk.rat) < err
max(f.rat)-min(f.rat) < err
max(idx.rat)-min(idx.rat) < err
max(idx2.rat)-min(idx2.rat) < err





# is the vcov matrix ok
vce <- c(cov(t(params(qmodel(pars(obj))[[2]]))))
vc <- vcov(qmodel(pars(obj))[[2]])[,,1,drop=TRUE]
# check only those that have correlations stronger than 50%
i <- c(abs(cov2cor(vc))>0.5)
vrat <- vce[i]/c(vc)[i]
max(vrat)-min(vrat) < 3*err

vce <- c(cov(t(params(vmodel(pars(obj))[[1]]))))
vc <- vcov(vmodel(pars(obj))[[1]])[,,1,drop=TRUE]
# check only those that have correlations stronger than 50%
i <- c(abs(cov2cor(vc))>0.5)
vrat <- vce[i]/c(vc)[i]
max(vrat)-min(vrat) < 4*err

vce <- c(cov(t(params(stkmodel(pars(obj))))))
vc <- vcov(stkmodel(pars(obj)))[,,1,drop=TRUE]
# check only those that have correlations stronger than 50%
i <- c(abs(cov2cor(vc))>0.5)
vrat <- vce[i]/c(vc)[i]
max(vrat)-min(vrat) < 8*err

# what about variance 
vce <- c(diag(cov(t(params(qmodel(pars(obj))[[2]])))))
vc <- diag(vcov(qmodel(pars(obj))[[2]])[,,1,drop=TRUE])
vrat <- vce/c(vc)
max(vrat)-min(vrat) < 3*err

vce <- c(diag(cov(t(params(vmodel(pars(obj))[[1]])))))
vc <- diag(vcov(vmodel(pars(obj))[[1]])[,,1,drop=TRUE])
vrat <- vce/c(vc)
max(vrat)-min(vrat) < 4*err

vce <- c(diag(cov(t(params(stkmodel(pars(obj)))))))
vc <- diag(vcov(stkmodel(pars(obj)))[,,1,drop=TRUE])
vrat <- vce/c(vc)
max(vrat)-min(vrat) < 6*err

#====================================================================
# more than one iter in vcov, coming from assessments with iters
#====================================================================

fit0 <-  a4aSCA(ple4, FLIndices(ple4.index), qmodel=list(~s(age, k=4)))

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

max(stk.rat)-min(stk.rat) < 3*err
max(f.rat)-min(f.rat) < 3*err
max(idx.rat)-min(idx.rat) < 3*err


