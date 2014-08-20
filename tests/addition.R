#====================================================================
# tests for predict
#====================================================================
library(FLa4a)
data(ple4)
data(ple4.indices)
data(ple4.index)
nits <- 2

#====================================================================
# abundance indices
#====================================================================
#--------------------------------------------------------------------
# both assessment types
#--------------------------------------------------------------------
# 1 x 1
fit0 <- sca(ple4, FLIndices(ple4.index), qmodel=list(~s(age, k=4)))
stk0 <- ple4 + fit0
fit <- sca(ple4, FLIndices(ple4.index), qmodel=list(~s(age, k=4)), fit="assessment")
stk <- ple4 + fit
all.equal(stk, stk0)

#--------------------------------------------------------------------
# N
#--------------------------------------------------------------------
# N x 1
set.seed(1234)
stk0 <- propagate(ple4, nits) + fit
dims(stk0)$iter==nits

# are the right slots being copyed
identical(c(stock.n(stk0)), c(stock.n(fit)[,,,,,rep(1,nits)]))
identical(c(catch.n(stk0)), c(catch.n(fit)[,,,,,rep(1,nits)]))
identical(c(harvest(stk0)), c(harvest(fit)[,,,,,rep(1,nits)]))

set.seed(1234)
stk <- propagate(ple4, nits) * fit
dims(stk)$iter==nits

# ref sim for comparison
fits <- simulate(fit, nits, 1234)

# are the right slots being copyed
identical(stock.n(stk), stock.n(fits))
identical(catch.n(stk), catch.n(fits))
identical(harvest(stk), harvest(fits))

# must be different because there's no simulation involved in stk0
!identical(stk, stk0)

# 1 x N
stk0 <- ple4 + fits
dims(stk0)$iter==nits

# are the right slots being copyed
identical(stock.n(stk0), stock.n(fits))
identical(catch.n(stk0), catch.n(fits))
identical(harvest(stk0), harvest(fits))

set.seed(1234)
stk <- ple4 * propagate(fit, nits)
dims(stk)$iter==nits

# are the right slots being copyed
identical(stock.n(stk), stock.n(fits))
identical(catch.n(stk), catch.n(fits))
identical(harvest(stk), harvest(fits))

# must be equal because simulations are controled in both cases
identical(stk, stk0)

# N x N
stk0 <- propagate(ple4, 2) + fits
dims(stk0)$iter==nits

# are the right slots being copyed
identical(stock.n(stk0), stock.n(fits))
identical(catch.n(stk0), catch.n(fits))
identical(harvest(stk0), harvest(fits))

set.seed(1234)
stk <- propagate(ple4, 2) * propagate(fit, 2)
dims(stk)$iter==nits

# are the right slots being copyed
identical(stock.n(stk), stock.n(fits))
identical(catch.n(stk), catch.n(fits))
identical(harvest(stk), harvest(fits))

# must be equal because simulations are controled in both cases
identical(stk, stk0)

#====================================================================
# biomass indices
#====================================================================

# creating idx 1
bioidx <- FLIndexBiomass(FLQuant(NA, dimnames=list(age="all", year=range(ple4)["minyear"]:range(ple4)["maxyear"])), name="bioidx")
index(bioidx) <- stock(ple4)*0.001
index(bioidx) <- index(bioidx)*exp(rnorm(index(bioidx), sd=0.1))
range(bioidx)[c("startf","endf")] <- c(0,0)

#--------------------------------------------------------------------
# both assessment types
#--------------------------------------------------------------------
fit0 <- sca(ple4, FLIndices(bioidx), qmodel=list(~1))
stk0 <- ple4 + fit0
fit <- sca(ple4, FLIndices(bioidx), qmodel=list(~1), fit="assessment")
stk <- ple4 + fit
all.equal(stk, stk0)

#--------------------------------------------------------------------
# N
#--------------------------------------------------------------------
# N x 1
set.seed(1234)
stk0 <- propagate(ple4, nits) + fit
dims(stk0)$iter==nits

# are the right slots being copyed
identical(c(stock.n(stk0)), c(stock.n(fit)[,,,,,rep(1,nits)]))
identical(c(catch.n(stk0)), c(catch.n(fit)[,,,,,rep(1,nits)]))
identical(c(harvest(stk0)), c(harvest(fit)[,,,,,rep(1,nits)]))

# ref sim for comparison
fits <- simulate(fit, nits, 1234, ple4)

# 1 x N
stk0 <- ple4 + fits
dims(stk0)$iter==nits

# are the right slots being copyed
identical(stock.n(stk0), stock.n(fits))
identical(catch.n(stk0), catch.n(fits))
identical(harvest(stk0), harvest(fits))

# N x N
stk0 <- propagate(ple4, 2) + fits
dims(stk0)$iter==nits

# are the right slots being copyed
identical(stock.n(stk0), stock.n(fits))
identical(catch.n(stk0), catch.n(fits))
identical(harvest(stk0), harvest(fits))

## must fail
##stk0 <- propagate(ple4, 2) + fit0

## N x 1
#set.seed(1234)
#stk0 <- propagate(ple4, 2) + fit
#set.seed(1234)
#stk <- propagate(ple4, 2) * fit
#all.equal(stk, stk0)

## 1 x N
#stk0 <- ple4 + simulate(fit, 2, 1234)
#stk <- ple4 * simulate(fit, 2, 1234)
#all.equal(stk, stk0)

#stk0 <- ple4 + simulate(fit, 2, 1234, ple4)
#stk <- ple4 * simulate(fit, 2, 1234, ple4)
#all.equal(stk, stk0)

## N x N
#stk0 <- propagate(ple4, 2) + simulate(fit, 2, 1234)
#stk <- propagate(ple4, 2) * simulate(fit, 2, 1234)
#all.equal(stk, stk0)

#====================================================================
# both indices
#====================================================================

#--------------------------------------------------------------------
# both assessment types
#--------------------------------------------------------------------
fit0 <- sca(ple4, FLIndices(bioidx, ple4.index), qmodel=list(~1, ~s(age, k=4)))
stk0 <- ple4 + fit0
fit <- sca(ple4, FLIndices(bioidx, ple4.index), qmodel=list(~1, ~s(age, k=4)), fit="assessment")
stk <- ple4 + fit
all.equal(stk, stk0)

#--------------------------------------------------------------------
# N
#--------------------------------------------------------------------
# N x 1
set.seed(1234)
stk0 <- propagate(ple4, nits) + fit
dims(stk0)$iter==nits

# are the right slots being copyed
identical(c(stock.n(stk0)), c(stock.n(fit)[,,,,,rep(1,nits)]))
identical(c(catch.n(stk0)), c(catch.n(fit)[,,,,,rep(1,nits)]))
identical(c(harvest(stk0)), c(harvest(fit)[,,,,,rep(1,nits)]))

# ref sim for comparison
fits <- simulate(fit, nits, 1234, ple4)

# 1 x N
stk0 <- ple4 + fits
dims(stk0)$iter==nits

# are the right slots being copyed
identical(stock.n(stk0), stock.n(fits))
identical(catch.n(stk0), catch.n(fits))
identical(harvest(stk0), harvest(fits))

# N x N
stk0 <- propagate(ple4, 2) + fits
dims(stk0)$iter==nits

# are the right slots being copyed
identical(stock.n(stk0), stock.n(fits))
identical(catch.n(stk0), catch.n(fits))
identical(harvest(stk0), harvest(fits))

## must fail
##stk0 <- propagate(ple4, 2) + fit0

## N x 1
#set.seed(1234)
#stk0 <- propagate(ple4, 2) + fit
#set.seed(1234)
#stk <- propagate(ple4, 2) * fit
#all.equal(stk, stk0)

## 1 x N
#stk0 <- ple4 + simulate(fit, 2, 1234)
#stk <- ple4 * simulate(fit, 2, 1234)
#all.equal(stk, stk0)

#stk0 <- ple4 + simulate(fit, 2, 1234, ple4)
#stk <- ple4 * simulate(fit, 2, 1234, ple4)
#all.equal(stk, stk0)

## N x N
#stk0 <- propagate(ple4, 2) + simulate(fit, 2, 1234)
#stk <- propagate(ple4, 2) * simulate(fit, 2, 1234)
#all.equal(stk, stk0)


