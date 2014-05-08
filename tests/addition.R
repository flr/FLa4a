#====================================================================
# tests for predict
#====================================================================
library(FLa4a)
data(ple4)
data(ple4.indices)
data(ple4.index)

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
# must fail
#stk0 <- propagate(ple4, 2) + fit0

# N x 1
set.seed(1234)
stk0 <- propagate(ple4, 2) + fit
set.seed(1234)
stk <- propagate(ple4, 2) * fit
all.equal(stk, stk0)

# 1 x N
set.seed(1234)
stk0 <- ple4 + simulate(fit, 2)
set.seed(1234)
stk <- ple4 * simulate(fit, 2)
all.equal(stk, stk0)

# N x N
set.seed(1234)
stk0 <- propagate(ple4, 2) + simulate(fit, 2)
set.seed(1234)
stk <- propagate(ple4, 2) * simulate(fit, 2)
all.equal(stk, stk0)

#====================================================================
# biomass indices
#====================================================================

# creating idx 1
bioidx <- FLIndex(FLQuant(NA, dimnames=list(age="all", year=range(ple4)["minyear"]:range(ple4)["maxyear"])), name="bioidx")
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
# must fail
#stk0 <- propagate(ple4, 2) + fit0

# N x 1
set.seed(1234)
stk0 <- propagate(ple4, 2) + fit
set.seed(1234)
stk <- propagate(ple4, 2) * fit
all.equal(stk, stk0)

# 1 x N
set.seed(1234)
stk0 <- ple4 + simulate(fit, 2)
set.seed(1234)
stk <- ple4 * simulate(fit, 2)
all.equal(stk, stk0)

set.seed(1234)
stk0 <- ple4 + simulate(fit, 2, ple4)
set.seed(1234)
stk <- ple4 * simulate(fit, 2, ple4)
all.equal(stk, stk0)

# N x N
set.seed(1234)
stk0 <- propagate(ple4, 2) + simulate(fit, 2)
set.seed(1234)
stk <- propagate(ple4, 2) * simulate(fit, 2)
all.equal(stk, stk0)

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
# must fail
#stk0 <- propagate(ple4, 2) + fit0

# N x 1
set.seed(1234)
stk0 <- propagate(ple4, 2) + fit
set.seed(1234)
stk <- propagate(ple4, 2) * fit
all.equal(stk, stk0)

# 1 x N
set.seed(1234)
stk0 <- ple4 + simulate(fit, 2)
set.seed(1234)
stk <- ple4 * simulate(fit, 2)
all.equal(stk, stk0)

set.seed(1234)
stk0 <- ple4 + simulate(fit, 2, ple4)
set.seed(1234)
stk <- ple4 * simulate(fit, 2, ple4)
all.equal(stk, stk0)

# N x N
set.seed(1234)
stk0 <- propagate(ple4, 2) + simulate(fit, 2)
set.seed(1234)
stk <- propagate(ple4, 2) * simulate(fit, 2)
all.equal(stk, stk0)


