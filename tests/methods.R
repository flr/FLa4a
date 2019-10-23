#====================================================================
# tests for accessors etc
#====================================================================

#library(FLa4a)
data(ple4)
data(ple4.index)

fmodel <- ~ s(age, k = 3)
qmodel <- list(~s(age, k = 3))
vmodel <- list(~1, ~1)
n1model <- ~ s(age, k = 3)
srmod <- ~ geomean(a = ~ year, CV = 0.5)

testfit <- sca(ple4, FLIndices(BTS = ple4.index),
           fmodel = fmodel,
           qmodel = qmodel,
           vmodel = vmodel,
           n1model = n1model,
           srmodel = srmod)

# accessors on the fit
testfit

pars(testfit)

stkmodel(testfit)

fmodel(testfit)
n1model(testfit)
rmodel(testfit)
srmodel(testfit)

qmodel(testfit)
vmodel(testfit)

qmodel(testfit)[1]
vmodel(testfit)[1:2]
vmodel(testfit)[-1]

qmodel(testfit)[[1]]
vmodel(testfit)[[1]]
vmodel(testfit)[[2]]

qmodel(testfit)[["BTS"]]
vmodel(testfit)[["catch"]]

#covariates(testfit)

m(testfit)
wt(testfit)
mat(testfit)

clock(testfit)
fitSumm(testfit)
summary(testfit)
stock.n(testfit)
harvest(testfit)
catch.n(testfit)
index(testfit)

logLik(testfit)

as(testfit, "FLSR")

