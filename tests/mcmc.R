# mcmc.R - DESC
# /mcmc.R

# Copyright European Union, 2019
# Author: Iago Mosqueira (EC JRC) <iago.mosqueira@ec.europa.eu>
#
# Distributed under the terms of the European Union Public Licence (EUPL) V.1.1.

library(FLa4a)

#====================================================================
# MCMC class and fit
#====================================================================

data(ple4)
data(ple4.indices)

mcmcobj <- new("SCAMCMC")
identical(getADMBCallArgs(mcmcobj), " -mcmc 10000 -mcsave 100")

obj <- a4aFitMCMC()
is(obj, "a4aFitMCMC")
is(slot(obj, "mcmc"), "SCAMCMC")
obj <- a4aFitSA(obj)
is(obj, "a4aFitSA")
obj <- a4aFitMCMC(obj)
is(obj, "a4aFitMCMC")
is(slot(obj, "mcmc"), "SCAMCMC")
obj <- a4aFitMCMC(obj, mcmc=SCAMCMC())
is(obj, "a4aFitMCMC")
is(slot(obj, "mcmc"), "SCAMCMC")

mc <- SCAMCMC()
fit0 <- FLa4a:::a4aInternal(ple4, ple4.indices, fit="MCMC", mcmc=mc)
fit00 <- sca(ple4, ple4.indices, fit="MCMC", mcmc=mc)
identical(dimnames(catch.n(fit0))[-6], dimnames(catch.n(fit00))[-6])
identical(dimnames(stock.n(fit0))[-6], dimnames(stock.n(fit00))[-6])
dim(catch.n(fit0))[6]==getN(mc)
dim(catch.n(fit00))[6]==getN(mc)
sum(!is.na(fit0@pars@stkmodel@vcov))==0
sum(!is.na(fit00@pars@stkmodel@vcov))==0
sum(unlist(lapply(lapply(lapply(fit0@pars@vmodel, vcov), is.na), "!")))==0
sum(unlist(lapply(lapply(lapply(fit00@pars@vmodel, vcov), is.na), "!")))==0
sum(unlist(lapply(lapply(lapply(fit0@pars@qmodel, vcov), is.na), "!")))==0
sum(unlist(lapply(lapply(lapply(fit00@pars@qmodel, vcov), is.na), "!")))==0

mc <- SCAMCMC(mcmc=as.integer(10000), mcsave=as.integer(200), mcu=TRUE)
fit1 <- FLa4a:::a4aInternal(ple4, ple4.indices, fit="MCMC", mcmc=mc)
fit11 <- sca(ple4, ple4.indices, fit="MCMC", wkdir="mcmcalt2", mcmc=mc)
identical(dimnames(catch.n(fit1))[-6], dimnames(catch.n(fit11))[-6])
identical(dimnames(stock.n(fit1))[-6], dimnames(stock.n(fit11))[-6])
dim(catch.n(fit1))[6]==getN(mc)
dim(catch.n(fit11))[6]==getN(mc)
sum(!is.na(fit1@pars@stkmodel@vcov))==0
sum(!is.na(fit11@pars@stkmodel@vcov))==0
sum(unlist(lapply(lapply(lapply(fit1@pars@vmodel, vcov), is.na), "!")))==0
sum(unlist(lapply(lapply(lapply(fit11@pars@vmodel, vcov), is.na), "!")))==0
sum(unlist(lapply(lapply(lapply(fit1@pars@qmodel, vcov), is.na), "!")))==0
sum(unlist(lapply(lapply(lapply(fit11@pars@qmodel, vcov), is.na), "!")))==0

# must fail, with hybrid mcsave must be 1
is(try(SCAMCMC(hybrid=TRUE)), "try-error")
mc <- SCAMCMC(mcmc=100, mcsave=1, hybrid=TRUE)
fit2 <- FLa4a:::a4aInternal(ple4, ple4.indices, fit="MCMC", mcmc=mc)
fit22 <- sca(ple4, ple4.indices, fit="MCMC", mcmc=mc)
identical(dimnames(catch.n(fit2))[-6], dimnames(catch.n(fit22))[-6])
identical(dimnames(stock.n(fit2))[-6], dimnames(stock.n(fit22))[-6])
dim(catch.n(fit2))[6]==getN(mc)
dim(catch.n(fit22))[6]==getN(mc)
sum(!is.na(fit2@pars@stkmodel@vcov))==0
sum(!is.na(fit22@pars@stkmodel@vcov))==0
sum(unlist(lapply(lapply(lapply(fit2@pars@vmodel, vcov), is.na), "!")))==0
sum(unlist(lapply(lapply(lapply(fit22@pars@vmodel, vcov), is.na), "!")))==0
sum(unlist(lapply(lapply(lapply(fit2@pars@qmodel, vcov), is.na), "!")))==0
sum(unlist(lapply(lapply(lapply(fit22@pars@qmodel, vcov), is.na), "!")))==0

# check seed's working
mc <- SCAMCMC(mcmc=10, mcsave=1, hybrid=TRUE, mcseed=123)
fit3 <- sca(ple4, ple4.indices, fit="MCMC", mcmc=mc)
fit33 <- sca(ple4, ple4.indices, fit="MCMC", mcmc=mc)
identical(fit3@stock.n, fit33@stock.n)

#====================================================================
# iter
#====================================================================

data(ple4)
data(ple4.index)
stk <- propagate(ple4, 2)

# iter a4aFit
fit0 <- sca(ple4, FLIndices(a=ple4.index), fit='MP')
fit <- sca(stk, FLIndices(a=ple4.index), fit='MP')
fit1 <- iter(fit, 1)
identical(fit1@harvest, fit0@harvest)
identical(fit1@stock.n, fit0@stock.n)
identical(fit1@catch.n, fit0@catch.n)
identical(fit1@index, fit0@index)

# iter a4aFitSA
fit0 <- sca(ple4, FLIndices(a=ple4.index))
fit <- sca(stk, FLIndices(a=ple4.index))
fit1 <- iter(fit, 1)
identical(fit1@harvest, fit0@harvest)
identical(fit1@stock.n, fit0@stock.n)
identical(fit1@catch.n, fit0@catch.n)
identical(fit1@index, fit0@index)

