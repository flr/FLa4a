#====================================================================
# tests for multiple sca
#====================================================================
library(FLa4a)
data(ple4)
data(ple4.index)

#====================================================================
# run fits sca MP
#====================================================================

qmods <- list(list(~s(age, k=6)))
fmods = list()
for(i in 1:6) {
  fmods[[paste0(i)]] <- as.formula(paste0("~te(age, year, k = c(6,", i+14,"), bs = 'tp') + s(age, k = 6)"))
}

#--------------------------------------------------------------------
# MP 
#--------------------------------------------------------------------

fits <- FLa4a:::multisca(FLStocks(ple4), list(FLIndices(ple4.index)), fmodel = fmods, qmodel=qmods, fit="MP")

is(fits, "a4aFits")
is(fits[[1]], "a4aFit")


#--------------------------------------------------------------------
# SA 
#--------------------------------------------------------------------

fits <- FLa4a:::multisca(FLStocks(ple4), list(FLIndices(ple4.index)), fmodel = fmods, qmodel=qmods, fit="assessment")

is(fits, "a4aFitSAs")
is(fits[[1]], "a4aFitSA")

#--------------------------------------------------------------------
# MCMC 
#--------------------------------------------------------------------

fits <- FLa4a:::multisca(FLStocks(ple4), list(FLIndices(ple4.index)), fmodel = fmods, qmodel=qmods, fit="MCMC", mcmc=SCAMCMC())

is(fits, "a4aFitMCMCs")
is(fits[[1]], "a4aFitMCMC")


