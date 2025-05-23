#====================================================================
# tests for multiple sca
#====================================================================
library(FLa4a)
data(ple4)
data(ple4.index)
j <- 6

#====================================================================
# run fits sca 
#====================================================================

qmods <- list(list(~s(age, k=6)))
fmods = list()
for(i in 1:j) {
  fmods[[paste0(i)]] <- as.formula(paste0("~te(age, year, k = c(6,", i+14,"), bs = 'tp') + s(age, k = 6)"))
}
stks <- FLStocks(ple4)
idxss <- list(FLIndices(ple4.index)) 

#--------------------------------------------------------------------
# MP 
#--------------------------------------------------------------------

fits <- scas(stks, idxss, fmodel = fmods, qmodel=qmods, fit="MP")

is(fits, "a4aFits")
is(fits[[1]], "a4aFit")


#--------------------------------------------------------------------
# SA 
#--------------------------------------------------------------------

fitsa <- scas(stks, idxss, fmodel = fmods, qmodel=qmods, fit="assessment")

is(fitsa, "a4aFitSAs")
is(fitsa[[1]], "a4aFitSA")

#--------------------------------------------------------------------
# MCMC 
#--------------------------------------------------------------------

fitsm <- scas(stks, idxss, fmodel = fmods, qmodel=qmods, fit="MCMC", mcmc=SCAMCMC())

is(fitsm, "a4aFitMCMCs")
is(fitsm[[1]], "a4aFitMCMC")

#====================================================================
# +
#====================================================================

#--------------------------------------------------------------------
# MP 
#--------------------------------------------------------------------

stks.mp <- stks + fits
length(stks.mp) == length(fits)
identical(stks.mp[[2]], stks[[1]] + fits[[2]])

#--------------------------------------------------------------------
# SA 
#--------------------------------------------------------------------

stks.sa <- stks + fitsa
length(stks.sa) == length(fitsa)
identical(stks.sa[[2]], stks[[1]] + fitsa[[2]])

#--------------------------------------------------------------------
# MCMC 
#--------------------------------------------------------------------

stks.mc <- stks + fitsm
length(stks.mc) == length(fitsm)
identical(stks.mc[[2]], stks[[1]] + fitsm[[2]])

#====================================================================
# run fits sca multiple stocks and indices
#====================================================================

stks[1:j] <- stks[1]
idxss[1:j] <- idxss[1] 

#--------------------------------------------------------------------
# MP 
#--------------------------------------------------------------------

fits <- scas(stks, idxss, fmodel = fmods, qmodel=qmods, fit="MP")

is(fits, "a4aFits")
is(fits[[1]], "a4aFit")


#--------------------------------------------------------------------
# SA 
#--------------------------------------------------------------------

fitsa <- scas(stks, idxss, fmodel = fmods, qmodel=qmods, fit="assessment")

is(fitsa, "a4aFitSAs")
is(fitsa[[1]], "a4aFitSA")

#--------------------------------------------------------------------
# MCMC 
#--------------------------------------------------------------------

fitsm <- scas(stks, idxss, fmodel = fmods, qmodel=qmods, fit="MCMC", mcmc=SCAMCMC())

is(fitsm, "a4aFitMCMCs")
is(fitsm[[1]], "a4aFitMCMC")

#====================================================================
# +
#====================================================================

#--------------------------------------------------------------------
# MP 
#--------------------------------------------------------------------

stks.mp <- stks + fits
length(stks.mp) == length(fits)
identical(stks.mp[[2]], stks[[1]] + fits[[2]])

#--------------------------------------------------------------------
# SA 
#--------------------------------------------------------------------

stks.sa <- stks + fitsa
length(stks.sa) == length(fitsa)
identical(stks.sa[[2]], stks[[1]] + fitsa[[2]])

#--------------------------------------------------------------------
# MCMC 
#--------------------------------------------------------------------

stks.mc <- stks + fitsm
length(stks.mc) == length(fitsm)
identical(stks.mc[[2]], stks[[1]] + fitsm[[2]])





