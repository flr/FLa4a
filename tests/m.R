#====================================================================
# tests for sca
#====================================================================
library(FLa4a)

#====================================================================
# test a4aM age range
#====================================================================

mod02 <- FLModelSim(model=~a, params=FLPar(a=0.2))
m1 <- a4aM(level=mod02)
range(m1, c("min","max")) <- c(0,7) # set the quant range
range(m1, c("minyear","maxyear")) <- c(2000, 2010) # set the year range
dnms <- dimnames(m(m1))
identical(dnms$quant, ac(0:7))
identical(dnms$year, ac(2000:2010))
range(m1, c("minmbar","maxmbar")) <- c(1,5)
identical(unname(m1@range[c("minmbar","maxmbar")]), c(1,5)) 

# mbar will change the age range accordingly
range(m1, c("minmbar","maxmbar")) <- c(0,10)
dnms <- dimnames(m(m1))
identical(dnms$quant, ac(0:10))

# age range will change mbar accordingly
range(m1, c("min","max")) <- c(1,5)
identical(unname(m1@range[c("minmbar","maxmbar")]), c(1,5)) 

# if min > max take min
range(m1, c("min","max")) <- c(5,1)
dnms <- dimnames(m(m1))
identical(dnms$quant, ac(5))

# if max < min take min
range(m1, c("min","max")) <- c(1,5)
range(m1, c("max")) <- c(0)
dnms <- dimnames(m(m1))
identical(dnms$quant, ac(1))

# if minmbar > maxmbar take minmbar
range(m1, c("min","max")) <- c(1,7)
range(m1, c("minmbar","maxmbar")) <- c(2,4)
range(m1, c("minmbar")) <- c(5)
identical(unname(m1@range[c("minmbar","maxmbar")]), c(5,5)) 

# if maxmbar < minmbar take minmbar
range(m1, c("min","max")) <- c(1,7)
range(m1, c("minmbar","maxmbar")) <- c(2,4)
range(m1, c("maxmbar")) <- c(1)
identical(unname(m1@range[c("minmbar","maxmbar")]), c(2,2)) 

