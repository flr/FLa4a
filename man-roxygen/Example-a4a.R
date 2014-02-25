#' @examples
#' data(ple4)
#' data(ple4.index)
#' 
#' #====================================================================
#' # fishing mortality by age and year (~separable)
#' # catchability at age without year trend
#' #====================================================================
#' 
#' fmodel <- ~factor(age) + factor(year)
#' qmodel <- list(~factor(age))
#' fit1 <- a4a(fmodel, qmodel, stock=ple4, indices=FLIndices(ple4.index))
#' 
#' #====================================================================
#' # fishing mortality as a smoother by age and year
#' # catchability at age without year trend
#' #====================================================================
#' 
#' fmodel <- ~ s(age, k=4) + s(year, k=10)
#' fit2 <- a4a(fmodel, qmodel, stock=ple4, indices=FLIndices(ple4.index))
#' 
#' #====================================================================
#' # fishing mortality as a smoother by age and year
#' # catchability at age without year trend
#' #====================================================================
#' 
#' fmodel <- ~ s(age, k=4) + s(year, k=10)
#' qmodel <- list(~s(age, k=4) + year)
#' fit3 <- a4a(fmodel, qmodel, stock=ple4, indices=FLIndices(ple4.index))
#' 
#' #====================================================================
#' # fishing mortality as a smoother by age and year
#' # catchability at age with year trend
#' #====================================================================
#' 
#' fmodel <- ~ te(age, year, k=c(4, 10))
#' qmodel <- list(~s(age, k=4))
#' fit4 <- a4a(fmodel, qmodel, stock=ple4, indices=FLIndices(ple4.index))
#' 
#' #====================================================================
#' # It's R/FLR !! Comparisons are easy :)
#' #====================================================================
#' 
#' fs <- FLQuants("log f ~ factor(age) + factor(year)"=fit1@harvest, "log f ~ s(age, k=4) + s(year, k=10)"=fit2@harvest, "log q ~ s(age, k=4) + year"=fit3@harvest, "log f ~ te(age, year, k=c(4, 10))"=fit4@harvest)
#' #wireframe(data~age*year|qname,data=as.data.frame(fs), screen = list(x = -90, y=-45), drape=T, layout=c(2,2), as.table=T, zlab="")
#' 
#' #====================================================================
#' # It's a statistical model 
#' #====================================================================
#' 
#' BIC(fit1, fit2, fit3)
#' 
#' #====================================================================
#' # f3 + smoother in recruitment 
#' #====================================================================
#' fmodel <- ~ s(age, k=4) + s(year, k=20)
#' qmodel <- list(~s(age, k=4))
#' rmodel <- ~s(year, k=20)
#' fit5 <- a4a(fmodel, qmodel, rmodel, ple4, FLIndices(ple4.index))
#' #xyplot(data~year, groups=qname, data=FLQuants(smooth=fit5@stock.n[1], nosmooth=fit3@stock.n[1]), type="l")
#' 
#' #====================================================================
#' # f3 + bevholt 
#' #====================================================================
#' rmodel <- ~ bevholt(CV=0.05)
#' fit6 <- a4a(fmodel, qmodel, rmodel, ple4, FLIndices(ple4.index))
#' #xyplot(data~year, groups=qname, data=FLQuants(smooth=fit6@stock.n[1], nosmooth=fit3@stock.n[1]), type="l")

