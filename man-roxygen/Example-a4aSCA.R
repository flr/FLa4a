#' @examples
#' data(ple4)
#' data(ple4.index)
#'
#' # fishing mortality by age and year (separable) AND catchability at age without year trend
#' fmodel <- ~factor(age) + factor(year)
#' qmodel <- list(~factor(age))
#' fit1 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#' 
#' # fishing mortality as a smoother by age and year (but still separable) AND catchability at age without year trend
#' fmodel <- ~ s(age, k=4) + s(year, k=10)
#' qmodel <- list(~factor(age))
#' fit2 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#' 
#' # fishing mortality as a smoother by age and year (but still separable) AND catchability as a smoother by age without year trend
#' fmodel <- ~ s(age, k=4) + s(year, k=10)
#' qmodel <- list(~s(age, k=4))
#' fit3 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#'
#' # fishing mortality as a smoother by age and year (but still separable) AND catchability as a smoother by age with year trend
#' fmodel <- ~ s(age, k=4) + s(year, k=10)
#' qmodel <- list(~s(age, k=4) + year)
#' fit4 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#' 
#' # It's a statistical model
#' BIC(fit1, fit2, fit3, fit4)
#'
#' # fishing mortality as a smoother by age and year with interactions (i.e. non-separable) AND catchability as a smoother by age without year trend
#' fmodel <- ~ te(age, year, k=c(4, 10))
#' qmodel <- list(~s(age, k=4))
#' fit5 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#' 
#' # fit3 + smoother in recruitment
#' fmodel <- ~ s(age, k=4) + s(year, k=20)
#' qmodel <- list(~s(age, k=4))
#' rmodel <- ~s(year, k=20)
#' fit6 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, srmodel=rmodel, ple4, FLIndices(ple4.index))
#' 
#' # fit3 + bevholt 
#' rmodel <- ~ bevholt(CV=0.05)
#' fit7 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, srmodel=rmodel, ple4, FLIndices(ple4.index))



