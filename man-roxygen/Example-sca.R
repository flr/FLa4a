#' @examples
#' data(ple4)
#' data(ple4.index)
#' fit1 <- sca(ple4, FLIndices(ple4.index))
#' plot(ple4 + fit1)
#' fit2 <- sca(ple4, FLIndices(ple4.index), fmodel=~factor(age) + factor(year))
#' plot(ple4 + fit2)
#' fit3 <- sca(ple4, FLIndices(ple4.index), fmodel=~factor(age) + factor(year), qmodel=list(~s(age, k=4)))
#' plot(ple4 + fit3)
#' fit4 <- sca(ple4, FLIndices(ple4.index), fmodel=~factor(age) + factor(year), qmodel=list(~s(age, k=4)), srmodel=~s(year, k=45))
#' plot(ple4 + fit4)
#' AIC(fit1, fit2, fit3, fit4)
#' BIC(fit1, fit2, fit3, fit4)
#' 

