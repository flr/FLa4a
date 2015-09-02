#' @examples
#' data(ple4)
#' data(ple4.index)
#'
#' # fit using the default submodels
#' fit1 <- sca(ple4, FLIndices(ple4.index))
#' plot(ple4 + fit1)
#'
#' # see default submodels (set through automated procedure)
#' sca(ple4, FLIndices(ple4.index), fit='assessment')
#'
#' # fishing mortality by age and year (separable)
#' fit2 <- sca(ple4, FLIndices(ple4.index), fmodel=~factor(age) + factor(year))
#' plot(ple4 + fit2)
#' wireframe(data~year*age, data=harvest(fit2), zlab="F")
#'
#' # fit2 + catcability as a smoother by age without year trend
#' fit3 <- sca(ple4, FLIndices(ple4.index), fmodel=~factor(age) + factor(year), qmodel=list(~s(age, k=4)))
#' plot(ple4 + fit3)
#'
#' # fit3 + srmodel as a smoother by year
#' fit4 <- sca(ple4, FLIndices(ple4.index), fmodel=~factor(age) + factor(year), qmodel=list(~s(age, k=4)), srmodel=~s(year, k=45))
#' plot(ple4 + fit4)
#'
#' AIC(fit1, fit2, fit3, fit4)
#' BIC(fit1, fit2, fit3, fit4)

