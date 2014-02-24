#' @examples
#' library(FLa4a)
#' data(ple4)
#' data(ple4.index)
#' fmodel <- ~factor(age) + factor(year)
#' qmodel <- list(~factor(age))
#' fit. <- a4a(stock=ple4, qmodel = qmodel, fmodel=fmodel, indices=FLIndices(ple4.index), fit ="assessment", wkdir="test")
#' flqs <- residuals(fit., ple4, FLIndices(idx=ple4.index))

