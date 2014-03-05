#' @examples
#' data(ple4)
#' data(ple4.index)
#' fmodel <- ~factor(age) + factor(year)
#' qmodel <- list(~factor(age))
#' fit1 <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, stock=ple4, indices=FLIndices(ple4.index))
#' #ToDo
