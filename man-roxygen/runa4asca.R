#' @examples
#' data(ple4)
#' data(ple4.index)
#' fmodel <- ~factor(age) + factor(year)
#' qmodel <- list(~s(age, k=4) + year)
#' rmodel <- ~s(year, k=20)
#' obj <-  a4aSCA(fmodel=fmodel, qmodel=qmodel, srmodel=rmodel, ple4, FLIndices(ple4.index))

