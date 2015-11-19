#' @examples
#' data(ple4)
#' data(ple4.index)
#'
#' obj <- sca(stock=ple4, indices=FLIndices(ple4.index), fit="assessment")
#' obj
#'
#' slotNames(obj)
#' clock(obj)
#' fitSumm(obj)
#'
#' flq <- stock.n(obj)
#' is(flq)
#' flq <- index(obj)
#' is(flq)
#'
#' logLik(obj)
#' AIC(obj)
#' BIC(obj)
#'
#' is(pars(obj))
