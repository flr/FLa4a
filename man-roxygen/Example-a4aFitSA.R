#' @examples
#' data(ple4)
#' data(ple4.index)
#' obj <- a4a(stock=ple4, indices=FLIndices(ple4.index), fit="assessment")
#' obj
#' clock(obj)
#' fitSumm(obj)
#' flq <- stock.n(obj)
#' is(flq)
#' flq <- catch.n(obj)
#' is(flq)
#' flq <- harvest(obj)
#' is(flq)
#' flq <- index(obj)
#' is(flq)
#' logLik(obj)
#' is(pars(obj))

