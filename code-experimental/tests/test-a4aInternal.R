context("a4aInternal")

test_that("a4aInternal runs with default settings without crashing",
{
  data(ple4)
  data(ple4.index)  
  FLa4a:::a4aInternal(stock=ple4, indices = FLIndices(ple4.index), fit="MP")
  FLa4a:::a4aInternal(stock=ple4, indices = FLIndices(ple4.index), , fit="assessment")
  FLa4a:::a4aInternal(stock=ple4, indices = FLIndices(ple4.index), , fit="setup")


})
