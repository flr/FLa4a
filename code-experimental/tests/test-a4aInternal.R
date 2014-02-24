context("a4aInternal")

test_that("a4aInternal runs with default settings without crashing",
{
  data(ple4)
  data(ple4.index)  

  a4aInternal(stock=ple4, indices = FLIndices(ple4.index))
})
