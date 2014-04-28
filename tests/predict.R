library(FLa4a)
data(ple4)
data(ple4.index)
fit1 <-  a4aSCA(ple4, FLIndices(ple4.index))
flqs <- predict(fit1)
sfrac <- mean(range(ple4.index)[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit1))*sfrac
lst <- dimnames(fit1@index[[1]])
lst$x <- stock.n(fit1)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit1)[[1]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[1]]))
all.equal(c(harvest(fit1)), c(flqs$stkmodel$harvest), tolerance=10e-4)
all.equal(c(stock.n(fit1)[,1]), c(flqs$stkmodel$ny1), tolerance=10e-4)
all.equal(c(stock.n(fit1)[1]), c(flqs$stkmodel$rec), tolerance=10e-4)


