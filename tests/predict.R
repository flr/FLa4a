#====================================================================
# tests for predict
#====================================================================
library(FLa4a)
data(ple4)
data(ple4.indices)
data(ple4.index)

#====================================================================
# abundance indices
#====================================================================
fit <-  a4aSCA(ple4, FLIndices(ple4.index), qmodel=list(~s(age, k=4)))
flqs <- predict(fit)

#--------------------------------------------------------------------
# check
#--------------------------------------------------------------------
sfrac <- mean(range(ple4.index)[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[1]])
lst$x <- stock.n(fit)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit)[[1]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[1]]), tolerance=10e-4)
all.equal(c(harvest(fit)), c(flqs$stkmodel$harvest), tolerance=10e-4)
all.equal(c(stock.n(fit)[,1]), c(flqs$stkmodel$ny1), tolerance=10e-4)
all.equal(c(stock.n(fit)[1]), c(flqs$stkmodel$rec), tolerance=10e-4)

#--------------------------------------------------------------------
# N
#--------------------------------------------------------------------
fit <-  simulate(fit, 2)
flqs <- predict(fit)
sum(unlist(lapply(flqs, is, "FLQuants")))==3
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[1]])
lst$x <- stock.n(fit)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit)[[1]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[1]]), tolerance=10e-4)
all.equal(c(harvest(fit)), c(flqs$stkmodel$harvest), tolerance=10e-4)
all.equal(c(stock.n(fit)[,1]), c(flqs$stkmodel$ny1), tolerance=10e-4)
all.equal(c(stock.n(fit)[1]), c(flqs$stkmodel$rec), tolerance=10e-4)

#====================================================================
# biomass indices
#====================================================================
bioidx <- FLIndexBiomass(FLQuant(NA, dimnames=list(age="all", year=range(ple4)["minyear"]:range(ple4)["maxyear"])), name="bioidx")
index(bioidx) <- stock(ple4)*0.001
index(bioidx) <- index(bioidx)*exp(rnorm(index(bioidx), sd=0.1))
range(bioidx)[c("startf","endf")] <- c(0,0)

# fitting the model
fit <- a4aSCA(ple4, FLIndices(bioidx), qmodel=list(~1))
flqs <- predict(fit)

#--------------------------------------------------------------------
# check
#--------------------------------------------------------------------
sfrac <- mean(range(bioidx)[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[1]])
lst$x <- quantSums(stock.n(fit)*exp(-Z)*stock.wt(ple4))
stkn <- do.call("trim", lst)
qhat <- index(fit)[[1]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[1]]), tolerance=10e-4)
all.equal(c(harvest(fit)), c(flqs$stkmodel$harvest), tolerance=10e-4)
all.equal(c(stock.n(fit)[,1]), c(flqs$stkmodel$ny1), tolerance=10e-4)
all.equal(c(stock.n(fit)[1]), c(flqs$stkmodel$rec), tolerance=10e-4)

#--------------------------------------------------------------------
# further checks
#--------------------------------------------------------------------
qmod <- qmodel(pars(fit))[[1]]
rngyear(qmod) <- NA
rngage(qmod) <- NA
obj <- predict(qmod)
all.equal(dimnames(obj), dimnames(FLQuant(quant="age")))

#--------------------------------------------------------------------
# N
#--------------------------------------------------------------------
fit <- simulate(fit, 2, stock=ple4)
flqs <- predict(fit)
sum(unlist(lapply(flqs, is, "FLQuants")))==3
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[1]])
lst$x <- quantSums(stock.n(fit)*exp(-Z)*stock.wt(ple4))
stkn <- do.call("trim", lst)
qhat <- index(fit)[[1]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[1]]), tolerance=10e-4)
all.equal(c(harvest(fit)), c(flqs$stkmodel$harvest), tolerance=10e-4)
all.equal(c(stock.n(fit)[,1]), c(flqs$stkmodel$ny1), tolerance=10e-4)
all.equal(c(stock.n(fit)[1]), c(flqs$stkmodel$rec), tolerance=10e-4)

#====================================================================
# both
#====================================================================

# fitting the model
fit <- a4aSCA(ple4, FLIndices(bioidx, ple4.index), qmodel=list(~1, ~s(age, k=4)))
flqs <- predict(fit)

#--------------------------------------------------------------------
# check
#--------------------------------------------------------------------
sfrac <- mean(range(bioidx)[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[1]])
lst$x <- quantSums(stock.n(fit)*exp(-Z)*stock.wt(ple4))
stkn <- do.call("trim", lst)
qhat <- index(fit)[[1]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[1]]))

sfrac <- mean(range(ple4.index)[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[2]])
lst$x <- stock.n(fit)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit)[[2]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[2]]), tolerance=10e-4)
all.equal(c(harvest(fit)), c(flqs$stkmodel$harvest), tolerance=10e-4)
all.equal(c(stock.n(fit)[,1]), c(flqs$stkmodel$ny1), tolerance=10e-4)
all.equal(c(stock.n(fit)[1]), c(flqs$stkmodel$rec), tolerance=10e-4)

#--------------------------------------------------------------------
# N
#--------------------------------------------------------------------
fit <- simulate(fit, 2, stock=ple4)
flqs <- predict(fit)
sum(unlist(lapply(flqs, is, "FLQuants")))==3

sfrac <- mean(range(bioidx)[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[1]])
lst$x <- quantSums(stock.n(fit)*exp(-Z)*stock.wt(ple4))
stkn <- do.call("trim", lst)
qhat <- index(fit)[[1]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[1]]))

sfrac <- mean(range(ple4.index)[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[2]])
lst$x <- stock.n(fit)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit)[[2]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[2]]), tolerance=10e-4)
all.equal(c(harvest(fit)), c(flqs$stkmodel$harvest), tolerance=10e-4)
all.equal(c(stock.n(fit)[,1]), c(flqs$stkmodel$ny1), tolerance=10e-4)
all.equal(c(stock.n(fit)[1]), c(flqs$stkmodel$rec), tolerance=10e-4)


#====================================================================
# several
#====================================================================

# fitting the model
fit <- a4aSCA(ple4, ple4.indices, fit="assessment", qmodel=list(~s(age, k=4), ~s(age, k=4), ~s(age, k=3)))
flqs <- predict(fit)

#--------------------------------------------------------------------
# check
#--------------------------------------------------------------------

sfrac <- mean(range(ple4.indices[[1]])[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[1]])
lst$x <- stock.n(fit)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit)[[1]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[1]]), tolerance=10e-4)

sfrac <- mean(range(ple4.indices[[2]])[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[2]])
lst$x <- stock.n(fit)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit)[[2]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[2]]), tolerance=10e-4)

sfrac <- mean(range(ple4.indices[[3]])[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[3]])
lst$x <- stock.n(fit)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit)[[3]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[3]]), tolerance=10e-4)
all.equal(c(harvest(fit)), c(flqs$stkmodel$harvest), tolerance=10e-4)
all.equal(c(stock.n(fit)[,1]), c(flqs$stkmodel$ny1), tolerance=10e-4)
all.equal(c(stock.n(fit)[1]), c(flqs$stkmodel$rec), tolerance=10e-4)

#--------------------------------------------------------------------
# N
#--------------------------------------------------------------------
fit <- simulate(fit, 2)
flqs <- predict(fit)
sum(unlist(lapply(flqs, is, "FLQuants")))==3

sfrac <- mean(range(ple4.indices[[1]])[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[1]])
lst$x <- stock.n(fit)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit)[[1]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[1]]))

sfrac <- mean(range(ple4.indices[[2]])[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[2]])
lst$x <- stock.n(fit)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit)[[2]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[2]]), tolerance=10e-4)

sfrac <- mean(range(ple4.indices[[3]])[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[3]])
lst$x <- stock.n(fit)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit)[[3]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[3]]), tolerance=10e-4)
all.equal(c(harvest(fit)), c(flqs$stkmodel$harvest), tolerance=10e-4)
all.equal(c(stock.n(fit)[,1]), c(flqs$stkmodel$ny1), tolerance=10e-4)
all.equal(c(stock.n(fit)[1]), c(flqs$stkmodel$rec), tolerance=10e-4)

#====================================================================
# S/R
#====================================================================
fit <-  a4aSCA(ple4, FLIndices(ple4.index), srmodel=~bevholt(CV=0.1))
flqs <- predict(fit)

#--------------------------------------------------------------------
# check
#--------------------------------------------------------------------
sfrac <- mean(range(ple4.index)[c("startf", "endf")])
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[1]])
lst$x <- stock.n(fit)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit)[[1]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[1]]), tolerance=10e-4)
all.equal(c(harvest(fit)), c(flqs$stkmodel$harvest), tolerance=10e-4)
all.equal(c(stock.n(fit)[,1]), c(flqs$stkmodel$ny1), tolerance=10e-4)
all.equal(c(stock.n(fit)[1]), c(flqs$stkmodel$rec), tolerance=10e-4)

#--------------------------------------------------------------------
# N
#--------------------------------------------------------------------
fit <-  simulate(fit, 2)
flqs <- predict(fit)
sum(unlist(lapply(flqs, is, "FLQuants")))==3
Z <- (m(ple4) + harvest(fit))*sfrac
lst <- dimnames(fit@index[[1]])
lst$x <- stock.n(fit)*exp(-Z)
stkn <- do.call("trim", lst)
qhat <- index(fit)[[1]]/stkn
all.equal(c(qhat), c(flqs$qmodel[[1]]), tolerance=10e-4)
all.equal(c(harvest(fit)), c(flqs$stkmodel$harvest), tolerance=10e-4)
all.equal(c(stock.n(fit)[,1]), c(flqs$stkmodel$ny1), tolerance=10e-4)
all.equal(c(stock.n(fit)[1]), c(flqs$stkmodel$rec), tolerance=10e-4)



